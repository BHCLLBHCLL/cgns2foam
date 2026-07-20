● 分析完成。我把源码（src/cht_case.py、src/cht_direct.py、src/writer.py、src/regions_config.py、src/topology.py）与生成
  的算例 cases/laptop_thermal_steady_scaled_v3_orig_BCs_fix_cht_air/
  及其求解日志做了交叉核对。以下按五个类别列出仍然存在的问题，标注 🔴 已确认的 bug / 🟡 物理建模缺口 / ⚪
  健壮性或易用性问题。

  ---
  1. 边界条件

  🔴 1.1 JSON 里的 "g" 被静默丢弃，浮升力方向错误

  - JSON 中写了 "g": [0, -9.81, 0]，但 regions_config.py 全文没有解析 "g" 键；write_cht_direct_case 虽有 gravity
  形参（cht_direct.py:477），convert_cht_direct 调用时（cht_direct.py:681）从不传它。
  - 结果：constant/g 永远落回默认值 (0 0 -9.81)。本案例风扇轴为 ±Y，重力应为
  −Y，实际算例的浮升力方向是错的，且无任何警告。

  🔴 1.2 热源单位错误：总功率 W 被当成 W/m³

  - JSON heat_sources 语义是"总功率（瓦）"，但 _fv_options_solid_heat（cht_case.py:366）写的是：
  volumeMode  specific;      ← OpenFOAM 语义为 W/m³
  explicit    1;
  - specific 模式下 explicit 1 = 1 W/m³。CPU 区域体积 ~1e-6 m³ 量级，实际注入功率 ~1e-6 W。日志佐证：跑满 500
  步后全场最高温只有 300.0031 K，热源几乎等于零。
  - 正确做法：volumeMode absolute（OpenFOAM
  按选中单元总体积归一），或转换器自行累加区域体积后换算功率密度。cht_case.py:369-371 的 docstring
  对此的解释本身就是错的。

  🔴 1.3 除 open* 以外的进/出口被静默改成绝热壁面

  - 场文件生成只按 patch 名字 做启发式判断（cht_case.py:485-596）：不是 ami_*、不含 _to_、不以 open 开头、不含 impeller
  → U 给 noSlip、T 给 zeroGradient。
  - CGNS 的 BCInflow/BCOutflow/BCFarfield 等类型信息在 _bc_type_to_foam（topology.py:388）里统一映射为 patch
  后就丢失了——一个叫 inlet 的进口会变成 patch 网格类型 + noSlip 速度，求解照常收敛，物理全错，零警告。
  - CGNS BC_t 下的 BCDataSet（Dirichlet/Neumann 数据：进口速度、流量、壁温、热流密度）reader
  完全不读，源算例携带的物理边界值全部丢弃。

  🔴 1.4 constraint patch（symmetryPlane / empty / wedge）在 CHT 0/ 场中类型不匹配，启动即 fatal

  - _field_U/_field_T/_field_p/_field_p_rgh 没有 symmetryPlane/empty 分支，一律落到 noSlip/zeroGradient/calculated。
  - 而 polyMesh 中这些 patch 类型是 symmetryPlane/empty，OpenFOAM 要求场类型与 constraint patch 类型一致——带对称面或 2D
  empty 的网格在 cht-direct 模式下启动即报类型不匹配错误。单块模式的
  _write_initial_field（writer.py:562-581）反而处理了这两种类型，CHT 路径是倒退。

  🟡 1.5 --cht（两阶段）模式的区域界面不会热耦合

  - mono 网格里流-固界面 patch 类型是 wall；createPatchDict 只转换同区域的 cyclicAMI
  对（cht_case.py:884-918），没有任何步骤把界面转成 mappedWall；splitMeshRegions 后界面仍是 wall，0.orig/T 给
  zeroGradient → 各区域热学上完全解耦，CHT 名存实亡。
  - 0.orig 场对 AMI patch 的判断还要求 CGNS BC 名里含 "ami"（_is_cyclic_ami_patch），不匹配时 createPatch 后 patch
  类型（cyclicAMI）与场类型（noSlip）冲突，启动报错。

  🟡 1.6 外壁面只能绝热

  机壳/盖板外表面一律 zeroGradient，无法配置 externalWallHeatFluxTemperature（对外自然对流+辐射）、固定热流或固定壁温。
  对笔记本散热这类"对外散热才是主通路"的场景，生成的 case 必须手工改 BC 才有物理意义。

  ⚪ 1.7 其他

  - MRF 全靠名字猜：叶轮识别靠 "impeller" in name（cht_direct.py:553），nonRotatingPatches 按 ami_*/open*/*_to_*
  猜（cht_case.py:421），旋转轴默认值按 zone 名含 "rotation1/rotation2"
  猜（regions_config.py:324）。命名一变就静默失效（叶片变固定壁面，风扇不转但照常收敛）。
  - 同一 OpenFOAM 区域内的固-固 STITCH 界面被跳过、落为普通 wall（cht_direct.py:176-178，代码注释自承 "leave as ordinary
  walls for now"）→ 区域内隔热。
  - matchTolerance 0.001、sampleMode nearestPatchFace 硬编码；default_exterior 兜底面默认
  wall，几何缝隙造成的未覆盖面会静默变成绝热壁。
  - open 口 T 的 inletValue（回流温度）硬编码 300 K。

  ---
  2. 初始条件

  🔴 2.1 初值全部硬编码，JSON 无入口

  - T0=300.0（cht_case.py:498）、p0=101325.0（cht_case.py:618）、U=(0 0 0)；_field_* 虽有 T0/p0
  形参，write_cht_direct_case 从不传值。无法按区域给定不同初温（如预热固体）。

  🟡 2.2 p_rgh 未做静水压初始化

  p 与 p_rgh 都是均匀 101325，没有按 p_rgh = p − ρ·g·h 初始化；叠加上 1.1
  的重力方向错误，初始场与体积力项双重不一致（稳态 SIMPLE 能消化，但白白增加前几十步残差）。

  🟡 2.3 湍流场整体缺失

  0/air 只有 T/U/p/p_rgh，没有 k/omega/epsilon/nut/alphat——与 laminar 自洽，但用户一旦切换 RAS
  模型，连字段文件都得纯手工补。

  ⚪ 2.4 不搬 CGNS FlowSolution

  源文件中的流场/温度场无法作为初值（文档 §6.4 已承认），无法热启动。

  ---
  3. 物性参数

  🔴 3.1 所有固体共用一套铝物性

  _DEFAULT_SOLID_MAT（cht_case.py:70-86：ρ=2719, Cp=871,
  κ=202.4，铝）被写进每一个固体区域（cht_direct.py:560）。证据：生成的
  constant/solid_region_Cu_block/thermophysicalProperties 就是铝参数——铜块（应 ρ≈8960,
  κ≈390）、CPU、塑料风扇壳、盖板全部变成铝。JSON 没有 materials 段，无法按区域区分。

  🟡 3.2 流体物性固定为 300 K 空气

  hConst + const 输运（μ=1.846e-5, Pr=0.706，cht_case.py:53-68），不可配；没有 Sutherland/多项式 Cp
  选项；合并流体区域名固定为 "air"，与工质解耦。

  🟡 3.3 湍流模型固定 laminar

  风扇强制对流 + 浮升力混合对流在物理上是湍流；turbulenceProperties 写死 simulationType
  laminar（cht_case.py:175），无模型选择入口。单块模式写的 simulationType RAS; RAS { RASModel laminar;
  }（writer.py:531）也是个四不像。

  🟡 3.4 热源能力受限

  - 每个固体区域只能有一个热源：next(...) 取第一个匹配（cht_direct.py:566-571），同区域第二个热源被静默丢弃；
  - 只能 selectionMode all，不能按 cellZone 局部施加；不能给流体区域加热源；
  - heat_sources 的 key 靠子串模糊匹配（regions_config.py:471-476），有误配风险。

  🟡 3.5 辐射固定关闭

  radiationModel none（cht_case.py:397）。自然对流笔记本外壳辐射占比可观，无配置入口。

  ---
  4. 离散格式

  🟡 4.1 对流项全一阶迎风

  div(phi,U/h/K) 全部 bounded Gauss upwind（cht_case.py:226-233），精度低、数值耗散大；没有"一阶启动 → 二阶 linearUpwind
  收尾"的分档模板或注释指引。

  🟡 4.2 非正交处理保守但不闭环

  - laplacianSchemes/snGradSchemes 用 limited 0.333（固体还是 0.33，不统一）——对 README 自己记录的高非正交网格（laptop
  案例 checkMesh fail 2 项）是稳妥选择，但 0.333 对温度场耗散不小；
  - 配套的非正交修正却是 nNonOrthogonalCorrectors 0（见 5.2），耗散与误差两头都占了。
  - gradSchemes 只有裸 Gauss linear，扭曲网格上没有 cellLimited Gauss linear 1 之类的限制器。

  ⚪ 4.3 单块模式的 fvSchemes 是不可运行的占位

  divSchemes { default none; }（writer.py:494）只够过 checkMesh；application
  UserSolver;（writer.py:453）同理。文档虽有说明，但"开箱即用"的宣传与"必须手改才能跑"的实际有落差。

  ---
  5. 求解参数

  🔴 5.1 residualControl { default 1e-5; } 是无效配置

  OpenFOAM 的 residualControl 按场名匹配，不存在名为 default 的场 → 没有任何收敛判据生效 → 永远跑满
  endTime。日志佐证：500 步全部跑完，全程无 "converged" 字样。流体、固体的 fvSolution
  都犯了同样的错（cht_case.py:284、cht_case.py:313）。应写显式场名（p_rgh/U/h）分级给容差。

  🟡 5.2 nNonOrthogonalCorrectors 0

  本网格 checkMesh 报高度非正交面（README 实测记录），p_rgh 方程无非正交修正有压力场误差风险；建议至少 1–2。固体
  fvSolution 同样没有。

  ⚪ 5.3 并行核数硬编码 8，且分散在三处

  system/decomposeParDict、system/<region>/decomposeParDict（cht_direct.py:504/545/564）、Allrun 的 -np
  8（cht_direct.py:460）各自硬编码；改一处漏两处就崩，也不按机器核数自适应。

  ⚪ 5.4 controlDict 参数不可达

  write_cht_direct_case(end_time=500, write_interval=50) 有形参（cht_direct.py:476），但
  convert_cht_direct/convert_file/CLI 整条链路都不传值（convert.py:44-51、__main__.py 无对应 flag）——永远是
  500/50。purgeWrite 0 导致稳态 10 份中间结果全留盘。

  ⚪ 5.5 其他

  - functions {} 为空：无 solverInfo/residuals/minMax 监控，只能肉眼翻 log；
  - 欠松弛因子写死（p_rgh 0.3、rho 0.05、U/h 0.3）；limitT（200–500 K）、limitU（100）硬编码；
  - 顶层 system/fvSchemes、system/fvSolution 写的是流体版，多区域求解器根本不读，属于误导性冗余；
  - region_type.get(foam_name, "fluid")（cht_direct.py:507）把未识别区域静默当流体；
  - resolve_mrf_entries 不校验 MRF cellZone 是否真的属于 air 区域网格，指错 zone 时到 OpenFOAM 运行期才报错。

  ---
  优先级建议

  ┌────────┬────────────────────────────┬───────────────────────────────────────────────────────────────────────────┐
  │ 优先级 │            问题            │                                   理由                                    │

  OpenFOAM 的 residualControl 按场名匹配，不存在名为 default 的场 → 没有任何收敛判据生效 → 永远跑满 endTime。日志佐证：500 步全部跑完，全程无 "converged" 字样。流体、固体的 fvSolution
  都犯了同样的错（cht_case.py:284、cht_case.py:313）。应写显式场名（p_rgh/U/h）分级给容差。

  🟡 5.2 nNonOrthogonalCorrectors 0

  本网格 checkMesh 报高度非正交面（README 实测记录），p_rgh 方程无非正交修正有压力场误差风险；建议至少 1–2。固体 fvSolution 同样没有。

  ⚪ 5.3 并行核数硬编码 8，且分散在三处

  system/decomposeParDict、system/<region>/decomposeParDict（cht_direct.py:504/545/564）、Allrun 的 -np 8（cht_direct.py:460）各自硬编码；改一处漏两处就崩，也不按机器核数自适应。

  ⚪ 5.4 controlDict 参数不可达

  write_cht_direct_case(end_time=500, write_interval=50) 有形参（cht_direct.py:476），但 convert_cht_direct/convert_file/CLI 整条链路都不传值（convert.py:44-51、__main__.py 无对应 flag）——永远是
  500/50。purgeWrite 0 导致稳态 10 份中间结果全留盘。

  ⚪ 5.5 其他

  - functions {} 为空：无 solverInfo/residuals/minMax 监控，只能肉眼翻 log；
  - 欠松弛因子写死（p_rgh 0.3、rho 0.05、U/h 0.3）；limitT（200–500 K）、limitU（100）硬编码；
  - 顶层 system/fvSchemes、system/fvSolution 写的是流体版，多区域求解器根本不读，属于误导性冗余；
  - region_type.get(foam_name, "fluid")（cht_direct.py:507）把未识别区域静默当流体；
  - resolve_mrf_entries 不校验 MRF cellZone 是否真的属于 air 区域网格，指错 zone 时到 OpenFOAM 运行期才报错。

  ---
  优先级建议

  ┌────────┬────────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────────────┐
  │ 优先级 │                              问题                              │                                           理由                                            │
  ├────────┼────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
  │ P0     │ 1.2 热源单位、1.1 重力被忽略、5.1 residualControl 无效         │ 三个都是"算完了但物理不对/没收敛判据"，且测试案例的日志已实证（T≤300.003 K、跑满 500 步） │
  ├────────┼────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
  │ P0     │ 1.3 进/出口静默变壁面、1.4 constraint patch 启动崩溃           │ 静默错误比报错更危险                                                                      │
  ├────────┼────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
  │ P1     │ 3.1 全固体铝物性（JSON 加 materials 段）、1.6 外壁面对流散热   │ 决定结果是否有物理意义                                                                    │
  ├────────┼────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
  │ P1     │ 1.5 --cht 两阶段模式界面不耦合                                 │ 要么补 mappedWall 转换，要么在文档/CLI 明确标记弃用、只推 --cht-direct                    │
  ├────────┼────────────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────┤
  │ P2     │ 5.2 非正交修正、4.1 二阶格式分档、5.3 核数参数化、2.1 初值可配 │ 收敛性与精度改进                                                                          │
  └────────┴────────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────────────┘

  如需，我可以先修 P0 的三项（热源改 volumeMode absolute、JSON 解析 "g" 并透传、residualControl 写显式场名），都是小改动且能立刻在现有案例日志里验证效果。