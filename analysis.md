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
  
  
  我已重新通读 `src/cht_case.py`、`src/cht_direct.py`、`src/regions_config.py`、JSON 配置以及实际生成的算例文件（`0/`、`constant/`、`system/`）。下面是对**当前代码/算例**的重新分析。

---

## 一、Claude 修正后已确认正确的部分 ✅

| 项 | 文件 | 状态 |
|---|---|---|
| `open` patch 类型 = `patch` | `constant/air/polyMesh/boundary:135` | ✅ |
| `open` 的 `p_rgh` = `prghTotalPressure` (p0=101325) | `0/air/p_rgh:87-92` | ✅ |
| `open` 的 `U` = `pressureInletOutletVelocity` | `0/air/U:80-84` | ✅ |
| `open` 的 `T` = `inletOutlet` (inletValue=300) | `0/air/T:101-106` | ✅ |
| `open` 的 `p` = `calculated` | `0/air/p:87-91` | ✅ |
| `#includeEtc "caseDicts/setConstraintTypes"` 已加到所有 0/ 场 | 各 0/ 文件 line 23 | ✅ |
| CPU 热源 `scalarSemiImplicitSource` + `volumeMode absolute` + `explicit 1.0`（=1W 总功率） | `system/laptop_3d_geom_solid_region_CPU/fvOptions` | ✅ |
| Cu_block 无 fvOptions（不发热） | `system/solid_region_Cu_block/` 无 fvOptions | ✅ |
| CPU 物性已覆盖为硅：ρ=2330, Cp=705, κ=148 | `constant/laptop_3d_geom_solid_region_CPU/thermophysicalProperties` | ✅ |
| `pRefCell`/`pRefValue` 已从 fvSolution 移除 | `system/air/fvSolution` | ✅ |
| 8 核并行 `decomposeParDict` + `Allrun` | `system/decomposeParDict`、`Allrun` | ✅ |
| 重力 g=(0,-9.81,0) 与 JSON 一致 | `constant/g` | ✅ |
| CHT 耦合 `turbulentTemperatureRadCoupledMixed` 双向 | `0/air/T`、`0/.../T` | ✅ |
| MRF `nonRotatingPatches` 含 `open`、不含 impeller | `constant/air/MRFProperties` | ✅ |

这些上一轮分析中的关键问题都已落地。

---

## 二、仍然存在的关键问题（直接导致 Zmin ~3.5 m/s）

### 🔴 问题 1：整个外边界只有**一个** `open` patch，且是 `prghTotalPressure` —— 这是 Zmin 高速的根因

`constant/air/polyMesh/boundary` 显示 `open` 只有 **1 个 patch、17746 个面**，即空气域的全部外边界（很可能是包围盒的 6 个面合并）都被当作同一个总压出口。

后果：
- `prghTotalPressure` 是"总压给定、方向自适应"的边界 —— 同一个 patch 上**一面在出流、另一面在入流**是合法的。
- 两个 MRF 风扇把内部空气搅动，气流撞到 Zmin（盒子底部）那部分 `open` 面就**出流**；为了质量守恒，同一 `open` patch 的其它面（顶部/侧面）就**入流**，入流速度被 `pressureInletOutletVelocity` 钳成 `(0 0 0)`。
- 出流区与入流区在同一 patch 上紧挨（尤其在 Zmin 边缘/角点），形成强烈剪切层 → **角点附近 U ≈ 3.5 m/s**，这正是你观察到的现象。
- 这不是网格问题，也不是离散格式问题，是**边界拓扑问题**：把"应当作为墙面/明确入口/明确出口"的多个面合并成了一个自适应总压面。

### 🔴 问题 2：没有明确的"入口"，算例本质是"开盒+总压边界+内部风扇"

当前唯一能让空气进出的就是这一个 `open`。风扇只是内部再循环驱动者。物理上等价于"笔记本浸在无限大静止空气里"，但边界又离风扇太近（Zmin 紧贴风扇出口），所以近场速度被边界吸入/排出主导，CHT 温升被边界来的 300K 冷风持续冲刷掉。

### 🟡 问题 3：`pressureInletOutletVelocity` 的 `value uniform (0 0 0)` 在入流区强制零速

入流区被钉死为 0 速度，与出流区的有限速度在同一 patch 上对接 → 局部加速度大、`limitU max=100` 不会触发但剪切层速度可达数 m/s。这是问题 1 的放大器。

### 🟡 问题 4：`inletOutlet` T=300K 把入流钉成环境温度

与问题 1 叠加：只要 `open` 同时有入流和出流，入流区不断灌入 300K 冷气，1W 热源的温升会被立刻吹散，CHT 看不到效果（与你上一轮"2000 步 T 几乎不变"一致）。

---

## 三、次要问题（不致错但影响精度/收敛）

### 🟢 问题 5：`div(phi,*)` 全部一阶 upwind
`system/air/fvSchemes:30-32` 仍是 `bounded Gauss upwind`。一阶迎风数值粘性大，会抹平速度梯度、让 Zmin 剪切层看起来"扩散"得更广。注释里已写了切换到 `linearUpwind` 的方法，但默认未启用。稳态收敛后应切换。

### 🟢 问题 6：`rho` 松弛因子 0.05 偏低
`system/air/fvSolution:48` `rho 0.05`。对 perfectGas + 强制对流，ρ 变化不大时 0.05 可接受，但会拖慢收敛。可提到 0.5–1.0。

### 🟢 问题 7：固体 `relaxationFactors equations { h 1; }` 无欠松弛
`system/.../fvSolution`。固体纯导热线性问题，h=1 没问题；但 CHT 耦合初期固体温度会震荡，可临时用 0.7 加速稳定后再回到 1。

### 🟢 问题 8：`limitT max=500` 偏低
1W 热源下不会触发，但若后续加大功率（如 20W CPU），硅芯片局部温度可能超过 500K 被截断。建议提到 600。

### 🟢 问题 9：MRF omega=100 rad/s（≈955 RPM）
对笔记本风扇偏低（实际 3000–5000 RPM）。若是有意测试可忽略；若想看真实流场，应按实际转速设到 300–500 rad/s。这不是 Zmin 高速的原因（即便 ω=100，叶尖速度 ~5 m/s 量级，与观察一致）。

---

## 四、根因小结与建议

**Zmin ~3.5 m/s 的根因**：`open` 这一个 patch 同时承担了"入口+出口"两种角色，且物理上把本应是**墙面**或**分开设定的入口/出口**的外边界全部合并成了一个自适应总压面。风扇驱动的内部环流在 Zmin 出流、在其它面入流，同一 patch 上的出入流剪切层产生了你看到的 3.5 m/s。

**建议的修正方向**（按优先级）：

1. **拆分 `open` patch**（最关键）：在 CGNS 端或后处理阶段，把外边界按几何面拆成 6 个独立 patch（如 `open_xmin / open_xmax / open_ymin / open_ymax / open_zmin / open_zmax`），然后：
   - 把**真正应作为入口**的面（如风扇对面、顶部）设为 `pressureInletOutletVelocity` + `inletOutlet T=300` + `prghTotalPressure`；
   - 把**应作为出口**的面（如风扇出风方向那侧）设为 `inletOutlet` + `pressureInletOutletVelocity` + `prghTotalPressure`；
   - 把**应是墙面**的面（如笔记本底面 Zmin 若是桌面侧）改回 `wall` + `noSlip` + `fixedFluxPressure` + `zeroGradient T`。
   这样出入流分离，剪切层消失。

2. **若几何上 Zmin 确实是开口**（不是桌面）：保持 `open` 但把 `pressureInletOutlet2. **若几何上 Zmin 确实是开口**（不是桌面）：保持 `open` 但把 `pressureInletOutletVelocity` 的 `value` 从 `(0 0 0)` 改成 `internalField` 或方向化初值，减少入流区强制零速造成的剪切。同时把 `inletOutlet` 的 `inletValue` 改成 `zeroGradient`（让出流区不被 300K 钳制），仅保留真正入口面为 300K。

3. **稳态收敛后**：把 `div(phi,U/h/K)` 从 `upwind` 切到 `bounded Gauss linearUpwind grad(...)`，去掉一阶数值粘性。

4. **加大热源到真实值**（如 20W）后再看 CHT 效果；1W 太小，温升会被边界入流完全淹没。

---

需要我直接帮你改代码吗？如果要，我建议先确认一件事：**Zmin 那个面在物理上到底是"桌面（应是 wall）"还是"开口（应是 open）"**？这决定了是走方向 1（拆 patch + 改 wall）还是方向 2（保留 open 但调 BC）。


  修正内容（已提交于 4dc6efb，README 最后两行测试清单微调待提交）

  P0（功能性 bug）

  ┌───────────────────────────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┬───────────────────────────────┐
  │             问题              │                                                                   修法                                                                   │             验证              │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───────────────────────────────┤
  │ 热源单位错（W 被当 W/m³）     │ _fv_options_solid_heat 改 volumeMode absolute，同区域多热源求和                                                                          │ 新案例 CPU fvOptions 已是     │
  │                               │                                                                                                                                          │ absolute                      │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───────────────────────────────┤
  │ JSON "g" 被忽略               │ regions_config 解析 "g"/"gravity" → 透传至 constant/g                                                                                    │ 新案例 value (0.0 -9.81 0.0)  │
  │                               │                                                                                                                                          │ ✓                             │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───────────────────────────────┤
  │ residualControl { default }   │ 流体写 p_rgh/U/h、固体写 h 显式场名                                                                                                      │ 已核验                        │
  │ 无效                          │                                                                                                                                          │                               │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───────────────────────────────┤
  │ 进/出口静默变壁面             │ CGNS bc_type 随 patch plan 透传（Patch.cgns_bc_type），BCInflow/BCOutflow/BCFarfield 等自动按总压开口生成，并在控制台 +                  │ open 口 U/p_rgh/T 正确        │
  │                               │ setup_report.json 的 warnings 中提示                                                                                                     │                               │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───────────────────────────────┤
  │ constraint patch 启动崩溃     │ _constraint_rule 让 symmetryPlane/empty/wedge/cyclic 场类型与网格 patch 自动一致                                                         │ 单测覆盖                      │
  └───────────────────────────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┴───────────────────────────────┘

  P1

  - materials 段：固体按区域覆盖 rho/Cp/kappa/molWeight，流体 air 覆盖 mu/Pr/Cp。新案例 Cu=8960/385/390、CPU=2330/705/148、Cover=塑料，不再全是铝。
  - external_convection 段：正则命中 patch → externalWallHeatFluxTemperature（mode coefficient）。
  - 弃用 --cht：CLI、convert_file(cht=...)、write_cht_case/_allrun_pre_full/_create_patch_ami 死代码全部移除；文档注明移除原因。

  P2

  - SIMPLE：nNonOrthogonalCorrectors 1 + momentumPredictor true（按你要求）。
  - n_procs/endTime/writeInterval/purgeWrite/initial_conditions 全部 JSON 化，decomposeParDict×3 与 Allrun -np 单一来源。
  - fvSchemes 附二阶 linearUpwind 注释模板；MRF cellZone 校验提前报错；region_type 兜底进 warnings。

  刷新图档与测试

  - cases/laptop_thermal_steady_scaled_v3_orig_BCs_fix_cht_air 已用新代码+新 JSON 重新生成（约 4 分钟），上述每项均在生成文件中逐一核验；setup_report.json 新增
  settings/materials_applied/opening_patches/warnings 元数据。
  - 单元测试 24/24 通过（新增 7 个 TestExtendedConfig 用例；test_box 3 个失败为预先存在的 LFS 数据缺失，与本次无关）。
  - README / AGENTS / TECHNICAL 已同步（--cht 移除说明、新 JSON 键表、BC 生成规则、求解设置说明）。

  可直接 cd cases/laptop_thermal_steady_scaled_v3_orig_BCs_fix_cht_air && ./Allrun 跑新案例——这次应能看到温度真正上升且 SIMPLE 按 residualControl 提前收敛。