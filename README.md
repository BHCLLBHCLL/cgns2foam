# cgns2foam

基于 `h5py` 的纯 Python 工具，把 CFD 中常用的 **CGNS（HDF5 编码）** 网格文件
转成 **OpenFOAM** 工程目录，开箱即用，无需安装 CGNS C/C++ 库。

> 仅依赖 `h5py` 与 `numpy`。  
> **目标 OpenFOAM 版本：openfoam.com 的 v2412**（ESI/OpenCFD 发行版）。
> 默认以 **binary polyMesh + ANSA 25.1 兼容头** 写出（`location ""`、
> `neighbour` 全长且边界面为 `-1`、ANSA 风格 `note` 与 banner 间距）。
> OpenFOAM v2412 与 ANSA 25.1 均可读取。若只需 OpenFOAM 原生布局，
> 使用 ``--openfoam-native``。

## 目录结构

```
src/                   # 转换器实现（Python 包，导入名 `src`）
├── reader.py          # 基于 h5py 的 CGNS 读取（CPEX-0001 / SIDS-to-HDF5）
├── topology.py        # NGON / NFACE → OpenFOAM polyMesh 拓扑与重排
├── couplings.py       # 扫描流-流 / 流-固 / 固-固耦合界面对
├── regions_config.py  # 同名 JSON 中的区域/材料/热源/重力/并行等配置
├── cht_case.py        # CHT 各文件模板（thermo / fv* / 0 场 / MRF）
├── cht_direct.py      # 一步 CGNS -> 多区域 CHT
├── writer.py          # OpenFOAM polyMesh / system / constant / 0 文件生成
├── convert.py         # 高层流水线（reader → topology → writer [→ CHT]）
└── __main__.py        # `python -m src` 命令行入口
tests/
├── validate.py        # 重新读取写出的二进制网格做对比
├── run_all.py         # 三个 case 的端到端跑通脚本（含 checkMesh）
├── test_box.py        # unittest：拓扑统计与 ANSA 头格式
├── test_bc_overlap.py # unittest：跨 zone 同名 BC 重叠裁剪
├── test_couplings.py  # unittest：区域类型与耦合报告
└── test_regions_config.py # unittest：JSON 区域/材料/热源/重力等配置解析
cases/                 # 测试 case + ANSA 产生的参考 OpenFOAM 工程
requirements.txt
docs/
└── TECHNICAL.md       # 详细技术文档（数据结构 / 算法 / 限制）
```

> 注：项目品牌名仍是 **cgns2foam**（README 标题、文件头 banner、日志
> 前缀等用户可见处保持不变），仅源码所在目录命名为 `src/`。Python 的
> 导入名因此是 `src`，CLI 通过 `python -m src` 调用。

## 安装

```bash
python3 -m pip install -r requirements.txt
```

无需编译 CGNS 原生库，只要环境有 HDF5 和 numpy 即可。

## 快速使用

```bash
# 在仓库根目录下运行（保证 `src` 包能被找到）：
# 基本用法：从 .cgns 文件生成同名 OpenFOAM 工程目录
python3 -m src path/to/case.cgns

# 指定输出目录
python3 -m src path/to/case.cgns /tmp/myCase

# 安静模式（脚本场景）
python3 -m src -q path/to/case.cgns /tmp/myCase

# OpenFOAM 原生 polyMesh（非 ANSA 回导场景）
python3 -m src --openfoam-native path/to/case.cgns /tmp/myCase

# 仅扫描体区域与耦合界面对（不写 polyMesh）
python3 -m src --scan path/to/case.cgns
python3 -m src --scan path/to/case.cgns --report couplings.json

# 一步到位：CGNS -> 多区域 chtMultiRegionSimpleFoam（无需 split）
python3 -m src --cht-direct path/to/case.cgns /tmp/myChtReady
```

### 耦合扫描与 CHT 模式

`--cht-direct` **必须**在 CGNS 旁提供同名 JSON（`mesh.cgns` → `mesh.json`），
用 `fluid_regions` / `solid_regions` 声明流体/固体区域（foam2thermal 风格
的 `regions` 声明也兼容）。

> **注**：旧的 `--cht`（mono + `splitMeshRegions` 两阶段）模式已移除——
> 该路径不会把流-固界面转成 `mappedWall`，生成的脚手架无法直接共轭耦合。
> 请统一使用 `--cht-direct`。

界面方法（按区域类型强制）：

| 耦合 | OpenFOAM 界面 |
|------|----------------|
| fluid–fluid | `cyclicAMI`（同 JSON 区域的多个 cellZone 会合并到同一 polyMesh） |
| fluid–solid / solid–solid | `mappedWall` |

| 开关 | 作用 |
|------|------|
| `--scan` | 读取 CGNS zone / BC，输出流-流、流-固、固-固耦合关系与界面对（有同名 JSON 则用其区域定义） |
| `--cht-direct` | **一步**写出 `constant/<region>/polyMesh`（JSON 区域合并）；流-流 `cyclicAMI`，流-固/固-固 `mappedWall` |
| `--report PATH` | 将扫描结果写为 JSON |
| `--solid-pattern` / `--fluid-pattern` | 无 JSON 时覆盖固体 / 流体 zone 命名规则（可重复） |

同名 JSON 最简格式（参考
`tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json`）：

```json
{
  "fluid_regions": [
    "laptop_3d_geom.air.air_domain",
    "FPHPARTS.rotation1",
    "FPHPARTS.rotation2"
  ],
  "solid_regions": [
    "laptop_3d_geom.fan1.case1",
    "solid_region.Cu_block",
    "solid_region.Cover"
  ]
}
```

每个字符串对应一个 CGNS zone。**所有 `fluid_regions` 放入单一区域
`air`**（`constant/air/polyMesh`，域内流-流界面为 `cyclicAMI`）；
每个固体 zone 各自一个区域（sanitized 名）。不另建名为 `fluid` 的区域。

同名 JSON 还支持以下可选键（完整示例见
`tests/laptop_thermal_steady_scaled_v3_orig_BCs_fix.json`）：

| 键 | 作用 |
|----|------|
| `"g"` / `"gravity"` | 重力矢量，写入 `constant/g`（默认 `(0 0 -9.81)`） |
| `"heat_sources"` | `{"<区域或 zone 名>": <总功率 W>}`；同区域多个键会**求和**；写为 `scalarSemiImplicitSource` + `volumeMode absolute`（OpenFOAM 按单元体积自动摊成 W/m³） |
| `"materials"` | 按区域覆盖物性：固体 `rho`/`Cp`/`kappa`/`molWeight`，流体 `air` 可覆盖 `mu`/`Pr`/`Cp`/`molWeight`（未给的项用内置默认：空气 / 铝） |
| `"external_convection"` | `{"patches": [正则...], "Ta": 300, "h": 8}`；命中的外壁面写 `externalWallHeatFluxTemperature`（`mode coefficient`） |
| `"initial_conditions"` | `{"T": 300.0, "p": 101325.0}` 全局初值 |
| `"n_procs"`（或 `"parallel": {"nProcs": N}`） | MPI 核数，同时写进 `decomposeParDict` 与 `Allrun`（默认 8） |
| `"endTime"` / `"writeInterval"` / `"purgeWrite"` | `controlDict` 覆盖（默认 500 / 50 / 0） |
| `"mrf_regions"` | 见下文 |

边界条件自动生成规则（`0/<region>/`）：

- 耦合界面：mappedWall → `compressible::turbulentTemperatureRadCoupledMixed`（T）/ `noSlip`（U）/ `fixedFluxPressure`（p_rgh）；cyclicAMI → 同名 constraint 类型。
- 开口：patch 名以 `open` 开头，**或** CGNS BC 类型为 `BCInflow`/`BCOutflow`/`BCFarfield` 等 → `prghTotalPressure` + `pressureInletOutletVelocity` + T `inletOutlet`（转换时会打印警告提醒校核数值）。
- `symmetryPlane` / `empty` / `wedge` / `cyclic` 等 constraint patch 自动写同名场类型。
- 其余壁面：U `noSlip`、T `zeroGradient`（绝热），可用 `"external_convection"` 改成对外散热。

求解设置要点（自动生成）：

- `SIMPLE`：`momentumPredictor true`、`nNonOrthogonalCorrectors 1`、显式 `residualControl`（流体 `p_rgh/U/h`，固体 `h`）。
- `div(phi,U/h/K)` 默认一阶 `bounded Gauss upwind` 保证起步稳健；`fvSchemes` 中附有二阶 `linearUpwind` 注释模板，收敛后可切换。

可选 MRF（写出 `constant/air/MRFProperties`，叶轮 patch 用 `movingWallVelocity`）：

```json
"mrf_regions": [
  {
    "cellZone": "FPHPARTS.rotation1",
    "origin": [-0.0678, -0.003, 0.081],
    "axis": [0.0, 1.0, 0.0],
    "omega": 100
  },
  {
    "cellZone": "FPHPARTS.rotation2",
    "origin": [0.0802, -0.003, 0.081],
    "axis": [0.0, -1.0, 0.0],
    "omega": 100
  }
]
```

`origin` 也可为 `"centroid"`（用该 zone 顶点坐标均值）。
MRF 的 `cellZone` 必须是 `fluid_regions` 中的 CGNS zone，否则转换即报错。

```bash
# --cht-direct（一步）
cd /tmp/myChtReady && ./Allrun
```

生成的目录结构为标准 OpenFOAM 工程（v2412 期望布局）：

```
myCase/
├── 0/                          # 初始条件（U、p、p_rgh）
├── constant/
│   ├── polyMesh/
│   │   ├── points              # binary vectorField（默认，ANSA 可回导）
│   │   ├── faces               # binary faceCompactList
│   │   ├── owner               # binary labelList
│   │   ├── neighbour           # binary labelList（长度 = nFaces，边界面 -1）
│   │   ├── boundary            # patch 字典（头 format binary，体 ASCII）
│   │   ├── cellZones           # 每个 CGNS Zone 对应一个 cellZone
│   │   │                       # （class regIOobject，v2412 期望的命名）
│   │   └── faceZones           # 空（与 ANSA 输出保持一致）
│   └── turbulenceProperties    # ANSA 头：location ""、format binary
└── system/
    ├── controlDict             # ANSA 头；writeCompression uncompressed
    ├── fvSchemes               # 最小化占位，v2412 启动时强制要求
    └── fvSolution              # 最小化占位
```

## 在 Python 里调用

```python
# 仓库根目录加入 sys.path 后即可：
from src import convert_file, scan_file, WriteOptions

mesh = convert_file("case.cgns", "out_dir", verbose=True)
print(mesh.n_cells, mesh.owner.size, len(mesh.patches))

# 仅扫描耦合
report = scan_file("case.cgns", report_path="couplings.json")
print(report.fluid_regions, len(report.couplings))

# 一步多区域 CHT
from src import convert_cht_direct
convert_cht_direct("case.cgns", "cht_ready")
# 或: convert_file("case.cgns", "cht_ready", cht_direct=True)

# OpenFOAM 原生 binary 输出
mesh = convert_file(
    "case.cgns", "out_dir",
    write_options=WriteOptions.openfoam_native(),
)
```

## 端到端测试

仓库自带 3 个测试 case：

| Case               | nPoints | nFaces  | nCells | nPatches | 说明                       |
|--------------------|---------|---------|--------|----------|----------------------------|
| `box_ansa`         | 8       | 18      | 6      | 1        | 6-tet 立方体（小烟测试）   |
| `tr03`             | 220 333 | 323 146 | 63 882 | 10       | 2 个 zone 的旋转机械网格   |
| `laptop_simplified`| 735 419 | 1 643 240 | 482 034 | 9      | 3 个 zone 的笔记本散热网格 |

跑全部 case（可选附带 OpenFOAM 的 `checkMesh`）：

```bash
# 仅做拓扑/统计对比
python3 tests/run_all.py

# 一并跑 checkMesh（需要先安装 openfoam2412，来源：https://dl.openfoam.com/）
#   curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
#   sudo apt-get install -y openfoam2412
python3 tests/run_all.py --with-checkmesh
```

测试脚本会按以下顺序探测 OpenFOAM v2412 的 `bashrc`：

```
/usr/lib/openfoam/openfoam2412/etc/bashrc   # apt 包默认路径
/opt/openfoam2412/etc/bashrc                # tarball 解压常见路径
/opt/OpenFOAM-v2412/etc/bashrc              # 手动编译常见路径
$FOAM_BASHRC                                # 环境变量覆盖
```

样例输出节选：

```
### Case: box_ansa
   nPoints / nFaces / nCells / patches 全部与 ANSA 参考完全一致
   checkMesh ours: Mesh OK
   checkMesh ref : Mesh OK

### Case: tr03
   nPoints=220333, nFaces=323146, nCells=63882  ← 全部与 ANSA 参考一致
   checkMesh ours: Failed 3 mesh checks.    (源 CGNS 自身的几何问题)
   checkMesh ref : Failed 3 mesh checks.    (同样的失败计数 → 源数据问题)

### Case: laptop_simplified
   nPoints=735419, nFaces=1643240, nCells=482034  ← 全部一致
   checkMesh ours: Failed 2 mesh checks.    (源 CGNS 自身的几何问题)
   checkMesh ref : Failed 2 mesh checks.
```

## 更多细节

- 转换算法、数据结构、ANSA 写出格式、BC 重叠裁剪、限制和扩展点详见
  [docs/TECHNICAL.md](docs/TECHNICAL.md)（§2.1 写出选项、§3.5 BC 裁剪）。
- 单元测试：
  ```bash
  python3 -m unittest tests.test_box tests.test_bc_overlap tests.test_couplings tests.test_regions_config -v
  ```
