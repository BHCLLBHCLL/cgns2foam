# cgns2foam

基于 `h5py` 的纯 Python 工具，把 CFD 中常用的 **CGNS（HDF5 编码）** 网格文件
转成 **OpenFOAM** 工程目录，开箱即用，无需安装 CGNS C/C++ 库。

> 仅依赖 `h5py` 与 `numpy`。  
> **目标 OpenFOAM 版本：openfoam.com 的 v2412**（ESI/OpenCFD 发行版）。
> 输出按 OpenFOAM 默认编译参数（32 位 label、双精度）的二进制格式生成，
> 并与 v2412 自带的 `checkMesh` 工具完全对接（已验证）。

## 目录结构

```
src/                   # 转换器实现（Python 包，导入名 `src`）
├── reader.py          # 基于 h5py 的 CGNS 读取（CPEX-0001 / SIDS-to-HDF5）
├── topology.py        # NGON / NFACE → OpenFOAM polyMesh 拓扑与重排
├── writer.py          # OpenFOAM polyMesh / system / constant / 0 文件生成
├── convert.py         # 高层流水线（reader → topology → writer）
└── __main__.py        # `python -m src` 命令行入口
tests/
├── validate.py        # 重新读取写出的二进制网格做对比
├── run_all.py         # 三个 case 的端到端跑通脚本（含 checkMesh）
└── test_box.py        # unittest 单元测试
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
```

生成的目录结构为标准 OpenFOAM 工程（v2412 期望布局）：

```
myCase/
├── 0/                          # 初始条件（U、p、p_rgh）
├── constant/
│   ├── polyMesh/
│   │   ├── points              # 二进制 vectorField
│   │   ├── faces               # 二进制 faceCompactList
│   │   ├── owner               # 二进制 labelList
│   │   ├── neighbour           # 二进制 labelList（长度 = nInternalFaces）
│   │   ├── boundary            # ASCII patch 字典
│   │   ├── cellZones           # 每个 CGNS Zone 对应一个 cellZone
│   │   │                       # （class regIOobject，v2412 期望的命名）
│   │   └── faceZones           # 空（与 ANSA 输出保持一致）
│   └── turbulenceProperties
└── system/
    ├── controlDict
    ├── fvSchemes               # 最小化占位，v2412 启动时强制要求
    └── fvSolution              # 最小化占位
```

## 在 Python 里调用

```python
# 仓库根目录加入 sys.path 后即可：
from src import convert_file

mesh = convert_file("case.cgns", "out_dir", verbose=True)
print(mesh.n_cells, mesh.owner.size, len(mesh.patches))
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

- 转换算法、数据结构、限制和扩展点详见 [docs/TECHNICAL.md](docs/TECHNICAL.md)。
- 单元测试：`python3 -m unittest tests.test_box -v`
