# cgns2foam 技术文档

本文档详述 **cgns2foam** 项目从 CGNS（HDF5）到 OpenFOAM `polyMesh` 的
完整转换过程：数据模型、关键算法、文件二进制布局、已知限制以及扩展
指引。源代码统一放在仓库根目录下的 `src/` 包里（Python 导入名为
`src`，CLI 通过 `python -m src` 调用），项目品牌名仍为 cgns2foam。

**目标 OpenFOAM 发行版**：openfoam.com 的 **v2412**（ESI/OpenCFD）。
我们没有对 openfoam.org 的 Foundation 版本（11/12/13 等）做兼容；
两者在以下几个地方有不同的 polyMesh I/O 期望（详见 §2 / §4）：

| 文件 / 配置项                       | openfoam.com v2412               | openfoam.org 13                  |
|-------------------------------------|----------------------------------|----------------------------------|
| `cellZones` 的 `class` 头           | `regIOobject`                    | `cellZoneList`                   |
| `faceZones` 的 `class` 头           | `regIOobject`                    | `faceZoneList`                   |
| `controlDict` 的 `writeCompression` | `on`/`off`                       | `compressed`/`uncompressed`      |
| `system/fvSchemes`/`fvSolution`     | **必需**（即使只跑 `checkMesh`） | 启动 `checkMesh` 时可以缺失      |

我们的 writer 完全按 v2412 的期望产出。如果以后需要兼容 Foundation
版本，把上述 3 处头部 / 关键字改一改即可。

**目标 OpenFOAM 发行版**：openfoam.com 的 **v2412**（ESI/OpenCFD）。
我们没有对 openfoam.org 的 Foundation 版本（11/12/13 等）做兼容；
两者在以下几个地方有不同的 polyMesh I/O 期望（详见 §2 / §4）：

| 文件 / 配置项                       | openfoam.com v2412               | openfoam.org 13                  |
|-------------------------------------|----------------------------------|----------------------------------|
| `cellZones` 的 `class` 头           | `regIOobject`                    | `cellZoneList`                   |
| `faceZones` 的 `class` 头           | `regIOobject`                    | `faceZoneList`                   |
| `controlDict` 的 `writeCompression` | `on`/`off`                       | `compressed`/`uncompressed`      |
| `system/fvSchemes`/`fvSolution`     | **必需**（即使只跑 `checkMesh`） | 启动 `checkMesh` 时可以缺失      |

我们的 writer 完全按 v2412 的期望产出。如果以后需要兼容 Foundation
版本，把上述 3 处头部 / 关键字改一改即可。

---

## 1. CGNS（HDF5）数据模型简述

CGNS 的 HDF5 编码遵循 CPEX 0001（SIDS-to-HDF5 mapping）。其核心规则：

| 概念        | HDF5 表示                                                   |
|-------------|-------------------------------------------------------------|
| CGNS 节点   | HDF5 group                                                  |
| 节点 *label*| 属性 `label`（如 `Zone_t`、`Elements_t`）                   |
| 节点 *name* | HDF5 link 名（group 在父节点里的名字）                      |
| 节点数据    | 子 dataset，名字固定为 `" data"`（**前面带一个空格**）       |
| 数据类型    | 属性 `type`（`MT`/`I4`/`R8`/`C1` …）                         |

也就是说，要读取某个 CGNS 节点的实际数据，写的是
`group[" data"][()]`，而不是 `group["data"][()]`。`cgns2foam`
在 `src/reader.py` 里封装了这一约定。

### 1.1 我们关心的节点类型

| Label                | 作用                                                         |
|----------------------|--------------------------------------------------------------|
| `CGNSBase_t`         | 顶层，载有 `CellDim`、`PhysDim`                              |
| `Zone_t`             | 网格区域（一个 zone = 一组共享拓扑/坐标的点+单元）            |
| `ZoneType_t`         | `"Unstructured"` / `"Structured"`（我们只支持 Unstructured） |
| `GridCoordinates_t`  | 子节点 `CoordinateX/Y/Z`（双精度）                           |
| `Elements_t`         | 一段元素，附带 `etype`、`ElementRange`、`ElementStartOffset`、`ElementConnectivity` |
| `ZoneBC_t` / `BC_t`  | 边界条件，含 `BCType`、`GridLocation`、`PointList`/`PointRange` |
| `FlowSolution_t`     | 流场结果（暂只用来标识网格，不写入 OpenFOAM 0 目录）          |

### 1.2 NGON_n（22）/ NFACE_n（23）

`cgns2foam` 主要处理 **多面体** 网格：

- **NGON_n** 段：列出 *面*，每个面由若干顶点构成；用 `ElementStartOffset`
  指明每张面的起始偏移，用 `ElementConnectivity` 存顶点 id（**1 基**）。
- **NFACE_n** 段：列出 *单元*，每个单元由若干面构成；连接性中的整数是
  **带符号** 的面 id，正号表示该面外法线指向单元外部，负号表示指向内部。

固定形状的单元段（`TETRA_4`、`PYRA_5`、`PENTA_6`、`HEXA_8`）也能被
读取并按 CGNS SIDS §11.2 的局部面定义自动展开成等价的 NGON/NFACE，详见
`src/topology.py::_ngon_from_fixed`。

---

## 2. OpenFOAM polyMesh 二进制文件布局

OpenFOAM 默认编译 `WM_LABEL_SIZE=32`、`WM_PRECISION_OPTION=DP`，所以：

- *label* = `int32` little-endian
- *scalar* = `float64` little-endian
- *vector* = 3 个连续的 `float64`

每个二进制 polyMesh 文件均由 **ASCII 文件头 + ASCII 计数 + 二进制裸数据**
组成（小端、连续打包，无任何填充）：

```
points        :  <header>\n N \n ( <N*3*8 bytes float64> ) \n
faces*        :  <header>\n (N+1) \n ( <(N+1)*4 bytes int32 offsets> )
                 <S> \n ( <S*4 bytes int32 connectivity> ) \n
owner         :  <header>\n N \n ( <N*4 bytes int32> ) \n
neighbour     :  <header>\n M \n ( <M*4 bytes int32> ) \n     (M = nInternalFaces)
```

`*` 处 `faces` 文件用的是 OpenFOAM 的 *CompactList<labelList>* 格式：
先写 `nFaces+1` 个 offsets，再写所有面顶点的扁平连接性。这与 ANSA 的输出格式
一致；OpenFOAM 自身的标准 `faceList` 也兼容此布局，由文件头 `class faceCompactList`
表明。

`boundary` 文件是 ASCII 字典；`cellZones` 是 ASCII 字典外加每个 zone 内嵌的
二进制 `List<label>`（同样是 `N\n(<N*4 bytes int32>)\n`）。

---

## 3. 转换算法

### 3.1 单个 zone 的拓扑构造（`_build_zone_topology`）

输入：NGON 与 NFACE 段。

1. **面数组**：直接复用 NGON 的 `ElementStartOffset` 和
   `ElementConnectivity`，把顶点 id 由 1 基转 0 基。
2. **owner / neighbour / flip 计算**（向量化、O(总面引用数)）：

   ```python
   abs_face  = |NFACE_conn| - first_face_id   # 0 基面 id
   sign      = sign(NFACE_conn)               # ±1
   cell_for_ref = repeat(cellId, NFACE.counts)
   # 按 face id 排序，所有引用同一面的 (cell, sign) 相邻
   order = argsort(abs_face, kind='stable')
   ```

   再用 `searchsorted` 定位每张面的第一条引用，按引用条数（1 或 2）分两类：

   - **边界面**（1 条引用）：owner = 该 cell；若 sign 为负则将面顶点反向，
     使法线指向单元外部（OpenFOAM 要求边界面法线由 owner 指向外部）。
   - **内部面**（2 条引用）：CGNS 中标 `+` 的 cell 我们称为 `pos_cell`，
     标 `-` 的 cell 为 `neg_cell`。OpenFOAM 约定 owner 取下标较小者，
     于是：
     - 若 `pos_cell < neg_cell`：owner = pos_cell, neighbour = neg_cell；
       面方向恰好由 owner 指向 neighbour，**无需翻转**。
     - 反之：owner = neg_cell, neighbour = pos_cell；面方向是从 neighbour
       指向 owner，需要 **翻转** 顶点顺序。

3. **边界条件**：仅采纳 `GridLocation == "FaceCenter"` 的 BC。把其
   `PointList` 中的 1 基面 id 转 0 基，过滤掉非边界面 id（即误把内部面
   写进 BC 的情况）。

### 3.2 多 zone 合并（`src/topology.py::build_mesh`）

cgns2foam 采用 **不做几何 stitching** 的合并策略：

- 顶点 / 单元 / 面分别拼接，分别加上跨 zone 偏移。
- **每个 CGNS zone → 一个 cellZone**，方便后续按子域操作（如 MRF）。
- 若同名 BC 在多个 zone 中出现，依次给后出现者加 `_1`、`_2` 后缀以避免
  patch 名冲突；这与 ANSA 在 `laptop_simplified` 中给 `impeller_1` /
  `impeller_2` 改名为 `impeller_11` / `impeller_21` 的做法等价。
- 同一 zone 内若同一面被多个 BC 引用（ANSA 在旋转机械里常这样标），
  按 BC 出现顺序 **先到先得**，避免一张面被分到两个 patch。
- 任何没被 BC 覆盖的边界面归入 `default_exterior` patch（类型 `wall`），
  ANSA 也使用相同策略。

> **限制**：跨 zone 真正重合的内部界面不会被自动识别并合并为内部面。
> 若有需要，请在 OpenFOAM 中用 `mergeMeshes` + `stitchMesh`，或在
> ANSA 中先做 mesh assembly 再导出 CGNS。我们的转换保证每个 zone 内部
> 的拓扑严格正确。

### 3.3 OpenFOAM 排序约束

OpenFOAM 要求：

1. 内部面排在边界面之前；
2. 内部面按 `(owner, neighbour)` 字典序升序（upper-triangular ordering）；
3. 边界面按 patch 连续聚集。

`build_mesh` 末段用一次稳定排序 + 一次 `argsort` 完成全部重排，并且
同步重排 `face_vertices`、`owner`、`neighbour` 与 `flip`，重新构建
偏移数组。

### 3.4 CGNS BC 类型到 OpenFOAM patch 类型的映射

| CGNS BCType                  | OpenFOAM patch type |
|------------------------------|---------------------|
| `BCWall` / 含 `wall`         | `wall`              |
| `BCSymmetryPlane` / 含 `symmetry` | `symmetryPlane`|
| 含 `axis`                    | `empty`             |
| 其它                          | `patch`             |

映射逻辑在 `src/topology.py::_bc_type_to_foam`，可按需扩展。

---

## 4. 生成的 OpenFOAM 文件汇总

| 路径                                  | 内容                                              |
|---------------------------------------|---------------------------------------------------|
| `constant/polyMesh/points`            | 顶点坐标（二进制 vectorField）                    |
| `constant/polyMesh/faces`             | 面顶点列表（二进制 faceCompactList）              |
| `constant/polyMesh/owner`             | 每个面所属 owner cell id（二进制 labelList）       |
| `constant/polyMesh/neighbour`         | 每个内部面的 neighbour cell id（二进制 labelList）|
| `constant/polyMesh/boundary`          | patch 字典（ASCII）                                |
| `constant/polyMesh/cellZones`         | 每个 CGNS zone 对应一个 cellZone（ASCII + 内嵌二进制 label list，`class regIOobject`）|
| `constant/polyMesh/faceZones`         | 空（`class regIOobject`，与 ANSA 输出保持一致）   |
| `constant/turbulenceProperties`       | `simulationType RAS;` 默认 `laminar`              |
| `system/controlDict`                  | `application UserSolver;` 占位、`writeCompression off;`（用户改成实际求解器）|
| `system/fvSchemes`                    | 最小占位（v2412 启动时必需）                       |
| `system/fvSolution`                   | 最小占位（v2412 启动时必需）                       |
| `0/U`                                 | 体积矢量场，默认 `(0 0 0)`、墙面 `fixedValue`     |
| `0/p`, `0/p_rgh`                      | 体积标量场，默认 `0`、墙面 `zeroGradient`         |

初始条件主要起“撑起 case 骨架”的作用，便于用户随后用编辑器或脚本
（如 `setFields`）覆盖实际场值。

---

## 5. 与 ANSA 25.1 输出的对比（OpenFOAM v2412）

`tests/run_all.py` 会把 cgns2foam 的输出和 `cases/*/*.zip` 里 ANSA
产生的参考工程做以下比较（拓扑不变量）：

- `nPoints` / `nFaces` / `nInternalFaces` / `nCells`
- 边界面总数 / 各 patch 面数之和

`--with-checkmesh` 会调用 OpenFOAM v2412 的 `checkMesh`。由于 ANSA 的
参考压缩包没有 `system/fvSchemes` / `fvSolution`（v2412 要求两者），
测试脚本会把我们生成的两个最小占位文件复制到参考目录之后再跑
`checkMesh`，确保两边在同一条件下比较。

实测结果（OpenFOAM v2412 patch 260127，2026-05）：

| Case               | 拓扑不变量 | `checkMesh` 我们 | `checkMesh` ANSA |
|--------------------|-----------|------------------|------------------|
| `box_ansa`         | 全部一致  | `Mesh OK`        | `Mesh OK`        |
| `tr03`             | 全部一致  | `Failed 3 mesh checks.` | `Failed 3 mesh checks.` |
| `laptop_simplified`| 全部一致  | `Failed 2 mesh checks.` | `Failed 2 mesh checks.` |

`tr03` / `laptop_simplified` 上 `checkMesh` 报错的 9 个面取向错误、
22 个高度非正交面、若干高扭曲面，**在 ANSA 的参考工程里以完全相同
的数量出现**，说明这些是 CGNS 源网格自身的几何问题，而非转换器引入。

差异点（不影响仿真）：

- ANSA 会用 `*_Moving` / `*_Static` 后缀手工拆分接口；本转换器统一
  使用 `_1` / `_2` 后缀做名字去重。
- ANSA 会把附在 `PSHELL` 卡上的面单独抽成 `Default_PSHELL_Property`
  patch；本转换器把它们一并归入 `default_exterior`。

如需匹配 ANSA 的命名风格，可在 `topology.py` 里扩展自定义 BC → patch
名映射。

---

## 6. 局限性

1. **不做几何 stitching**：跨 zone 几何重合的界面不会自动并成内部面；
   多区域情况下 `checkMesh` 会报告“多个不连通的网格区域”。
2. **不支持 Structured zone**：只支持 `Unstructured` 类型的 zone，结构
   网格请先用 cgnsTools 中的 `convertCGNS` 或类似工具转成 NGON/NFACE。
3. **只处理单个 `CGNSBase_t`**：多 base 的文件会主动抛错。
4. **不写 FlowSolution → 时间步**：CGNS 中的 `FlowSolution_t` 仅被识别
   而不被搬运到 OpenFOAM 的时间目录。如需迁移流场，可基于 `reader.py`
   再加 50 行代码完成 cell-centered 标量/矢量场的写入。
5. **只用 32-bit label**：超过 ≈21 亿单元/面/点的网格需要切换到
   `WM_LABEL_SIZE=64` 编译的 OpenFOAM；改写代价是把 `writer.py` 里
   `int32` 改为 `int64`。
6. **仅针对 openfoam.com v2412**：openfoam.org 11/12/13 等 Foundation
   版本对 `cellZones` / `faceZones` 头部、`writeCompression` 关键字、
   以及是否需要 `fvSchemes`/`fvSolution` 期望不同；如需兼容请按 §0
   的对照表调整 `writer.py` 中的相应字段。

---

## 7. 扩展点速查

- **加新的 BC 类型映射**：编辑 `src/topology.py::_bc_type_to_foam`。
- **改 patch 命名/合并策略**：编辑 `src/topology.py::build_mesh` 内的
  去重 / 默认 patch 段。
- **写入流场初值**：在 `src/writer.py` 里仿照 `_write_initial_field`
  实现 `nonuniform List<scalar>` 二进制写入，并在 `src/convert.py` 串
  起来。
- **支持固定形状元素的二维（面）段**：扩展
  `src/topology.py::_FIXED_CELL_FACES` 以及 `_ngon_from_fixed`。

---

## 8. 命令行 / Python API 速查

```bash
# 在仓库根目录下运行，确保 `src` 包能被 import
python3 -m src <in.cgns> [out_dir] [-q|--quiet]
python3 tests/run_all.py [--with-checkmesh] [--out-root /tmp/out]
python3 -m unittest tests.test_box -v
```

```python
from src import read_cgns, convert_file
case = read_cgns("in.cgns")          # 仅读取，得到 CGNSCase
mesh = convert_file("in.cgns", "out")
# mesh.points / mesh.face_offsets / mesh.face_vertices /
# mesh.owner / mesh.neighbour / mesh.n_internal_faces /
# mesh.n_cells / mesh.patches / mesh.cell_zones
```
