# MIR 工程说明

## 概要

当前仓库实现的是一个 **第一阶段 MIR 原型**：

- 输入：Silo mixed-material 数据
- 输出：VTK `UNSTRUCTURED_GRID`
- 重点能力：节点重采样、共享点缓存、二材料 `zoo` 主路径、VTK 分析/可视化

这不是完整的 VisIt 风格 Zoo MIR。当前实现仍然只对二材料混合 zone 做强支持。

## 目录结构

- [main.py](/d:/agentsAI/meshing/main.py)
  - 根目录唯一入口
- [mir_mesh.py](/d:/agentsAI/meshing/src/mir_mesh.py)
  - Silo/PDB reader、mesh 检测、zone 迭代
- [mir_reconstruct.py](/d:/agentsAI/meshing/src/mir_reconstruct.py)
  - 节点重采样、共享点缓存、三种重建模式
- [mir_utils.py](/d:/agentsAI/meshing/src/mir_utils.py)
  - 几何辅助函数
- [analyze_vtk.py](/d:/agentsAI/meshing/test/analyze_vtk.py)
  - VTK 可读性、体积、拓扑定量分析
- [visualize_vtk.py](/d:/agentsAI/meshing/test/visualize_vtk.py)
  - VTK 可视化
- [MIR_ROADMAP.md](/d:/agentsAI/meshing/MIR_ROADMAP.md)
  - 后续开发路线图

## 当前支持的模式

### `--mode zoo`

当前默认模式。

做法：

- 先做 `cell -> node` 体积分数重采样
- 用共享点缓存复用边交点
- 对二材料 mixed zone 构造面一致的 12-tet 辅助分解
- 在 tet 顶点上用 `phi = VF(mat_a) - VF(mat_b)` 做局部分割

特点：

- 比 `fast` 平滑
- 比旧的局部独立切分更一致
- 仍有 non-manifold 表面问题

### `--mode plane`

旧基线路径。

做法：

- 用全局材料质心估计 split normal
- 在 tetrahedralized zone 上按体积分数做平面裁剪

特点：

- 体积误差较小
- 在 `ucd3d.silo` 上拓扑更干净
- 在 `rect3d.silo` 上非常慢

### `--mode fast`

rectilinear-only 速度基线。

做法：

- 仅沿主轴切分 mixed voxel

特点：

- 很快
- 体积几乎精确
- 界面明显方块化

## 当前实现边界

已实现：

- 轻量 Silo mixed-material 读取
- rectilinear mesh
- UCD 8-node hex mesh
- `cell -> node` 重采样
- shared-point cache
- `zoo / plane / fast`
- in-memory VTK writer
- VTK 定量分析和可视化

未实现：

- `>2` 材料完整 zoo 重建
- 完整 zoo case-table
- Isovolume 兜底
- 体积分数迭代校正
- species / mixed variable 传播
- 完整拓扑稳定的 watertight 多材料界面

## 常用命令

重建：

```powershell
& 'D:\programs\conda\python.exe' main.py data\rect3d.silo rect3d_zoo.vtk --mode zoo
& 'D:\programs\conda\python.exe' main.py data\rect3d.silo rect3d_fast_new.vtk --mode fast
& 'D:\programs\conda\python.exe' main.py data\ucd3d.silo ucd3d_zoo.vtk --mode zoo
& 'D:\programs\conda\python.exe' main.py data\ucd3d.silo ucd3d_plane_new.vtk --mode plane
```

分析：

```powershell
& 'D:\programs\conda\python.exe' test\analyze_vtk.py rect3d_zoo.vtk --reference-silo data\rect3d.silo
& 'D:\programs\conda\python.exe' test\analyze_vtk.py ucd3d_zoo.vtk --reference-silo data\ucd3d.silo
```

可视化：

```powershell
& 'D:\programs\conda\python.exe' test\visualize_vtk.py rect3d_zoo.vtk --interfaces
& 'D:\programs\conda\python.exe' test\visualize_vtk.py rect3d_fast_new.vtk --interfaces
```

## 已验证结果

### `rect3d.silo`

`--mode zoo`：

- 553101 cells
- 124376 points
- 最大总域相对体积误差约 `1.93e-3`
- `boundary_edges = 0`
- `non_manifold_edges` 非零

`--mode fast`：

- 49818 cells
- 260489 points
- 体积几乎精确
- 外观明显方块化

`--mode plane`：

- 当前环境下全量运行超过 4 小时仍未完成

### `ucd3d.silo`

`--mode zoo`：

- 20500 cells
- 5845 points
- material 1 / 4 总域相对体积误差约 `2.49e-3`
- `boundary_edges = 0`
- `non_manifold_edges` 非零

`--mode plane`：

- 31420 cells
- 121846 points
- 体积误差接近数值零
- `boundary_edges = 0`
- `non_manifold_edges = 0`

## 工程判断

当前版本的定位是：

**一个可运行、可分析、可扩展的 MIR 第一阶段原型。**

它已经具备继续扩展所需的基本结构，但算法核心还没有达到完整、稳健、一般多材料的 Zoo MIR 水平。

后续工作请直接参考 [MIR_ROADMAP.md](/d:/agentsAI/meshing/MIR_ROADMAP.md)。
