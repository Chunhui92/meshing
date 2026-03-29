# MIR 工程说明

## 概要

当前仓库实现的是一个 **zoo-only 的第一阶段 MIR 原型**：

- 输入：Silo mixed-material 数据
- 输出：VTK `UNSTRUCTURED_GRID`
- 当前主线：共享点、节点体积分数字段、逐材料顺序裁剪、局部体积分数校正、VTK 分析/可视化

这不是完整的 VisIt 风格 Zoo MIR，也不是最终 watertight 的生产级多材料重建器；当前目标是建立一个稳定、可扩展、可验证的 `zoo` 主干。

## 目录结构

- [main.py](/d:/agentsAI/meshing/main.py)
  - 根目录唯一入口，当前只暴露 `zoo` 模式
- [src/mir_mesh.py](/d:/agentsAI/meshing/src/mir_mesh.py)
  - Silo/PDB reader、mesh 检测、zone 迭代
- [src/mir_reconstruct.py](/d:/agentsAI/meshing/src/mir_reconstruct.py)
  - 节点体积分数重采样、共享点缓存、zoo-only 重建、局部体积分数校正
- [src/mir_utils.py](/d:/agentsAI/meshing/src/mir_utils.py)
  - 几何辅助函数
- [test/analyze_vtk.py](/d:/agentsAI/meshing/test/analyze_vtk.py)
  - VTK 可读性、体积、拓扑、阈值化 pass/fail 分析
- [test/visualize_vtk.py](/d:/agentsAI/meshing/test/visualize_vtk.py)
  - VTK 可视化
- [MIR_ROADMAP.md](/d:/agentsAI/meshing/MIR_ROADMAP.md)
  - 当前 zoo-only 路线图

## 当前 `zoo` 实现

### 1. 材料预处理

- 每个 zone 先按 `--min-fraction` 过滤小体积分数材料
- 对剩余材料重新归一化
- 按体积分数从大到小排序

### 2. 几何表示

- zone 仍然使用确定性的 face-consistent tet 辅助分解
- 保留三类点：
  - corner nodes
  - edge intersections
  - zone-local helper points
- 共享点缓存保留，但缓存键已扩展为 `edge + clip context`

### 3. 重建路径

- 二材料 zone 与多材料 zone 共用同一条 `zoo` 路径
- 对当前 residual 逐个材料做顺序裁剪
- 每一步使用节点体积分数字段构造标量场
- 最后一个材料接收 residual
- 不再回退到 `plane` 或 `fast`

### 4. 体积分数校正

- 每次裁剪都做一维二分校正
- 目标是在不改变拓扑框架的前提下，让当前材料片段体积逼近目标体积分数
- 当前暴露参数：
  - `--volume-correction-iters`
  - `--volume-correction-tol`

## 当前实现边界

已实现：

- 轻量 Silo mixed-material 读取
- rectilinear mesh
- UCD 8-node hex mesh
- multiblock UCD 8-node hex mesh（通过 block 级 `mesh1/zl1/mat1` 聚合）
- `cell -> node` 体积分数重采样
- zoo-only 重建
- 共享点缓存
- 局部体积分数校正
- VTK 输出
- VTK 定量分析和可视化

未实现：

- 完整 VisIt 风格 zoo case-table
- 完整 `>2` 材料生产级 watertight 输出保证
- 2D mesh 上的 zoo 重建（例如 `torus_mat.silo` 这类仅可做材料分布检查的数据）
- Isovolume fallback
- species / mixed variable 传播
- 大规模性能优化、并行化、GPU

## 常用命令

重建：

```powershell
& 'D:\programs\conda\python.exe' main.py data\rect3d.silo rect3d_zoo.vtk
& 'D:\programs\conda\python.exe' main.py data\ucd3d.silo ucd3d_zoo.vtk
& 'D:\programs\conda\python.exe' main.py data\rect3d.silo rect3d_smoke.vtk --max-zones 200
```

分析：

```powershell
& 'D:\programs\conda\python.exe' test\analyze_vtk.py rect3d_zoo.vtk --reference-silo data\rect3d.silo --max-boundary-edges 0 --max-rel-total-error 1e-3
& 'D:\programs\conda\python.exe' test\analyze_vtk.py ucd3d_zoo.vtk --reference-silo data\ucd3d.silo --max-boundary-edges 0 --max-rel-total-error 1e-3
```

可视化：

```powershell
& 'D:\programs\conda\python.exe' test\visualize_vtk.py rect3d_zoo.vtk --interfaces
```

材料分布检查：

```powershell
& 'D:\programs\conda\python.exe' test\inspect_silo_materials.py data\multi_ucd3d.silo
& 'D:\programs\conda\python.exe' test\inspect_silo_materials.py data\torus_mat.silo
```

## 已验证结果

### `rect3d.silo` smoke (`--max-zones 200`)

- 8426 cells
- 56482 points
- `boundary_edges = 0`
- 最大总体积相对误差约 `2.58e-06`
- `non_manifold_edges` 仍非零

### `ucd3d.silo` smoke (`--max-zones 200`)

- 3060 cells
- 7426 points
- `boundary_edges = 0`
- 最大总体积相对误差约 `6.72e-06`
- `non_manifold_edges` 仍非零

### 人工 3 材料 mock case

- 3 材料顺序裁剪路径可执行
- 能生成可读 VTK
- 当前仓库样例数据本身只包含 1/2 材料 zone，因此 `>2` 材料回归仍需要补专门测试集

### 新补数据检查

- `multi_ucd3d.silo`
  - 可由当前 multiblock loader 读取
  - 全量 `53176` 个 zone 中，`3808` 个是二材料 zone
  - 不包含 `>2` 材料 zone
- `torus_mat.silo`
  - 当前 3D zoo 重建器不支持其 mesh 类型
  - 但材料块统计显示：`840` 个单材料、`16` 个二材料、`44` 个三材料 zone
  - 可作为后续 `>2` 材料材料分布检查样例，但不能直接用于当前 3D 重建验证

## 工程判断

当前版本的定位是：

**一个可运行、可分析、可扩展、zoo-only 的 MIR 第一阶段原型。**

它已经具备继续扩展所需的主干结构，但仍需要继续降低 non-manifold 问题、补多材料专门回归样例，并逐步接近更完整的 Zoo MIR。
