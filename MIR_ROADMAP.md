# MIR 后续开发路线图

## 目标

在当前第一阶段原型基础上，把 MIR 逐步升级为：

- 支持一般多材料混合 zone
- 拓扑更稳定
- 体积误差更低
- 验证更自动化
- 大规模数据下更可控

当前已具备：

- `cell -> node` 重采样
- shared-point 缓存
- `zoo / plane / fast` 三条路径
- `src/` 与 `test/` 分层
- VTK 分析和可视化脚本

## 开发顺序

建议严格按下面顺序推进：

1. 完善当前二材料 `zoo` 内核
2. 支持 `>2` 材料 zone
3. 加入体积分数校正
4. 升级回归验证
5. 再考虑性能和大规模数据

不要先做并行/GPU，再回头补算法正确性。

## Phase 1: 完善二材料 `zoo`

目标：

- 降低 `non_manifold_edges`
- 减少退化单元
- 减少局部拓扑伪影

建议工作：

- 把当前 tet-only 分割升级为更完整的 case-table / zoo 分解
- 明确区分 corner nodes、face intersections、内部辅助点
- 补 edge interpolation 的退化保护和重复点合并规则

验收：

- `boundary_edges` 保持为 `0`
- `non_manifold_edges` 显著下降
- 体积误差不明显恶化

## Phase 2: 支持 `>2` 材料 zone

目标：

- 不再把 `>2` 活跃材料 zone 直接回退到 `plane`

建议工作：

- 在 `src/mir_reconstruct.py` 中加入逐材料递归裁剪
- 为中间 shape 保存顶点 ids、材料场值、shape 类型
- 先支持 3 材料和 4 材料情况
- 用 `top-k` 控制复杂度

验收：

- 3/4 材料人工样例可稳定输出
- 不依赖 `plane` fallback
- 输出 VTK 可读取、可视化、可分析

## Phase 3: 体积分数校正

目标：

- 降低当前固定 `phi = 0` 带来的总体积误差

建议工作：

- 加入阻尼式迭代校正
- 校正目标放在体积分数/阈值层，而不是直接改拓扑
- 每轮统计总体积误差和可选的局部误差

验收：

- 总体积误差进一步下降
- 不明显增加拓扑问题和退化单元

## Phase 4: 回归测试升级

目标：

- 把分析脚本升级成可判定 pass/fail 的回归测试工具

建议工作：

- 在 `test/analyze_vtk.py` 中加入阈值参数
- 固定 smoke/regression 数据集
- 让脚本在阈值超标时返回失败状态

验收：

- 能稳定复现当前样例结果
- 可作为后续改动的回归闸门

## Phase 5: 性能与大规模数据

目标：

- 在算法正确性基本稳定后，再处理时间和内存问题

建议工作：

- 评估 chunked writer
- 评估 edge cache 内存增长
- 考虑两阶段输出：计数后分配再写入
- 若必要，再考虑 CPU 并行、prefix-sum 输出、GPU 路线

当前不建议优先做：

- 在 `zoo` 仍有明显 non-manifold 问题时直接做 GPU 化
- 在没有 `>2` 材料支持前大量投入性能优化

## 建议测试集

- 当前样例：
  - `data/rect3d.silo`
  - `data/ucd3d.silo`
- 小规模 smoke case：
  - `--max-zones` 子集
- 人工构造 case：
  - 二材料平面切分
  - 同心层
  - 三材料交汇
  - 薄层/小夹杂

## 当前最值得做的一步

优先把当前二材料 `zoo` 路径从“共享点 + tet-only 近似分割”升级为更完整的 case-table / zoo 分解。

原因：

- 它同时决定拓扑质量、`>2` 材料扩展方式、体积误差上限，以及是否值得继续做性能优化。
