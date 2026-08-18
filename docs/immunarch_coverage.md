# immunarch 功能覆盖审计

本审计以 immunarch 当前官方[函数索引](https://immunarch.com/reference/)及其
[分析教程](https://immunarch.com/articles/)为准。`immune` 不复制 immunarch
代码；目标是在 Scanpy 风格 API 下，用成熟 Python 库覆盖同类科学问题。

状态说明：`已覆盖`表示存在可用的 Python API；`部分覆盖`表示核心能力可用但缺少
immunarch 的某些模式；`未覆盖`表示当前没有对应功能。

| immunarch 功能族 | 状态 | `immune` API / 差距 |
|---|---|---|
| MiXCR、10x、AIRR 数据读取 | 已覆盖 | `iu.io.read_mixcr/read_10x/read_airr`；未覆盖 immunarch 的全部旧格式 |
| 基础统计：克隆数、测序量、长度分布 | 已覆盖 | `iu.tl.bulk_summary`、`iu.tl.spectratype` |
| Shannon、Simpson、Chao1、丰富度 | 已覆盖 | `iu.tl.alpha_diversity`，统计由 scikit-bio 提供 |
| Pielou、Gini、inverse Simpson、Hill numbers | 已覆盖 | `iu.tl.alpha_diversity` 默认指标与 `iu.tl.hill_diversity` 多 q profile，均调用 scikit-bio |
| D50 / DXX | 已覆盖 | `iu.tl.coverage_diversity`、`iu.tl.dxx` 与 `iu.pl.coverage_diversity` |
| 稀释曲线 | 已覆盖 | `iu.tl.rarefaction`，基于 scikit-bio subsampling |
| clonality 标量 | 已覆盖 | `iu.tl.repertoire_metrics` |
| rank/proportion clonality 分箱及逐克隆注释 | 已覆盖 | `iu.tl.clonality_rank/clonality_proportion` 及两个 annotate API |
| V/D/J/C gene usage | 已覆盖 | `iu.tl.segment_usage`，支持加权、非加权和标准化 |
| gene family/allele 层级、歧义基因策略 | 已覆盖 | `iu.tl.gene_usage` 支持 family/segment/allele 与 first/exclude/combine |
| CDR3 spectratype | 已覆盖 | `iu.tl.spectratype`、`iu.pl.spectratype` |
| 样本间 Jaccard/Bray-Curtis 等距离 | 已覆盖 | `iu.tl.beta_diversity`、`iu.pl.repertoire_overlap` |
| shared/public clonotype 数、Jaccard 相似度 | 已覆盖 | `iu.tl.public_overlap` 支持 intersection/Jaccard/overlap coefficient |
| Tversky、cosine、Morisita、incremental overlap | 未覆盖 | 可优先调用 SciPy/scikit-bio，Morisita 需单独验证实现 |
| public repertoire 构建、组合共享统计 | 部分覆盖 | `iu.tl.public_repertoire` 已提供 cohort table；精确样本组合统计待补 |
| clonotype 时间追踪 | 已覆盖 | `iu.tl.timecourse`、`iu.tl.longitudinal_expansion`、轨迹图 |
| 通用过滤、下采样和重采样 | 已覆盖 | `iu.pp.filter_repertoire/downsample_repertoire` |
| 单细胞 paired-chain 与细胞群分析 | 已覆盖且更全面 | Scirpy、Scanpy、MuData、Milo、scvi-tools 适配层 |
| bulk CDR3 序列距离与聚类 | 未覆盖 | 单细胞可用 Scirpy；bulk 尚无 tcrdist3/Levenshtein 聚类入口 |
| k-mer、PFM/PPM/PWM、sequence logo | 未覆盖 | 可接成熟 Python 序列/Logo 软件，不建议手写核心算法 |
| VDJdb、McPAS、PIRD/TBAdb 注释 | 未覆盖 | 需要数据库下载、版本记录、精确/模糊匹配适配层 |
| CDR3 氨基酸理化性质 | 未覆盖 | 可接现有蛋白序列描述符库 |
| BCR germline、SHM、clonal lineage、系统发育树 | 未覆盖 | 应接 Immcantation/Change-O、IgPhyML 或 Dandelion 等成熟工具 |
| PCA/MDS/聚类等 repertoire-level 后分析 | 部分覆盖 | 已有距离矩阵；缺少统一 PCoA/PCA、聚类和 metadata 检验 API |
| 差异丰度 | 已覆盖并增强 | bulk 使用 PyDESeq2；单细胞使用 Pertpy Milo |
| publication-ready 绘图 | 部分覆盖 | bulk、clonality、DXX/Hill、public 和 Scirpy 图已覆盖；缺 circos/logo 图 |

## 建议实现顺序

1. **队列与特征工程**：public 组合共享统计、metadata 自动合并、PCoA/PCA、组间
   PERMANOVA/ANOSIM、ML-ready sample feature table。
2. **成熟软件适配**：tcrdist3、VDJdb/McPAS、k-mer/sequence logo、氨基酸性质。
3. **BCR 专线**：优先连接 Dandelion 或 Immcantation 生态，不在 `immune` 内
   从头实现 germline 重建、SHM 和谱系树。

## 结论

`immune` 当前已覆盖 immunarch 最常用的 bulk 基础分析和本轮核心增强，同时提供
更完整的单细胞 RNA/TCR 联合分析。后续重点转向队列统计、序列相似性/抗原注释
与 BCR 专线，而不是重复实现成熟算法。
