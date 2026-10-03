# immune：定位、架构与接口设计

状态：单细胞 RNA/VDJ 生态的基础版本已实现，2026 年 10 月。

本文是 [英文设计文档](DESIGN.md) 的中文版本，说明包的功能定位、数据模型、模块、接口、统计解释和扩展方向。接口名称与代码保持一致；尚未验证或尚未实现的能力会单独说明。

## 1. 包的定位

immune 分析 bulk 与单细胞适应性免疫受体组库。当前重点是联合研究单细胞转录组和 VDJ：识别克隆，描述克隆的细胞表型，追踪不同时间的检出与扩增，并将 MiXCR bulk 克隆关联到单细胞状态。

它包括 PhenoTrack 风格的克隆/状态分布与独立近似流量分析，以及 CloneTrack 风格的纵向扩增分析；保留 Python 原生 bulk 统计和 immunarch 的重要功能。包内不提供 VDJtools 命令行互操作模块。

## 2. 与其他包的职责边界

| 包/后端 | 职责 |
| --- | --- |
| cellscope | RNA 质控、归一化、降维、聚类、注释、功能评分、pseudobulk 与表达推断 |
| immune | 受体链、克隆、免疫组库统计、克隆/表型关系、纵向追踪、bulk 匹配与克隆分组后的联合分析 |
| Scirpy | AIRR 读取、链索引与质控、精确克隆型、受体相似性聚类及受体可视化 |

immune 从原生数据对象读取 cellscope 的结果，并构造受体定义的细胞分组。The application calls cellscope directly for RNA analysis; neither package depends on the other.

Scirpy 已有的成熟算法通过其公共接口调用，不重新实现底层算法。Generic RNA adapters have been removed from immune; use cellscope directly.

## 3. 通用数据模型

### 3.1 单细胞：MuData + AnnData

| 存放位置 | 内容与归属 |
| --- | --- |
| `data.mod["gex"]` | RNA AnnData，由 cellscope/Scanpy 分析 |
| `gex.layers["counts"]` | 原始非负整数计数 |
| `gex.obs` | RNA 质控、`cell_type`、`cell_state`、评分 |
| `gex.obsm` | 表达降维结果 |
| `data.mod["airr"]` | 兼容 Scirpy 的 AIRR AnnData |
| `airr.obsm["airr"]` | 可变长的完整受体链记录 |
| `airr.obsm["chain_indices"]` | Scirpy 选择分析链的索引 |
| `airr.obs` | 链质控、克隆标识、克隆大小与扩增状态 |
| `data.obs` | 公共样本/供体/文库/条件/时间信息、带模态前缀的注释、`has_gex` 和 `has_airr` |
| `data.uns["immune"]` | 联合结果表、参数、后端来源信息 |

不另造替代 AnnData/MuData 的单细胞容器。对象保存使用原生 `.h5ad` 和 `.h5mu`。

### 3.2 bulk 与研究管理

已有 `ImmuneProject` 保留研究管理及表格流程：

| 字段 | 内容 |
| --- | --- |
| `bulk_chains` | 标准化 bulk 链表 |
| `single_cell_chains` | 标准化的旧式单细胞链表 |
| `sample_metadata` | 样本层面的研究信息 |
| `links` | 显式的组库关联关系 |

bulk 和单细胞保持不同的观察轴。bulk 行表示测量到的受体重排，read/template/UMI 计数不能当成虚构的细胞。

已有标准化表格仍是有效输入。`io.read_10x` 返回原有表格格式；`io.read_10x_vdj` 返回原生 AIRR AnnData，两者不混用。

### 3.3 细胞身份与模态对齐

RNA 和 AIRR 都使用 `library_id:barcode` 作为细胞标识，保留原始 `barcode`、生物学 `sample_id` 和 `donor_id`。成对 RNA/VDJ 数据的文库标识必须相同。不同样本可重复使用原始 10x barcode，但不能共享组合后的细胞标识。

默认外连接保留 RNA-only 和 VDJ-only 细胞。联合表型分析只使用具有所需注释的细胞，并记录分母。内连接需要显式指定 `join="inner"`。

同一细胞的公共元数据冲突会报错。更新注释时显式发布 `gex:cell_state`、`airr:clone_id` 等带前缀的全局列，兼容 MuData ≥0.4 的行为。

## 4. 受体与克隆的定义

Scirpy 的 AIRR/Awkward 表示保留完整链记录。链索引选择 productive 分析链，但不会销毁其他原始记录。

精确克隆型使用核苷酸序列 identity。受体相似性分组通过 Scirpy 的独立接口 `define_clonotype_clusters` 完成，不能把相似性簇与严格克隆型混为一谈。

`tl.define_clonotypes` 默认在同一供体内，使用配对受体臂和全部双链规则进行匹配，因此默认要求非缺失的供体注释。研究跨供体序列共享时，必须有意识地设置 `scope=None`。

Scirpy 数值型簇 ID 只在当前分析对象内有效。应将对象与克隆定义一起保存，不能直接把不同独立分析中的数字 ID 当成同一个克隆。

序列类型、距离度量、供体范围、受体臂、双链规则及基因约束记录在 `uns["immune"][key]` 中。

完整配对克隆身份与 bulk 单链匹配键是不同概念。一条 TRB 可以对应多个具有不同 TRA 的配对克隆；匹配结果保留所有关联及歧义，不擅自简化成一对一。

## 5. 公共模块与接口

```python
import immune as iu
```

| 模块 | 职责 |
| --- | --- |
| `iu.io` | 标准化/原生组库读取与保存 |
| `iu.pp` | 数据约定、过滤、链质控、对象组装 |
| `iu.tl` | 组库、表型、追踪、bulk 匹配与联合表达分析 |
| `iu.pl` | 受体与联合分析绘图 |
| `iu.get` | 结果、链记录与注释提取 |
| `iu.datasets` | 可复现的配对供体示例 |

下列表格中的 `...` 表示省略了部分参数，不是可以直接执行的完整调用。

### 5.1 输入输出与预处理

| 接口 | 行为 |
| --- | --- |
| `io.read_mixcr`、`io.read_airr`、`io.read_10x` | 已有标准化 bulk/单细胞链读取器 |
| `io.read_10x_vdj(path, library_id=..., sample_id=..., donor_id=...)` | 调用 Scirpy 读取 CSV/JSON，保留原始链 |
| `io.read_airr_anndata(path, library_id=..., ...)` | 原生 AIRR 读取，使用组合细胞标识 |
| `io.read_h5ad`、`io.read_h5mu`、`io.write` | 原生对象持久化 |
| `pp.to_scirpy(chains)` | 将已有标准化链表转成 AIRR AnnData |
| `pp.merge_with_transcriptome(gex, airr, join="outer")` | 检查一致性并组装 MuData |
| `pp.scirpy_qc(data, ...)` | Scirpy 链索引与链质控 |
| `pp.CloneDefinition(...)` | 不可变的单链匹配键定义，用于 bulk 关联 |
| `pp.clone_abundance`、`pp.filter_repertoire`、`pp.downsample_repertoire` | 表格丰度聚合、过滤与下采样 |

### 5.2 单细胞与联合分析

| 接口 | 输出或用途 |
| --- | --- |
| `tl.define_clonotypes(data, sequence="nt", metric="identity", scope="donor_id", ...)` | 原生精确配对克隆型，或显式指定的相似性簇 |
| `tl.clone_summary(data, sample_col="sample_id", clone_key="clone_id")` | 每个样本/克隆的细胞数量与频率 |
| `tl.clonal_expansion(data, min_cells=2)` | 样本内克隆大小及扩增状态注释 |
| `tl.phenotype_composition(data, phenotype_col="cell_state")` | 克隆/状态计数、克隆内比例与样本比例 |
| `tl.phenotype_diversity(data)` | 状态丰富度、Shannon 熵与有效状态数 |
| `tl.clone_state_enrichment(data)` | 样本内探索性 Fisher 富集检验及 BH 校正 |
| `tl.phenotypic_flux(data, from_sample=..., to_sample=...)` | 共享克隆状态分布的 L1 变化 |
| `tl.phenotype_flow(data, from_sample=..., to_sample=...)` | 独立近似下的克隆状态流出/流入 |
| `tl.track_clones(data, sample_metadata=...)` | 按供体分开的克隆计数、频率与检出轨迹 |
| `tl.longitudinal_expansion(data, sample_metadata=...)` | 供体内基线/随访扩增检验，时间选择校正与全局 BH |
| `tl.match_bulk(data, bulk_chains, sample_pairs=..., definition=...)` | 细胞/完整克隆/单链/bulk 关联表，保留计数单位与歧义 |
| `tl.annotate_bulk_matches(data)` | 对唯一选定的比较结果生成逐细胞 bulk 注释 |

原生表型结果保存在 `data.uns["immune"][key]["table"]`，`params` 记录分母及分析设置。已有表格型表型分析输入和纵向检验仍保留兼容性。原生辅助接口默认使用 `gex` 和 `airr`，必要时可显式指定模态。

RNA expression analysis uses cellscope directly on receptor-labeled observations.
immune contains no expression aggregation or inference wrappers.

### 5.3 bulk 免疫组库统计

现有接口覆盖：

- alpha/beta 多样性、Hill 多样性曲线。
- 覆盖度与 DXX：达到指定累计丰度比例所需的克隆数。
- 克隆性比例/排名、rank-abundance。
- 基因片段和 V/J 使用频率、CDR3 长度谱。
- 公共克隆/组库重叠、稀释曲线。
- 基于 PyDESeq2 的样本层面克隆差异丰度。

成熟多样性方法由 scikit-bio 提供。比较不同流程或实验时，需要明确 `CloneDefinition` 与计数单位，不能直接混合 reads、templates、UMIs 和细胞数。

### 5.4 绘图、提取与示例

`pl` 保留 Scirpy 的受体网络、扩增、基因使用、组库重叠等适配接口，以及已有 bulk Matplotlib 绘图。

新增接口：

- `pl.phenotype_composition`：克隆/状态组成条形图。
- `pl.phenotype_flow`：状态流量热图。
- `pl.clone_embedding`：在 RNA 降维图上叠加受体注释，不永久改变 RNA 注释。

新接口返回原生绘图对象，文件保存需要显式执行。

| 提取接口 | 用途 |
| --- | --- |
| `get.chain_table` | 展开完整 AIRR 链记录，保留细胞元数据 |
| `get.airr` | 调用 Scirpy，访问选定的索引链字段 |
| `get.obs_df` | 对齐 RNA 与受体注释 |
| `get.result` | 返回已保存结果表的副本 |

`datasets.toy_multimodal()` 包含两个供体、前后两个时间样本、配对受体、共享 TRB 但 TRA 不同的克隆，以及 RNA-only 细胞。`datasets.toy_bulk()` 提供对应 bulk 重排。

可执行的完整流程见 [RNA–VDJ–bulk 示例](examples/rna_vdj_workflow.py)。

## 6. 统计解释与使用限制

1. 单细胞克隆大小按细胞计数，不按受体链数或 UMI 数计数。
2. 表型组成比例以同时具有克隆和状态注释的细胞为分母。
3. 表型流量估计克隆层面状态分布的变化，不是对真实细胞进行的谱系追踪，也不能直接解释成实际分化轨迹。
4. Fisher 状态富集是样本内探索性推断；生物学组间效应需要供体/样本层面的模型。
5. 克隆计数为零或缺失表示在当前阈值下未检出，不证明克隆灭绝。测序深度为零的样本不能用于扩增检验。
6. 纵向检验在各供体内独立运行；随访时间点选择先进行 Bonferroni 校正，再对受检的供体/克隆进行 BH 校正。继承的 modified-Fisher 统计是 CloneTrack 风格方法，不是通用重复测量组间效应模型。
7. 原生纵向样本元数据要求每个供体/时间点对应一个样本，技术重复文库应提前聚合。
8. bulk 匹配默认按相同样本 ID 关联。跨时间或跨组织关联需要明确的 `sample_pairs`；两侧可用的供体信息会检查冲突。
9. 单链匹配可能有歧义。同一细胞有多个 bulk 比较或多条匹配链时，不能静默压缩成一个注释。
10. 表达推断通过 cellscope/PyDESeq2，使用原始 counts 聚合及样本层面设计。细胞不能作为供体的独立重复。

## 7. 依赖、复用与兼容性

immune base dependencies are NumPy, Pandas and SciPy. Native receptor objects use `immune[singlecell]`; RNA inference is installed separately through cellscope. Repertoire plotting and statistical backends remain optional.

不启动 VDJtools 或 R 组库包命令。可借鉴 immunarch 的功能组织和分析思路，但执行流程统一为 Python。

设计和后端接口参考：

- [Scirpy 受体数据结构](https://scirpy.scverse.org/en/latest/data-structure.html)。
- [Scirpy 克隆分析源码](https://github.com/scverse/scirpy/blob/main/src/scirpy/tl/_clonotypes.py)。
- [Scirpy](https://github.com/scverse/scirpy)、[Scanpy](https://github.com/scverse/scanpy)、[MuData](https://github.com/scverse/mudata)。
- [decoupler](https://github.com/saezlab/decoupler-py)、[PyDESeq2](https://github.com/owkin/PyDESeq2)。

本次实现没有复制第三方源码文件。以后复制源码需要检查许可证并保留相关声明。本地 PhenoTrack/CloneTrack 参考代码保持独立且未修改。

## 8. 扩展原则与后续方向

ATAC、蛋白等模态可使用额外 AnnData 对象加入 MuData。空间 spot 如果不是单个细胞，需要建立 spot–cell/clone 映射，不能仅凭无关 barcode 字符串进行对齐。

初版尚未验证或尚未实现的方向包括：

- 抗原特异性数据库查询。
- 更完整的 BCR 突变和谱系分析。
- Sankey 等更丰富的可视化。
- 层次化供体模型与邻域差异丰度。
- 大数据性能优化。

新能力应继续通过统一模块访问结果，并遵守明确的数据约定。当前优先稳定 RNA/VDJ 联合分析，不把未来扩展写成已经完成的能力。
