# Pipeline现状核查（基于当前仓库与`results/`产物）

本报告基于当前仓库代码与`run_all.py`已生成的`data_raw/`、`data_processed/`、`results/`内容。

## 1) 三个目录存储内容

- `data_raw/`：原始输入文件（GEO矩阵、counts、normalized CPM、GMT基因集）。
- `data_processed/`：预处理中间结果（manifest、过滤后表达矩阵、log2CPM、清洗后的immune基因集、Ensembl-Symbol映射缓存）。
- `results/`：最终输出（`tables/`统计结果、`figures/`图像、`logs/`运行日志）。

## 2) 是否“完成了生信分析”

结论：
- **流程在技术上跑完了（12/12脚本完成）**。
- **但生物学结论链条没有闭环，属于“部分完成”**。

关键证据：
- 主流程日志显示全部DONE。
- 但关键下游文件为空：`inflammatory_core_genes.csv`、`progressive_inflammatory_genes.csv`、`lasso_selected_genes.csv`、`signature_scores.csv`、`gsea_summary_matrix.csv`、`ipa_input_selected_genes.csv`、`gene_immune_correlations.csv`均为0行。
- 根因在于`deg_primary_healthy_vs_npdr_pdr_dme.csv`中显著基因数为0，导致后续候选基因、LASSO、富集等链式清空。

## 3) results/figures打不开或只有标题的根因

- **打不开**：来自`03_qc_and_pca.py`，直接把文本`placeholder ...`写进了`.png/.pdf`文件（这不是合法图像/PDF）。
- **只有标题**：来自`10_make_figures.py`，即便matplotlib可用，也只绘制了`ax.text(..., name)`标题占位图，并未读取统计表绘制真实图。

因此问题并非单点bug，而是脚本设计本身就是“占位实现”。

## 4) 图像问题是代码问题还是生信问题

两者都有，但主因是代码层：
- 代码层：绘图脚本明确是placeholder逻辑，不会生成真实图。
- 生信层：关键中间结果为空（例如候选基因链路），即使替换绘图代码，也会有多张图无数据可画。

## 5) 完成论文级分析仍缺什么

至少缺以下真实实现：
1. `03_qc_and_pca.py`：真实QC/PCA统计与绘图（而非文本placeholder）。
2. `10_make_figures.py`：按`results/tables/*.csv`真实作图，而非标题占位图。
3. `07_lasso_signature.py`：当前是fallback pseudo-lasso和固定AUC占位；需改为真实LASSO+OOF ROC。
4. `08_enrichment_analysis.py`：当前是“单基因是否属于pathway”的简化占位，非标准GSEA。
5. `04_differential_expression.py`：当前统计策略需复核（显著基因为0导致下游全部断链）。

