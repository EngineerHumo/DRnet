import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# 1. 加载临床信息 (需要你自己准备)
# metadata 应该包含 SampleID (如 sample_09) 和 Group (Control, DR) 两列
try:
    meta = pd.read_csv('metadata.csv', index_col='SampleID')
except FileNotFoundError:
    raise Exception("严重错误：必须先创建包含 SampleID 和 Group 两列的 metadata.csv，并仅保留 Macula 样本。")

# 2. 读取 Counts 和 CPM 矩阵
# 注意：你的counts文件带有 ensemblID 和样本列
counts_df = pd.read_csv('GSE160306_human_retina_DR_totalRNA_counts.txt', sep='\t', index_col=0)
cpm_df = pd.read_csv('GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt', sep='\t', index_col=0)

# 3. 严格对齐样本 (使用 metadata 中的样本 ID 过滤和排序矩阵)
valid_samples = meta.index.intersection(counts_df.columns)
if len(valid_samples) == 0:
    raise ValueError("SampleID 无法匹配，请检查 metadata.csv 的命名是否与列名一致。")

# 过滤并排序矩阵
counts_df = counts_df[valid_samples]
cpm_df = cpm_df[valid_samples]
meta = meta.loc[valid_samples]

# 4. 对 CPM 进行 log2(CPM + 1) 转换，用于下游可视化、降维和 LASSO
# 检查是否已经转换过，避免重复 log
if cpm_df.max().max() > 50:
    log_cpm = np.log2(cpm_df + 1)
else:
    log_cpm = cpm_df