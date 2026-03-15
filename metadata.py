import pandas as pd
import re

def build_filtered_metadata(matrix_file, counts_file, output_file):
    print(f"正在解析 GEO 临床数据文件: {matrix_file} ...")
    
    sample_data = {}
    
    # 1. 逐行解析 Series Matrix 文件
    with open(matrix_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # 提取样本名称 (通常对应 sample_xx)
            if line.startswith('!Sample_title'):
                titles = line.split('\t')[1:]
                titles = [t.strip('"') for t in titles]
                for i, title in enumerate(titles):
                    sample_data[i] = {'SampleID': title}
            
            # 提取样本临床特征
            elif line.startswith('!Sample_characteristics_ch1'):
                chars = line.split('\t')[1:]
                chars = [c.strip('"') for c in chars]
                for i, char in enumerate(chars):
                    if i in sample_data and ':' in char:
                        key, value = char.split(':', 1)
                        key = key.strip().lower()
                        value = value.strip()
                        # 将特征动态写入字典
                        sample_data[i][key] = value

    # 转化为 DataFrame
    df_meta_full = pd.DataFrame.from_dict(sample_data, orient='index')
    
    if df_meta_full.empty:
        raise ValueError("未能从 Series Matrix 中提取到有效信息，请检查文件格式。")

    print(f"提取到原始样本总数: {len(df_meta_full)}")

    # 2. 清洗数据并建立标准化标签
    # 提取 Region 和 Group，GEO 的 key 可能有细微不同，这里做一定的容错匹配
    region_col = [col for col in df_meta_full.columns if 'region' in col]
    disease_col = [col for col in df_meta_full.columns if 'disease' in col or 'diagnosis' in col]
    
    if not region_col or not disease_col:
        print("当前提取的列名：", df_meta_full.columns.tolist())
        raise KeyError("未能找到 'region' 或 'disease state' 相关字段，请人工检查 matrix 文件中的确切属性名称。")

    region_key = region_col[0]
    disease_key = disease_col[0]

    # 将疾病状态映射为你的代码需要的 Control 和 DR
    def map_group(disease_state):
        state = str(disease_state).lower()
        if 'normal' in state or 'healthy' in state or 'control' in state:
            return 'Control'
        elif 'diabetic retinopathy' in state or 'dr' in state:
            return 'DR'
        return 'Other'

    df_meta_full['Group'] = df_meta_full[disease_key].apply(map_group)
    df_meta_full['Region_Clean'] = df_meta_full[region_key].apply(lambda x: str(x).lower())

    # 3. 严格执行规划：过滤出黄斑区 (Macula) 且标签明确的样本
    print("正在执行过滤：仅保留黄斑区 (Macula) 样本...")
    df_filtered = df_meta_full[
        (df_meta_full['Region_Clean'].str.contains('macula')) & 
        (df_meta_full['Group'].isin(['Control', 'DR']))
    ].copy()

    print(f"过滤后符合条件的 Macula 样本数: {len(df_filtered)}")

    # 4. 对齐 counts 矩阵
    # 验证提取的 SampleID 是否真的存在于 counts_file 的列名中
    try:
        counts_df = pd.read_csv(counts_file, sep='\t', nrows=0)
        counts_columns = set(counts_df.columns)
    except Exception as e:
        raise FileNotFoundError(f"无法读取 counts 文件以验证列名: {e}")

    # 检查重合度
    valid_samples = []
    for sid in df_filtered['SampleID']:
        # 有时候 GEO title 会包含冗余信息，这里做一个精准匹配
        if sid in counts_columns:
            valid_samples.append(sid)
        else:
            # 尝试通过正则找 sample_xx
            match = re.search(r'(sample_\d+)', sid, re.IGNORECASE)
            if match and match.group(1).lower() in counts_columns:
                valid_samples.append(match.group(1).lower())

    df_filtered['SampleID_Matched'] = valid_samples

    # 5. 组装最终可用的 metadata
    final_metadata = df_filtered[['SampleID_Matched', 'Group']].rename(columns={'SampleID_Matched': 'SampleID'})
    
    # 保存结果
    final_metadata.to_csv(output_file, index=False)
    print(f"\n成功生成标准化临床信息文件: {output_file}")
    print(f"最终纳入差异分析队列: Control 组 {len(final_metadata[final_metadata['Group']=='Control'])} 例，DR 组 {len(final_metadata[final_metadata['Group']=='DR'])} 例。")
    print("现在，你可以安全地运行前面提供的生信主流程代码了。")

if __name__ == "__main__":
    matrix_file = 'GSE160306_series_matrix.txt'
    counts_file = 'GSE160306_human_retina_DR_totalRNA_counts.txt'
    output_file = 'metadata.csv'
    
    build_filtered_metadata(matrix_file, counts_file, output_file)