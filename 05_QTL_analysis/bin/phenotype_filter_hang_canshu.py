import pandas as pd
import numpy as np
import argparse

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description='清理表型数据，将超出阈值的点替换为NA。')

# 添加参数
parser.add_argument('--input_file_path', type=str, required=True, help='输入文件的路径。')
parser.add_argument('--output_file_path', type=str, required=True, help='输出文件的路径。')
parser.add_argument('--threshold', type=float, default=3, help='标准差的阈值，默认为3。')
parser.add_argument('--na_rep', type=str, default='NA', help='替换NA的字符串，默认为"NA"。')

# 解析参数
args = parser.parse_args()

# 加载数据集，假设文件是CSV格式
df = pd.read_table(args.input_file_path, index_col=0)

# 对于每一行（表型），计算离散点
for phenotype in df.index:
    mean_value = df.loc[phenotype].mean()
    std_value = df.loc[phenotype].std()
    # 计算上下阈值
    upper_bound = mean_value + args.threshold * std_value
    lower_bound = mean_value - args.threshold * std_value
    # 找到离散点并替换为NA
    df.loc[phenotype, (df.loc[phenotype] > upper_bound) | (df.loc[phenotype] < lower_bound)] = np.nan

# 显示替换后的DataFrame
print(df)

# 保存清洗后的数据集
df.to_csv(args.output_file_path, sep='\t', na_rep=args.na_rep)

