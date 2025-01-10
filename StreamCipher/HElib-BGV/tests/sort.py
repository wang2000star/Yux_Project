import numpy as np
import pandas as pd

def read_data(file_path):
    # data = pd.read_csv(file_path, delim_whitespace=True)
    # return data
    return pd.read_csv(file_path, sep='\s+')

def remove_outliers(data, value_columns, z_threshold=3):
    for col in value_columns:
        col_zscore = (data[col] - data[col].mean()) / data[col].std()
        data = data[(col_zscore < z_threshold) & (col_zscore > -z_threshold)]
    return data

def calculate_grouped_averages(data, group_columns, value_columns):
    grouped = data.groupby(group_columns)
    averages = grouped[value_columns].mean().reset_index()
    return averages

# def write_averages_to_file(headers, averages, output_file_path, column_width=20):
#     with open(output_file_path, 'w') as file:
#         # 写入指标名称
#         header_str = "".join(f"{header:<{column_width}}" for header in headers)
#         file.write(header_str + "\n")
#         # 写入平均值
#         for _, row in averages.iterrows():
#             avg_str = "".join(f"{value:<{column_width}.6f}" for value in row)
#             file.write(avg_str + "\n")
def write_to_file(headers, data, file_path):
    with open(file_path, 'w') as f:
        # 写入标题
        f.write(' '.join(f'{header:<10}' for header in headers) + '\n')
        # 写入数据
        for row in data.itertuples(index=False, name=None):
            f.write(' '.join(f'{str(item):<10}' for item in row) + '\n')

def main():
    input_file_path = 'test_context.txt'
    output_file_path = 'test_context_sort.txt'
    
    data = read_data(input_file_path)

    #对数据进行排序
    data = data.sort_values(by=['p', 'nslots', 'Qbits', 'c'])
    

    
    headers = list(data.columns)
    write_to_file(headers, data, output_file_path)
    print(f"计算结果已输出到 {output_file_path}")

    # # 读取txt文件
    # df = pd.read_csv(output_file_path, sep='\s+')
    # # 将DataFrame保存为Excel文件
    # df.to_excel(excel_file_path, index=False)

if __name__ == "__main__":
    main()