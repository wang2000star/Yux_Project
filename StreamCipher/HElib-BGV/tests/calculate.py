import numpy as np
import pandas as pd

def read_data(file_path):
    data = pd.read_csv(file_path, delim_whitespace=True)
    return data

def remove_outliers(data, value_columns, z_threshold=3):
    for col in value_columns:
        col_zscore = (data[col] - data[col].mean()) / data[col].std()
        data = data[(col_zscore < z_threshold) & (col_zscore > -z_threshold)]
    return data

def calculate_grouped_averages(data, group_columns, value_columns):
    grouped = data.groupby(group_columns)
    averages = grouped[value_columns].mean().reset_index()
    return averages

def write_averages_to_file(headers, averages, output_file_path, column_width=20):
    with open(output_file_path, 'w') as file:
        # 写入指标名称
        header_str = "".join(f"{header:<{column_width}}" for header in headers)
        file.write(header_str + "\n")
        # 写入平均值
        for _, row in averages.iterrows():
            avg_str = "".join(f"{value:<{column_width}.6f}" for value in row)
            file.write(avg_str + "\n")

def main():
    input_file_path = 'test_Yus_p_C64_ClientAndServer6.txt'
    output_file_path = 'test_Yus_p_C64_ClientAndServer6_calculate.txt'
    excel_file_path = 'test_Yus_p_C64_ClientAndServer6_calculate.xlsx'
    
    data = read_data(input_file_path)
    
    # 指定分组列和计算平均值的列
    group_columns = ['Nr', 'p', 'nslots', 'bits', 'c']
    value_columns = [col for col in data.columns if col not in group_columns]
    
    # 去掉离群点
    data = remove_outliers(data, value_columns)
    
    averages = calculate_grouped_averages(data, group_columns, value_columns)
    
    headers = list(averages.columns)
    write_averages_to_file(headers, averages, output_file_path)
    print(f"计算结果已输出到 {output_file_path}")

    # 读取txt文件
    df = pd.read_csv(output_file_path, sep='\s+')
    # 将DataFrame保存为Excel文件
    df.to_excel(excel_file_path, index=False)

if __name__ == "__main__":
    main()