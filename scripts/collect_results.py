import pandas as pd
import sys
import os


model = sys.argv[1]
run = sys.argv[2]

repo_root = os.path.dirname(os.path.abspath(__file__))
result_dir = os.path.join(repo_root, '../results/'+model+'/run'+run)

print(result_dir)

# List all CSV files in the directory
csv_files = [f for f in os.listdir(result_dir) if f.endswith('.csv')]

# List to store each DataFrame
dataframes = []

# Read each file and append to the list
for file in csv_files:
    file_path = os.path.join(result_dir, file)
    df = pd.read_csv(file_path)
    dataframes.append(df)

# Concatenate all DataFrames
combined_df = pd.concat(dataframes, ignore_index=True)
combined_df = combined_df.sort_values(by=['index'])
print(len(combined_df))
# Write the combined DataFrame to a new CSV file

file_name = "result"+run+".csv"
result_dir = os.path.join(repo_root, '../results/'+model)
file_path = os.path.join(result_dir, file_name)
combined_df.to_csv(file_path, index=False)
