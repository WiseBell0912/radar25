import os
import pandas as pd
import re

# 엑셀 파일들이 있는 디렉터리 (현재 스크립트 기준)
folder = os.path.dirname(os.path.abspath(__file__))

# 파일 패턴에 맞는 것만 선택 + 정렬
excel_files = [
    f for f in os.listdir(folder)
    if f.endswith('.xlsx') and 'data_' in f and '_KR.xlsx' in f
]

# 정렬: YYYYMM 추출 후 정렬
def extract_yyyymm(filename):
    match = re.search(r'_(\d{6})_KR\.xlsx$', filename)
    return int(match.group(1)) if match else float('inf')

excel_files = sorted(excel_files, key=extract_yyyymm)

# 병합
merged_df = pd.DataFrame()
for file in excel_files:
    df = pd.read_excel(os.path.join(folder, file))
    merged_df = pd.concat([merged_df, df], ignore_index=True)

# 저장
output_path = os.path.join(folder, 'GYEONGPO_merged_2018_2023.xlsx')
merged_df.to_excel(output_path, index=False)
print(f"✅ 병합 완료: {output_path}")
