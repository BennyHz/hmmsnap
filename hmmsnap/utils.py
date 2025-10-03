# hmmcall/utils.py
# 标准染色体列表（两种格式）
STD_CHR_WITH_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
STD_CHR_WITHOUT_CHR = [f"{i}" for i in range(1, 23)] + ["X", "Y"]
STD_CHR_MAPPING = {f"{i}": f"chr{i}" for i in range(1, 23)}
STD_CHR_MAPPING.update({"X": "chrX", "Y": "chrY"})

GENOME_SIZES = {
    "hg38": 2652463870,
    "hg19": 2451960000,
    "mm10": 2150570000,
    "mm39": 2150570000
}