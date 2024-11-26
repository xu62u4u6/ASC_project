import pandas as pd

target_col = [
        "Hugo_Symbol",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2"
    ]


        
#def merge_mafs(maf_list: list[MAF]):
#    return MAF(pd.concat(maf_list, axis=0, ignore_index=True))

def read_maf(maf_path, case_ID, preffix="",suffix=""):
    maf = pd.read_csv(maf_path, skiprows=1, sep="\t")
    maf["case_ID"] = f"{preffix}{case_ID}{suffix}"
    maf.index = maf.loc[:, target_col].apply(lambda row: "|".join(row.astype(str)), axis=1)
    return maf

def filter_maf(maf: pd.DataFrame, 
               filter_mutation_list: list=["Silent", "3'UTR", "5'UTR", "IGR", "Intron", "RNA"]):
    return maf[~maf.Variant_Classification.isin(filter_mutation_list)]

def merge_mutations(column):
    """
    Merge mutations in a column:
    - If all mutations are the same type, return that type
    - If there are different types, return "multihit"
    - If no mutations, return False
    
    Args:
        column: pandas Series containing mutation types
    Returns:
        str or bool: merged mutation type or False
    """
    if column.any():
        # Get unique non-False mutation types
        unique_mutations = column[column != False].unique()
        if len(unique_mutations) > 1:
            return "multihit"
        elif len(unique_mutations) == 1:
            return unique_mutations[0]
    return False

def all_case_maf_to_matrix(all_case_maf):
    merged_mutations = all_case_maf.groupby(["case_ID", "Hugo_Symbol"])['Variant_Classification'].agg(
            merged_type=merge_mutations
        ).reset_index()

    mutation_matrix = merged_mutations.pivot(
        columns='case_ID',
        index='Hugo_Symbol',
        values='merged_type'
    ).fillna(False)
    return mutation_matrix

def calculate_frequency(df: pd.DataFrame) -> pd.DataFrame:
    return (df != False).sum(axis=1) / df.shape[1]

def sort_by_variant_frequency(merged_df: pd.DataFrame, sample_type_list: list = None) -> pd.DataFrame:
    """
    按照變異出現頻率對DataFrame進行排序

    Parameters:
    -----------
    merged_df : pd.DataFrame
        合併後的變異DataFrame

    Returns:
    --------
    pd.DataFrame
        按照變異頻率排序後的DataFrame
    """

    # 將頻率添加為新的列，便於查看
    result_df = merged_df.copy()
    if sample_type_list != None:
        for sample_type in sample_type_list:
            sub_df = merged_df.loc[:, merged_df.columns.str.contains(f"_{sample_type}")]
            result_df[f"{sample_type}_freq"] = calculate_frequency(sub_df)
    result_df["all_freq"] = calculate_frequency(merged_df)

    # 按照頻率降序排列
    result_df = result_df.sort_values('all_freq', ascending=False)
    return result_df