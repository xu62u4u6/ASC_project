import pandas as pd

target_col = [
        "Hugo_Symbol",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2"
    ]
 

def read_maf(maf_path, case_ID, preffix: str="",suffix: str="") -> pd.DataFrame:
    maf = pd.read_csv(maf_path, skiprows=1, sep="\t")
    maf["case_ID"] = f"{preffix}{case_ID}{suffix}"
    maf.index = maf.loc[:, target_col].apply(lambda row: "|".join(row.astype(str)), axis=1)
    maf = filter_vaild_variant_classfication(maf) # only keep vaild mutations
    return maf

def filter_vaild_variant_classfication(maf: pd.DataFrame):
    '''
    GDC Reference:
    https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/
    MAF file fields
    '''

    vaild_variant_classfication = [
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Silent",
        "Splice_Site",
        "Translation_Start_Site",
        "Nonstop_Mutation",
        "3'UTR",
        "3'Flank",
        "5'UTR",
        "5'Flank",
        "IGR",
        "Intron",
        "RNA",
        "Targeted_Region"
    ]
    return maf[maf.Variant_Classification.isin(vaild_variant_classfication)]

def filter_non_synonymous_variant_classfication(maf: pd.DataFrame):
    nonsynonymous_types = [
        "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
        "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
        "Translation_Start_Site", "Nonstop_Mutation"
    ]
    return maf[maf.Variant_Classification.isin(nonsynonymous_types)]

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
            return "Multi_Hit"
        elif len(unique_mutations) == 1:
            return unique_mutations[0]
    return False

def all_case_maf_to_pivot_table(all_case_maf: pd.DataFrame) -> pd.DataFrame: 
    pivot_table =  all_case_maf.pivot_table(
                        values="Variant_Classification",
                        index="Hugo_Symbol",
                        columns="case_ID",
                        aggfunc=merge_mutations
                        ).fillna(False)
    return pivot_table

def calculate_frequency(df: pd.DataFrame) -> pd.Series:
    return (df != False).sum(axis=1) / df.shape[1]

def sort_by_variant_frequency(merged_df: pd.DataFrame, sample_type_list: list = None) -> pd.DataFrame:
    """
    sort genes by variant_frequency

    Parameters:
    -----------
    merged_df : pd.DataFrame
    sample_type_list: 

    Returns:
    --------
    sorted_df: pd.DataFrame
    """

    result_df = merged_df.copy()
    # add sample type
    if sample_type_list != None:
        for sample_type in sample_type_list:
            sub_df = merged_df.loc[:, merged_df.columns.str.contains(f"_{sample_type}")]
            result_df[f"{sample_type}_freq"] = calculate_frequency(sub_df)
    result_df["all_freq"] = calculate_frequency(merged_df)

    # sort by freq
    result_df = result_df.sort_values('all_freq', ascending=False)
    return result_df

def binary_sort_key(column: pd.Series) -> int: 
    # binary column to int  
    binary_str = "".join(column.astype(int).astype(str))
    return int(binary_str, 2)

def sort_samples(pivot_table, top=10, freq_columns=["all_freq"]):
    """
    Sort samples by convert column to binary 
    """
    # check if freq columns exist
    if not all(col in pivot_table.columns for col in freq_columns):
        raise ValueError("Some freq_columns are not present in the pivot_table")
    
    tmp_pivot_table = pivot_table.drop(columns=freq_columns)
    binary_pivot_table = tmp_pivot_table != False
    sort_order = (binary_pivot_table.head(top)
                  .apply(binary_sort_key, axis=0)
                  .sort_values(ascending=False)  
                  .index)                        
    
    # sort by order
    sorted_samples = tmp_pivot_table[sort_order]
    
    # concat freq columns
    sorted_table = pd.concat([sorted_samples, pivot_table[freq_columns]], axis=1)
    
    return sorted_table