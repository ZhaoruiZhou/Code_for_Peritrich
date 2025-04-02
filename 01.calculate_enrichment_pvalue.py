# -*- coding: UTF-8 -*-
# by zzr 2023.11.9 
# for calculate the enrichment of specific genes within gene clusters


#################################******************###############################
import pandas as pd
from scipy.stats import hypergeom

##################################################################################
def read_file1(file1):
    df = pd.read_csv(file1, sep='\t', header=None, names=['Gene', 'Cluster'])
    return df

def read_file2(file2):
    with open(file2) as f:
        subset_genes = set(line.strip() for line in f)
    return subset_genes

def compute_enrichment(df, subset_genes):
    total_genes = len(df)
    subset_in_file1 = df[df['Gene'].isin(subset_genes)]
    cluster_counts = df['Cluster'].value_counts().to_dict()
    cluster_subset_counts = subset_in_file1['Cluster'].value_counts().to_dict()
    
    results = []
    for cluster, M in cluster_counts.items():
        n = cluster_subset_counts.get(cluster, 0) 
        N = len(subset_genes)  
        K = M
        
        p_value = hypergeom.sf(n-1, total_genes, K, N)
        enrichment_score = (n / N) / (K / total_genes) if K > 0 else 0
        
        results.append([cluster, n, K, enrichment_score, p_value])
    
    result_df = pd.DataFrame(results, columns=['Cluster', 'Subset_Count', 'Cluster_Size', 'Enrichment_Score', 'P_Value'])
    result_df = result_df.sort_values(by='P_Value')
    return result_df

############################################################################################
file1 = "mapping_file.txt"  #The first column is the gene name, and the second column is the cluster number
file2 = "gene_id.txt"     # gene id
############################################################################################

df = read_file1(file1)
subset_genes = read_file2(file2)
result_df = compute_enrichment(df, subset_genes)
print(result_df)
