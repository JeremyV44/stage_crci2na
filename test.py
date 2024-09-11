import pandas as pd 
import pickle 
import anndata
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

with open('286.pkl', 'rb') as f:
    df = pickle.load(f)

    filtered_data = df[df.obs['Condition'] == '4T1-ctrl']
    new_filtered = filtered_data.to_df()
    new_filtered.T.to_csv('filtered_data.csv')
