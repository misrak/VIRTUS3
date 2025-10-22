import pandas as pd

file = './data/NC_007605.1_CDS.fa'
file_out = './data/NC_007605.1_CDS.tgMap.tsv'

with open(file, 'r') as f:
    txt = f.readlines()

txt = [x for x in txt if x.startswith('>')]

tx = [x.split()[0][1:] for x in txt]
gene = [x.split()[1].replace('[gene=', '').replace(']', '') for x in txt]

pd.DataFrame([tx, gene]).T.to_csv(file_out, header=None, index=None, sep='\t')