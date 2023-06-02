import pandas as pd
#reading in orthology tables from ensembl database
ensembl_human=pd.read_csv("ens_to_hg.tsv",sep="\t") #contains Ensembl Gene stable ID and HGNC symbol
ensembl_zfish=pd.read_csv("ens_to_zfin.tsv",sep="\t") #contains Ensembl Zebrafish gene stable ID and ZFIN symbol
ensembl_human_to_zfish=pd.read_csv("ensg_to_ensdarg.ensembl105.tsv",sep="\t") #contains orthology ENSG and ENSDARG from ensembl 
oma_human_to_zfish=pd.read_csv("oma_enshg_ensdarg.tsv",sep="\t",header=None) #contains orthology ENSG and ENSDARG from oma orthology browser
alliance_human_to_zfish=pd.read_csv("ORTHOLOGY-ALLIANCE_COMBINED.tsv",skiprows=15,sep="\t") #contains orthology symbols from human and zfin

aerts_tbl=pd.read_csv("v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl",sep="\t") #human .tbl file

#create an ensembl orthologue table
human_zf = ensembl_human.merge(ensembl_human_to_zfish, on='Gene stable ID', how='left')
ensembl_zfish=ensembl_zfish.rename(columns={'Gene stable ID':'Zebrafish gene stable ID'})
zf_zf = ensembl_zfish.merge(ensembl_human_to_zfish, on='Zebrafish gene stable ID', how='left')
human_zf_grouped = human_zf.groupby('HGNC symbol').agg({'Gene stable ID': lambda x: tuple(x), 'Zebrafish gene stable ID': lambda x: tuple(x)})
human_zf_grouped["hg_symbol"]=human_zf_grouped.index
zf_zf_grouped = zf_zf.groupby('ZFIN symbol').agg({'Gene stable ID': lambda x: tuple(x), 'Zebrafish gene stable ID': lambda x: tuple(x)})
zf_zf_grouped["zf_symbol"]=zf_zf_grouped.index
human_zf_mapping = human_zf_grouped.merge(zf_zf_grouped, on='Gene stable ID',how="left")
ens_human_zf_mapping=human_zf_mapping.dropna(subset=['hg_symbol', 'zf_symbol'])

#create an oma orthology table
oma_human_to_zfish=oma_human_to_zfish.rename(columns={0:'ens_hg',1:'ens_zf',2:'n_orth',3:'score'})
oma_human_to_zfish=oma_human_to_zfish.merge(ensembl_human,left_on="ens_hg",right_on='Gene stable ID')
oma_human_to_zfish=oma_human_to_zfish.merge(ensembl_zfish,left_on="ens_zf",right_on='Zebrafish gene stable ID')
oma_human_to_zfish=oma_human_to_zfish.dropna(subset=['HGNC symbol', 'ZFIN symbol'])
oma_human_to_zfish=oma_human_to_zfish.rename(columns={'HGNC symbol':'hg_symbol','ZFIN symbol':'zf_symbol'})
#prepare alliance db
alliance_human_to_zfish=alliance_human_to_zfish[alliance_human_to_zfish["Gene1SpeciesName"]=="Homo sapiens"] 
alliance_human_to_zfish=alliance_human_to_zfish[alliance_human_to_zfish["Gene2SpeciesName"]=="Danio rerio"]
alliance_human_to_zfish=alliance_human_to_zfish.rename(columns={'Gene1Symbol':'hg_symbol','Gene2Symbol':'zf_symbol'})
#combine orthologies
gene_orthologies=pd.concat([ens_human_zf_mapping.iloc[:,[2,4]],alliance_human_to_zfish.iloc[:,[1,5]],oma_human_to_zfish.iloc[:,[5,7]]])
gene_orthologies=gene_orthologies.drop_duplicates()
#change gene_names
aerts_tbl=aerts_tbl.merge(gene_orthologies,left_on="gene_name",right_on='hg_symbol')
aerts_tbl['gene_name']=aerts_tbl['zf_symbol']
aerts_tbl=aerts_tbl.drop(columns=['hg_symbol', 'zf_symbol'])
#save final file
aerts_tbl.to_csv("motifs-v10-zfish.tbl")
