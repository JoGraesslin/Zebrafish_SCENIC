# Zebrafish_SCENIC
## Short Instruction on how to create a zfish .tbl file for SCENIC
The goal is to rename the gene names from human to zebrafish. We use the .tbl file with motifs published by the Aerts Lab ```v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl```

Download human to ens here https://www.ensembl.org/biomart/martview/. Get ENS to human symbols, ENSDARG to zfin symbols, and the orthology tables (at dataset)
Download alliance database [
](https://www.alliancegenome.org/downloads#orthology)

Download oma database https://omabrowser.org/oma/home/ [
](https://www.alliancegenome.org/downloads#orthology)

Change the filepath accordingly and run script create_zfish_tbl.py
