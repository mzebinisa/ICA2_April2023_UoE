**** OQSIL - Online Protein Quality aSsesIng tooL 
 
 This tool was developed as part of the coursework for the Bioinformatics programming and systems management module at the Bioinformatics postgraduate program at the University of Edinburgh in 2022-2023.

 **Introduction**

OQSIL is a generic Python-based tool designed to assist biologists and protein-focused researchers in conducting comprehensive protein sequence analyses. It allows users to download protein FASTA sequences from NCBI for specific organisms or taxa, then performs a range of quantitative and qualitative analyses. These include amino acid count, sequence alignment, conservation plotting, similarity scoring, motif identification, and physicochemical property analysis of individual proteins.


 Features

- Download protein FASTA sequences: Automatically retrieve protein sequences from NCBI for a given protein family and taxonomic group.
- Amino acid and sequence analysis: Perform analyses like amino acid count, sequence length distribution, and sequence alignment.
- Conservation plotting: Visualize the level of conservation between protein sequences.
- Motif identification: Scan protein sequences for known motifs from the PROSITE database.
- Physicochemical analysis: Calculate and plot various physicochemical properties of proteins.

 Prerequisites

To use OQSIL, you will need the following:

1. MSC8 server: The program is designed to run on the MSC8 server. Make sure you are working in this environment.
2. Required software: If you don't have access to the MSC8 server, you need to install the following tools:
   - Entrez utilities (for downloading sequences from NCBI)
   - Clustal Omega (for sequence alignment)
   - EMBOSS (for motif and physicochemical analysis)

These tools can typically be installed on most systems via package managers or by downloading directly from their respective websites.

 Setup Instructions

1. Download the file oqsil.py or tar genome-qc-uoe-pg-project.tar.gz to your working directory.

Once the file is downloaded, extract the files: using tar -xvf genome-qc-uoe-pg-project.tar.gz
Read through and follow the manual to perform analysis.
