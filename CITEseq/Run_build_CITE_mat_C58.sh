#!/bin/bash
#SBATCH -t 6-18:00:00
#SBATCH --mem=60000M 
#SBATCH -J c58_cite_processing
#SBATCH -p himem 
#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH -o %x-%j.out


cd /cluster/projects/macparland/TA/LiverMap2.0/RawData/C58_3pr_CITE_ADT
READ1=/cluster/projects/macparland/TA/LiverMap2.0/RawData/C58_3pr_CITE_ADT/C58_CITE_ADT_S9_L008_R1_001.fastq
READ2=/cluster/projects/macparland/TA/LiverMap2.0/RawData/C58_3pr_CITE_ADT/C58_CITE_ADT_S9_L008_R2_001.fastq

TRIPLET=/cluster/projects/macparland/TA/LiverMap2.0/RawData/C58_3pr_CITE_ADT/C58_CITE_ADT_L008_Triplet.out

SCRIPT1=/cluster/home/tandrews/scripts/LiverMap2.0/CITEseq/Make_Triplets_ADT_CITE_existingbarcodes.pl
SCRIPT2=/cluster/home/tandrews/scripts/LiverMap2.0/CITEseq/make_sparse_mat.R
BARCODE=/cluster/projects/macparland/TA/LiverMap2.0/RawData/C58_RESEQUENCED_JUL_20/raw_feature_bc_matrix/C58_3pr_gt250_barcodes.txt


perl $SCRIPT1 $READ1 $READ2 $BARCODE > $TRIPLET

/cluster/tools/software/R/3.5.0/Rscript $SCRIPT2 $TRIPLET C58_CITE_ADT_gt250_all_Matrix.rds
