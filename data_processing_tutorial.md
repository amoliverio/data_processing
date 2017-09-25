#Data Processing in the Fierer Lab: Sequence Variants (zOTUs) or Clustering at 97% identity (OTUs)

**Updated: 09/25/2017**  
**Questions: angela.oliverio43@gmail.com**

##Notes:

This pipeline for processing raw sequence data into sequence variants or OTUs with USEARCH and UNOISE3/UPARSE. This document is an updated version of data processing tutorial [here](https://github.com/leffj/data-tutorials/blob/master/amplicon_data_processing_tutorial/amplicon_data_processing-16S.md). If this is your first time working on the microbe server please go [here](https://github.com/leffj/data-tutorials/blob/master/amplicon_data_processing_tutorial/amplicon_data_processing-16S.md) for an overview of how to login onto the server etc.:  

**###Step 1: Get tutorial data and check mapping file (e.g. metadata)**

Demo mapping file is available here: http://fiererlab.org/?p=516

First, check the mapping file to make sure it is formatted correctly (QIIME command):

	validate_mapping_file.py -m Demo_16S_MappingFile.txt -o checkout/

	less checkout/Demo_16S_MappingFile.log


**###Step 2: Make directory for sequence data**

**Important** Raw data is deposited into `/BioFrontiers`, and should be copied into `/data/shared/`, where all Fierer lab raw sequence data are be stored long term. 

For example, to see the data for this tutorial:

	ls /data/shared/2014_02_03_data_tutorial/

Note the naming convention for directories in `/data/shared`: YYYY\_MM\_DD\_Name

###Step 2a: If your fragments were sequenced through the adapters
**This applies to you if your sequences were 200bp or longer**
**Otherwise, skip to 'Step 2b: Demultiplexing'**

The issue is that the DNA fragments have synthetic adapters on both ends due to the sequencing process. If the sequences are long enough, they will extend through the adapter on the far (3') end of the fragment. It is important to trim these adapter sequences so that you are analyzing real biological data without artificial sequences.

To do this, I recommend using a program called 'cutadapt'. For this program to work, it needs to know the adapter sequences you expect to find in each read. Below are the primer/adapter sequences associated with the primers we typically use. However, it is always a good idea to check to make sure these are the correct sequences used in your library prep.

	16S rRNA gene forward primer sequence: TTACCGCGGCKGCTGGCACACAATTACCATAGTGTAGATCTCGGTGGTCGCCGTATCATT
	16S rRNA gene reverse primer sequence: ATTAGAWACCCBDGTAGTCCGGCTGACTGACT
	ITS forward primer sequence: TTACTTCCTCTAAATGACCAAGCCGTGTAGATCTCGGTGGTCGCCGTATCATT
	ITS reverse primer sequence: GCATCGATGAAGAACGCAGCATCTGACTGACT

Use the following command as a base to to trim adapters. This example assumes a mixed 16S/ITS run and removes both sets of adapters simultaneously. Note that we remove reverse adapter/primer sequences from the forward reads and vice versa for the reverse reads. This will not work with the tutorial dataset since the reads are not long enough. It takes approx. 20 min on a MiSeq dataset.

	cutadapt -a 16S_rev_primer=<16S reverse adapter/primer sequence> -a ITS_rev_primer=<ITS reverse adapter/primer sequence> -A 16S_fwd_primer=<16S forward adapter/primer sequence> -A ITS_fwd_primer=<ITS forward adapter/primer sequence> -O 1 -o <PATH/TO/R1/OUTPUT.fq> -p <PATH/TO/R2/OUTPUT.fq> Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz 1> <PATH/TO/cutadapt.log>

You can now use the output files from cutadapt with the normal raw index reads file to do demultiplexing and all other downstream processing steps.

###Step 2b: Demultiplexing

	prep_fastq_for_uparse_paired.py -i /data/shared/2014_02_03_data_tutorial/Undetermined_S0_L001_R1_001_t.fastq.gz -r /data/shared/2014_02_03_data_tutorial/Undetermined_S0_L001_R2_001_t.fastq.gz -b /data/shared/2014_02_03_data_tutorial/Undetermined_S0_L001_I1_001_t.fastq.gz -m Demo_16S_MappingFile.txt -o demultiplexed_seqs/ -c

**NOTE:** If you are dealing with sequences that are already demultiplexed by the Illumina software and are distributed in different files (one for each sample), see [here](https://github.com/leffj/data-tutorials/blob/master/amplicon_data_processing_tutorial/preparing_already_demultiplexed_data.md).

###Step 2c: Merging paired reads

Merge paired end reads using this command. Note that settings are just guesses and can be altered. **Do not** merge reads unless there is substantial overlap. **Do not** merge blindly or you could bias your data.

	usearch10 -fastq_mergepairs demultiplexed_seqs/demultiplexed_seqs_1.fq -reverse demultiplexed_seqs/demultiplexed_seqs_2.fq -fastq_minovlen 16 -fastq_minmergelen 200 -fastqout demultiplexed_seqs/demultiplexed_seqs_merged.fq -report demultiplexed_seqs/merge_rpt.txt
	less demultiplexed_seqs/merge_rpt.txt

###Step 2d: Strip primers after merging

**NOTE:** In some cases you may need to remove primers after merging, which can be done with the following command:

trim forward primer:
	
	cutadapt -g ^fwd_primer demultiplexed_seqs/demultiplexed_seqs_merged.fq > demultiplexed_seqs_merged_trimmedf.fq

trim reverse primer [note: reverse sequence is r-complement]
	
	cutadapt -a rev_primer demultiplexed_seqs/demultiplexed_seqs_merged_trimmedf.fq > demultiplexed_seqs_merged_trimmedfr.fq

**###Step 3: Prepare sequences for de novo database creation**

####Check quality of sequences
	
	usearch10 -fastq_eestats2 demultiplexed_seqs/demultiplexed_seqs_merged.fq -output demultiplexed_seqs/demultiplexed_seqs_merged_eestats2.log -length_cutoffs 200,300,10
	less demultiplexed_seqs/demultiplexed_seqs_merged_eestats2.log

####Conduct quality filtering

The -fastq_trunclen is an optional flag to truncate sequences at specific length, should be based off of eestats2 output!

	usearch10 -fastq_filter demultiplexed_seqs/demultiplexed_seqs_merged.fq -fastaout seqs_filt.fa -fastq_maxee 1.0 -fastq_trunclen 240

####Dereplicate sequences

	usearch10 -fastx_uniques seqs_filt.fa -fastaout uniques.fa -sizeout -relabel Uniq

**###Step 4: Make zOTUs or OTUs and create de novo database (e.g. rep_set), then filter against an existing public database to remove highly divergent sequences** 

####To create zOTUS (note unoise command has abundance threshold (-minsize), default is 8):
	
	usearch10 -unoise3 uniques.fa -zotus rep_set_zotus.fa -tabbedout unoise3.txt
	
	usearch10 -usearch_global rep_set_zotus.fa -db /db_files/gg_files/gg_13_8_otus/rep_set/99_otus.fasta -id 0.75 -strand both -matched rep_set_zotus_filt.fa
	
	-fastx_relabel rep_set_zotus_filt.fa -prefix 'OTU_' -fastaout rep_set_zotus_filt_relabeled.fa -keep_annots


####To cluster OTUs at 97% similarity (again, there is an abundance threshold (-minsize), default is 2):

	usearch10 -cluster_otus uniques.fa -otus rep_set97.fa -relabel 'OTU_'
	
	usearch10 -usearch_global rep_set97.fa -db /db_files/gg_files/gg_13_8_otus/rep_set/97_otus.fasta -id 0.75 -strand both -matched rep_set97_filt.fa


**###Step 5: Map the raw/demultiplexed (fasta formatted) sequences to the de novo database and build the OTU table**

#### for zOTUs

	usearch10 -otutab demultiplexed_seqs/demultiplexed_seqs_merged.fq -zotus rep_set_zotus_filt_relabeled.fa -otutabout zotutab.txt -mapout zmap.txt

#### for OTUs
	
	usearch10 -otutab demultiplexed_seqs/demultiplexed_seqs_merged.fq -otus rep_set97_filt.fa -otutabout otutab97.txt -mapout map.txt

**NOTE:** **The rest of the pipeline will be the same for zOTUs or OTUs, I've just done zOTUs below to demonstrate. USEARCH manual recommends making tables for both to compare, even if you are sure you will use only one or the other method downstream.**

**###Step 6: Add taxonomic classifications**

The goal of this step is to provide taxonomic classifications for each OTU. We do this using the RDP classifier with the GreenGenes database.

	biom convert -i zotutab.txt -o zotutab.biom --table-type 'OTU table' --to-json

	assign_taxonomy.py -m rdp -i rep_set_zotus_filt_relabeled.fa -o rdp_assigned_taxonomy_zotus -c 0.5 -t <taxonomy database filepath> -r <rep set filepath> --rdp_max_memory 10000

	biom add-metadata -i zotutab.biom --observation-metadata-fp rdp_assigned_taxonomy_zotus/rep_set_zotus_filt_relabeled_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy -o zotutab_wTax.biom
	
Here are the different databases to use for taxonomic classification:

	16S Greengenes 99% cutoff
	taxonomy (-t): /db_files/gg_files/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt
	rep set (-r): /db_files/gg_files/gg_13_8_otus/rep_set/99_otus.fasta 
	
	16S Greengenes 97% cutoff
	taxonomy (-t): /db_files/gg_files/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
	rep set (-r): /db_files/gg_files/gg_13_8_otus/rep_set/97_otus.fasta 

	18S SILVA 99% cutoff
	taxonomy (-t): /db_files/silva_files/Silva119_release/taxonomy_eukaryotes/97/taxonomy_97_7_levels_18S.txt
	rep set (-r): /db_files/silva_files/Silva119_release/rep_set_eukaryotes/97/Silva_119_rep_set97_18S.fna
	
	18S Protist Ribosomal Reference (PR2)  
	taxonomy (-t):	/db_files/pr2_files/2016_02/entire_database/qiime_gb203_taxo.txt
	rep set (-r): /db_files/pr2_files/2016_02/entire_database/mothur_qiime_gb203.fasta
	
	ITS UNITE 97% cutoff
	/db_files/UNITE_files/sh_qiime_release_s_02.03.2015/sh_refs_qiime_ver7_97_s_02.03.2015.fasta
	/db_files/UNITE_files/sh_qiime_release_s_02.03.2015/sh_taxonomy_qiime_ver7_97_s_02.03.2015.txt

**###Step 7: Remove chloroplasts/mitochondria**

It is a good idea to check for and remove chloroplast and mitochondria sequences in most sample types:

	filter_taxa_from_otu_table.py -i zotutab_wTax.biom -o zotutab_wTax_noChloroMito.biom -n c__Chloroplast,f__mitochondria

**###Step 8: Odds and ends**

Convert OTU table back to text format if desired:

	biom convert -i zotutab_wTax_noChloroMito.biom -o zotutab_wTax_noChloroMito.txt --to-tsv --header-key=taxonomy

Get stats on OTU table - # OTUs, # samples, etc:

	biom summarize-table -i zotutab_wTax_noChloroMito.txt -o zotutab_wTax_noChloroMito_smry.txt
	usearch10 -otutab_stats zotutab_wTax_noChloroMito.txt -output report.txt
