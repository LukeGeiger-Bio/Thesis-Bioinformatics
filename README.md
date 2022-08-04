# Thesis-Bioinformatics
2021-2023 Princeton Senior Thesis Project; RNA-seq scripts, DESeq2 code for R, Orthofinder, WGCNA, scRNA-seq(?)

Authored by Luke Geiger, supported by (and much obliged to) Ian Traniello, advised by Dr. Sarah Kocher and Dr. Catherine Peña. Funded through the Fred Fox Class of 1939 Fund, the Princeton University Dept. of Ecology and Ev. Biology, and the Princeton Neuroscience Institute's Lambert Award.

# File Metadata
## All_Metadata.txt
Comma-delimited metadata file for all runs, not grouped by individual animal. 
## DESeq2_Pena_2020_Social_Isolation
RStudio code for combining runs by animal, generating a cts and coldata file for DESeq2, running DESeq2, and select graphical/data outputs.
## RNA_Seq_Accession_Metadata
Global metadata file for all runs and the experiment in general.
## RNA-seq_Pipeline
SBASH job scripts for running pre-DESeq2 RNA-seq analysis steps on a computing cluster, from .FASTA to count matrices. 
## cts_sum_m
Processed counts matrix; download and re-upload into R environment to remove factor levels brought on by DESeq2 steps on a computing cluster. 
## GeneOntology_All.R
R code for performing directional pathway enrichment analysis on mouse, fly, bumblebee, and finch samples. Dotplot, ridgeplot, barplot. 

# Data Availability
## Murine Data
GSE146472, Adolescent Social Isolation Reprograms the Medial Amygdala: Transcriptome and Sex Differences in Reward
Deena M. Walker, Xianxiao Zhou, Aarthi Ramakrishnan, Hannah M. Cates, Ashley M. Cunningham, Catherine J. Peña, Rosemary C. Bagot, Orna Issler, Yentl Van der Zee, Andrew P. Lipschultz, Arthur Godino, Caleb J. Browne, Georgia E. Hodes, Eric M. Parise, Angélica Torres-Berrio, Pamela J. Kennedy, Li Shen, Bin Zhang, Eric J. Nestler
bioRxiv 2020.02.18.955187; doi: https://doi.org/10.1101/2020.02.18.955187

# Tool Citations
## DESeq2:

Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

Tutorial URL: 
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#acknowledgments

## ffq

Gálvez-Merchán, Á., et al. (2022). Metadata retrieval from sequence databases with ffq. bioRxiv 2022.05.18.492548.

Github URL:
https://github.com/pachterlab/ffq

## fastQC

Simon Andrews

Github URL:
https://github.com/s-andrews/FastQC

## STAR aligner

Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.

Github URL:
https://github.com/alexdobin/STAR

## featureCounts/subread package

Yang Liao, Gordon K. Smyth, Wei Shi, featureCounts: an efficient general purpose program for assigning sequence reads to genomic features, Bioinformatics, Volume 30, Issue 7, 1 April 2014, Pages 923–930, https://doi.org/10.1093/bioinformatics/btt656

http://subread.sourceforge.net/
