This project describes the bioinformatic reproductible workflow associated to the paper "Multiple *Wolbachia* subpopulations co-occur in single *Culex pipiens* mosquitoes".

# Anvi'o metagenomic workflows and binning to reconstruct Wolbachia genomes

## Download data 

To explore the *Wolbachia* populations in single *Culex pipiens* midguts, you can download the six midgut metagenomes stored in ENA with the accession number PRJEB56379.

The sample names and correspoding ENA accession IDs are reported here:
| Sample name | ENA sample accession |
|:-----------:|:--------------------:|
| M03         | ERS13507996          |
| M07         | ERS13507997          |
| M11         | ERS13507998          |
| M12         | ERS13507999          |
| M01         | ERS13509473          |
| M09         | ERS13509474          |


The fastq.gz file names and corresponding ENA accession IDs are reported here: 

| Sample name | R1                                               | R2                                              | Read ENA accession |
|:-----------:|:------------------------------------------------:|:-----------------------------------------------:|:------------------:|
| M03         | CBX_AEOSDA_1_1_HT2J2BBXX.12BA150_clean.fastq.gz  | CBX_AEOSDA_1_2_HT2J2BBXX.12BA150_clean.fastq.gz | ERR10300695        |
| M03         | CBX_AEOSDA_2_1_HT2J2BBXX.12BA150_clean.fastq.gz  | CBX_AEOSDA_2_2_HT2J2BBXX.12BA150_clean.fastq.gz | ERR10300696        |
| M07         | CBX_AFOSDA_1_1_HT2J2BBXX.12BA151_clean.fastq.gz  | CBX_AFOSDA_1_2_HT2J2BBXX.12BA151_clean.fastq.gz | ERR10300697        |
| M07         | CBX_AFOSDA_2_1_HT2J2BBXX.12BA151_clean.fastq.gz  | CBX_AFOSDA_2_2_HT2J2BBXX.12BA151_clean.fastq.gz | ERR10300698        |
| M11         | CBX_AGOSDA_1_1_HT2J2BBXX.12BA152_clean.fastq.gz  | CBX_AGOSDA_1_2_HT2J2BBXX.12BA152_clean.fastq.gz | ERR10300699        |
| M11         | CBX_AGOSDA_2_1_HT2J2BBXX.12BA152_clean.fastq.gz  | CBX_AGOSDA_2_2_HT2J2BBXX.12BA152_clean.fastq.gz | ERR10300700        |
| M12         | CBX_AHOSDA_1_1_HT2J2BBXX.12BA153_clean.fastq.gz  | CBX_AHOSDA_1_2_HT2J2BBXX.12BA153_clean.fastq.gz | ERR10300701        |
| M12         | CBX_AHOSDA_2_1_HT2J2BBXX.12BA153_clean.fastq.gz  | CBX_AHOSDA_2_2_HT2J2BBXX.12BA153_clean.fastq.gz | ERR10300702        |
| M01         | CBX_AIOSDA_2_1_HT2J2BBXX.12BA154_clean.fastq.gz  | CBX_AIOSDA_2_2_HT2J2BBXX.12BA154_clean.fastq.gz | ERR10300703        |
| M09         | CBX_AJOSDA_2_1_HT2J2BBXX.12BA155_clean.fastq.gz  | CBX_AJOSDA_2_2_HT2J2BBXX.12BA155_clean.fastq.gz | ERR10300704        |


## Anvi'o installation

You can install anvi'o following https://anvio.org/install/. 
We used the anvio-7 version for this study. 

## Configuration files 

### Sample files 
[samples_midguts.txt](files/metadata/samples_midguts.txt)

[samples_ovaries.txt](files/metadata/samples_ovaries.txt)

### Config files
[megahit_midguts.txt](files/metadata/megahit_midguts.json)

[megahit_ovaries.txt](files/metadata/megahit_ovaries.json)

### Script files
[run_metagenomics_midguts.sh](files/script/run_metagenomics_midguts.sh)

[run_metagenomics_ovaries.sh](files/script/run_metagenomics_ovaries.sh)

## Run workflow

Overall, the anvi'o metagenomics workflow includes the following steps: 

1) Quality control of metagenomic short reads using illumina-utils.
2) Assembly of quality-filtered short reads using MEGAHIT.
3) Taxonomical (GTDB) and functionnal (COG, KOfam) assignment of contigs using anvi'o anvi-scg-taxonomy, anvi-run-kegg-kofams and anvi-run-ncbi-cogs programs
4) Reads recruitment from all samples (all-against-all) on each sample assembly using Bowtie2.
5) Profile the BAM files by relate them to the contigs databases using the anvi-profile program. The program computes the coverage by nucleotide position, the variants (at nucleotide, codon and amino-acid levels) and the structural variants (indel or insertions).
6) Merge the resulting single anvi'o profile databases.

Command lines to run workflow (on SGE cluster): 

```
# midguts
qsub -cwd -V -N qsub_midguts -q infinit.q -pe thread 2 -R y megahit_midguts.sh

# ovaries
qsub -cwd -V -N qsub_ovaries -q infinit.q -pe thread 2 -R y megahit_midguts.sh
```

The complete anvi'o metagenomics workflow results on the following folders: 
* `00_LOGS`: Log files for each processus
* `01_QC`: Quality-filtered short metagenomic reads and final statistics.
* `02_FASTA`: Megahit assembly of each sample (FASTA files).
* `03_CONTIGS`: Anvi’o contigs databases for each assembly with taxonomic and functionnal annoations
* `04_MAPPING`: BAM files from Bowtie2 read recruitment on each assembly (all-against-all)
* `05_ANVIO_PROFILE`: Anvi’o single profiles for each sample
* `06_MERGED`: Anvi’o merged profile databases for each assembly

All the anvi'o merged profile databases from midgut and ovary samples are available at [10.5281/zenodo.7183277](10.5281/zenodo.7183277).


## Assembly

We used the anvi-display-contig-stats to generate a report including all our samples.

```
anvi-display-contigs-stats 03_CONTIGS/*-contigs.db --report-as-text --output-file CONTIGS-STATS.txt 
```

Output for midguts:

  contigs_db	          | Culex_M01 |Culex_M03  |Culex_M07  | Culex_M09 | Culex_M11 | Culex_M12 | 
|:---------------------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
  Total Length	        | 360227583 | 249300854 | 344361265 |150118498  | 613485459 | 3365464   |
  Num Contigs	          | 175034	  | 133816	  | 169207	  | 92326	    | 237602	  | 2230      |
  Num Contigs > 100 kb	| 0	        | 0	        | 0	        | 0	        | 0	        | 0         |
  Num Contigs > 50 kb	  | 0	        | 0	        | 0	        | 0	        | 4	        | 0         |
  Num Contigs > 20 kb	  | 6	        | 3	        | 7	        | 0	        | 49	      | 0         |
  Num Contigs > 10 kb	  | 277 	    | 89	      | 282	      | 11	      | 1569	    | 1         |
  Num Contigs > 5 kb	  | 5667	    | 2281	    | 5201	    | 603	      | 20463	    | 17        |
  Num Contigs > 2.5 kb	| 41027	    | 22940	    | 38340	    | 8605	    | 88045	    | 145       |
  Longest Contig	      | 28329	    | 48513	    | 36009	    | 15727	    | 75699	    | 16613     |
  Shortest Contig	      | 1000	    | 1000	    | 1000	    | 1000	    | 1000	    | 1000      |
  Num Genes (prodigal)	| 265946	  | 194792	  | 256186	  | 123803	  | 441013	  | 2683      |
  L50	                  | 51752	    | 42561	    | 50264	    | 32193	    | 63529	    | 779       |
  L75	                  | 101140	  | 80672	    | 98167	    | 58676	    | 127361	  | 1443      |
  L90	                  | 141809    |	110425    |	137369    |	77946     |	184732    |	1902      |
  N50	                  | 2233	    | 1942	    | 2200	    | 1619	    | 3043	    | 1431      |
  N75	                  | 1505	    | 1391	    | 1487	    | 1254	    | 1906	    | 1153      |
  N90	                  | 1178	    | 1139	    | 1173  	  | 1092	    | 1343	    | 1053      |
  Archaea_76	          | 95	      | 93	      | 115	      | 91	      | 187	      | 0         |
  Bacteria_71	          | 74	      | 74	      | 85	      | 69	      | 200	      | 1         |
  Protista_83	          | 79	      | 79	      | 81	      | 66	      | 99	      | 1         |
  Ribosomal_RNA_12S	    | 0	        | 0	        | 0	        | 0	        | 0	        | 0         |
  Ribosomal_RNA_16S	    | 0	        | 0	        | 0	        | 0	        | 3	        | 0         |
  Ribosomal_RNA_18S	    | 1	        | 1	        | 1	        | 1	        | 1	        | 2         |
  Ribosomal_RNA_23S	    | 0	        | 0	        | 0	        | 0	        | 5	        | 0         |
  Ribosomal_RNA_28S	    | 1	        | 1	        | 1	        | 1	        | 2	        | 1         |
  Ribosomal_RNA_5S	    | 0	        | 0	        | 0	        | 0	        | 0	        | 0         |
  Transfer_RNAs	        | 5989	    | 3682	    | 5582	    | 1794	    | 13694     | 50        |
  eukarya (Protista_83)	| 1	        | 1	        | 1	        | 1	        | 1	        | 0         |
  archaea (Archaea_76)	| 1	        | 1	        | 1	        | 1	        | 0	        | 0         |
  bacteria (Bacteria_71)|	1	        | 1	        | 1	        | 0	        | 3	        | 0         |


Output for ovaries: 

  contigs_db	          | Culex_O03 |	Culex_O07	| Culex_O11	| Culex_O12 |
|:---------------------:|:---------:|:---------:|:---------:|:---------:|
  Total Length	        | 618074637	| 238540242	| 599730856	| 550583510 |
  Num Contigs	          | 247594	  | 131080	  | 243215	  | 236095    |
  Num Contigs > 100 kb	| 0	        | 0	        | 0	        | 0         |
  Num Contigs > 50 kb	  | 3	        | 3	        | 3	        | 3         |
  Num Contigs > 20 kb	  | 38	      | 22	      | 45	      | 42        |
  Num Contigs > 10 kb	  | 1264	    | 118	      | 1186    	| 849       |
  Num Contigs > 5 kb	  | 18599	    | 1535	    | 17386	    | 13707     |
  Num Contigs > 2.5 kb	| 87141	    | 19639	    | 83577	    | 72556     |
  Longest Contig	      | 87423	    | 87270	    | 75610	    | 87423     |
  Shortest Contig	      | 1000	    | 1000	    | 1000	    | 1000      |
  Num Genes (prodigal)	| 445444	  | 186211	  | 432072	  | 395881    |
  L50	                  | 67596	    | 43369	    | 66639	    | 65780     |
  L75	                  | 134695	  | 80292	    | 132760	  | 130997    |
  L90	                  | 193940	  | 108729	  | 190943	  | 187157    |
  N50	                  | 2895	    | 1890	    | 2845	    | 2636      |
  N75	                  | 1847	    | 1389	    | 1816	    | 1706      |
  N90                   | 1323	    | 1142	    | 1312	    | 1266      |
  Archaea_76	          | 159	      | 125	      | 161	      | 161       |
  Bacteria_71	          | 164	      | 136	      | 165	      | 154       |
  Protista_83	          | 98	      | 79	      | 93	      | 86        |
  Ribosomal_RNA_12S	    | 0	        | 0	        | 0	        | 0         |
  Ribosomal_RNA_16S	    | 1	        | 1	        | 1	        | 1         |
  Ribosomal_RNA_18S	    | 1	        | 1	        | 1	        | 1         |
  Ribosomal_RNA_23S	    | 3	        | 1	        | 1	        | 1         |
  Ribosomal_RNA_28S	    | 2	        | 1	        | 1	        | 1         |
  Ribosomal_RNA_5S	    | 0	        | 0	        | 0	        | 0         |
  Transfer_RNAs	        | 14242	    | 1681	    | 13443	    | 11832     |
  bacteria (Bacteria_71)|	1	        | 2	        | 2	        | 1         |
  eukarya (Protista_83)	| 1	        | 1	        | 1	        | 1         |
  archaea (Archaea_76)	| 0	        | 1	        | 0	        | 1         |


## Binning

We used the anvi-cluster-contigs anvi'o program (with CONCOCT) to reconstruct bacterial genomes. We tried some different numbers of cluster for each sample to make the binning.

### Midguts

Overall, we reconstructed a bacterial genome in only 1 midgut sample (M11) using a number of cluster = 5 for CONCOCT.

```
anvi-cluster-contigs -p 06_MERGED/Culex_M11/PROFILE.db \
                     -c 03_CONTIGS/M11-contigs.db \
                     -C CONCOCT_5 \
                     --driver CONCOCT \
                     --clusters 5 \
                     --just-do-it
```

We used the anvi-estimate-genome-completeness to estimate the completeness and redundancy of reconstructed genomes based on single-copy core genes.

```
anvi-estimate-genome-completeness -c 03_CONTIGS/Culex_M11-contigs.db -p 06_MERGED/Culex_M11/PROFILE.db -C CONCOCT_5
```

| bin_name | domain   | confidence | % completion | % redundancy | num_splits | total length |
|:--------:|:--------:|:----------:|:------------:|:------------:|:----------:|:------------:|
| Bin_4    | BLANK    | 0.9        | 0            | 0            | 76,194     | 142,877,516  |
| Bin_3    | BACTERIA | 0.9        | 91.55        | 5.63         | 53,089     | 103,399,384  |
| Bin_2    | EUKARYA  | 0.5        | 43.37        | 6.02         | 64,968     | 286,105,977  |
| Bin_1    | BLANK    | 0.7        | 0            | 0            | 16,323     | 25,781,103   |
| Bin_0    | EUKARYA  | 0.4        | 42.17        | 16.87        | 26,863     | 55,150,479   |


### Ovaries

Overall, we reconstructed a bacterial genome in all ovary samples using different number of clusters for each sample.

#### O03
```
anvi-estimate-genome-completeness -c 03_CONTIGS/Culex_O03-contigs.db -p 06_MERGED/Culex_O03/PROFILE.db -C CONCOCT_6
```

| bin_name | domain   | confidence | % completion | % redundancy | num_splits | total length |
|:--------:|:--------:|:----------:|:------------:|:------------:|:----------:|:------------:|
| Bin_1    | BACTERIA | 0.8        | 92.96        | 7.04         | 64,324     | 100,851,143  |
| Bin_0    | EUKARYA  | 0.6        | 57.83        | 22.89        | 30,387     | 67,553,370   |
| Bin_3    | BLANK    | 0.6        | 0            | 0            | 69,569     | 185,522,201  |
| Bin_2    | ARCHAEA  | 0.3        | 38.16        | 11.84        | 36,463     | 184,851,012  |
| Bin_5    | BLANK    | 1          | 0            | 0            | 35,485     | 60,816,468   |
| Bin_4    | BLANK    | 1          | 0            | 0            | 11,165     | 18,273,443   |

#### O07

```
anvi-estimate-genome-completeness -c 03_CONTIGS/Culex_O07-contigs.db -p 06_MERGED/Culex_O07/PROFILE.db -C CONCOCT_6
```
| bin_name | domain   | confidence | % completion | % redundancy | num_splits | total length |
|:--------:|:--------:|:----------:|:------------:|:------------:|:----------:|:------------:|
| Bin_2    | BLANK    | 0.9        | 0            | 0            | 37,977     | 55,880,364   |
| Bin_5    | BACTERIA | 1          | 92.96        | 4.23         | 24,216     | 37,601,680   |
| Bin_1    | EUKARYA  | 0.6        | 65.06        | 10.84        | 15,400     | 38,222,279   |
| Bin_3    | BLANK    | 0.9        | 0            | 0            | 32,903     | 78,923,862   |
| Bin_4    | BLANK    | 1          | 0            | 0            | 3,988      | 5,927,922    |
| Bin_0    | BLANK    | 1          | 0            | 0            | 16,401     | 21,783,135   |

#### O11
```
anvi-estimate-genome-completeness -c 03_CONTIGS/Culex_O11-contigs.db -p 06_MERGED/Culex_O11/PROFILE.db -C CONCOCT_8
```

| bin_name | domain   | confidence | % completion | % redundancy | num_splits | total length |
|:--------:|:--------:|:----------:|:------------:|:------------:|:----------:|:------------:|
| Bin_3    | BLANK    | 0.9        | 0            | 0            | 42,901     | 128,367,800  |
| Bin_1    | BACTERIA | 0.9        | 91.55        | 5.63         | 30,194     | 47,193,763   |
| Bin_6    | BLANK    | 1          | 0            | 0            | 6,732      | 10,992,557   |
| Bin_4    | BLANK    | 1          | 0            | 0            | 31,073     | 51,503,590   |
| Bin_5    | BLANK    | 0.6        | 0            | 0            | 52,147     | 101,954,135  |
| Bin_2    | EUKARYA  | 0.6        | 57.83        | 16.87        | 15,879     | 40,566,635   |
| Bin_7    | BLANK    | 1          | 0            | 0            | 30,092     | 48,678,430   |
| Bin_0    | ARCHAEA  | 0.3        | 42.11        | 10.53        | 34,004     | 170,275,946  |


#### O12
```
anvi-estimate-genome-completeness -c 03_CONTIGS/Culex_O12-contigs.db -p 06_MERGED/Culex_O12/PROFILE.db -C CONCOCT_5
```

| bin_name | domain   | confidence | % completion | % redundancy | num_splits | total length |
|:--------:|:--------:|:----------:|:------------:|:------------:|:----------:|:------------:|
| Bin_2    | ARCHAEA  | 0.2        | 30.26        | 14.47        | 57,616     | 222,612,692  |
| Bin_4    | BACTERIA | 0.9        | 91.55        | 15.49        | 72,612     | 122,214,342  |
| Bin_0    | EUKARYA  | 0.6        | 57.83        | 20.48        | 20,919     | 60,123,079   |
| Bin_1    | BLANK    | 0.9        | 0            | 0            | 20,391     | 31,423,540   |
| Bin_3    | BLANK    | 1          | 0            | 0            | 64,359     | 114,004,857  |


## Taxonomy 

All the bacterial reconstructed genomes were assigned to *Wolbachia* (GTDB). 

## Refinement

We used the anvi-refine program to manually refine our reconstructed genomes. To do that, we removed all the splits that do not correspond to *Wolbachia*.

After refinement using the anvi-refine program, we obtained five *Wolbachia* MAGs:

| Wolbachia MAG | % completion | % redundancy | num_splits | total length | % GC content |
|:-------------:|:------------:|:------------:|:----------:|:------------:|:------------:|
| M11           | 91.55        | 0            | 138        | 1,331,260    | 34.2         |
| O11           | 91.55        | 0            | 119        | 1,290,070    | 34.2         |
| O03           | 91.55        | 0            | 73         | 1,164,954    | 33.8         |
| O07           | 91.55        | 0            | 143        | 1,340,038    | 34.4         |
| O12           | 91.55        | 0            | 75         | 1,181,440    | 33.9         |

The FASTA files for individual bin are available at [10.5281/zenodo.7183304](10.5281/zenodo.7183304).


# Explore the co-occurence of Wolbachia subpopulations accross metagenomes using the references-mode

## Configuration files

### Fasta file
[ref_mode_fasta.txt](files/metadata/ref_mode_fasta.txt)

### Sample file
[ref_mode_samples.txt](files/metadata/ref_mode_samples.txt)

### Config file
[ref_mode_config.json](files/metadata/ref_mode_config.json)

## Run the anvi'o metagenomics workflow using "references-mode"

The references-mode includes the same steps described above but it also allows to map all the metagenomes (midguts and ovaries) on each *Wolbachia* MAG we reconstructed. This mapping is essential to study the metapangenomics and the variability among *Wolbachia* subpopulations in our metagenomes. 

To run it:
[run_metagenomics_ref_mode.sh](files/script/run_metagenomics_ref_mode.sh)

The final anvi'o merged profile.db from references-mode are available at [10.5281/zenodo.7183324](10.5281/zenodo.7183324).


## Variability and (Meta)pangenomics

To study the variability and metapangenomics, you can run these scripts:

* [1-run-pangenome.sh](files/script/snv_metapan/1-run-pangenome.sh): Prepare *Wolbachia* MAGs contigs.db, Wolbachia reference genomes and run pangenomes.
* [2-run-blast.sh](files/script/snv_metapan/2-run-blast.sh): Download the MLST + wsp gene sequences and blast them against *Wolbachia* MAGs contigs.db. 
* [3-generate-coverage-data-for-GC.Rmd](files/script/snv_metapan/3-generate-coverage-data-for-GC.Rmd): Generate coverage values of each metagenome and MLST+wsp genes at Gene Clusters (GCs) level.
* [4-add-coverage-data-to-pan.sh](files/script/snv_metapan/4-add-coverage-data-to-pan.sh): Add coverage values previously generated in pangenomes. 
* [5-generate-variability-tables-and-interactive-snv.Rmd](files/script/snv_metapan/5-generate-variability-tables-and-interactive-snv.Rmd): Generate variability tables (Single Nucleotide Variants, Single Codon Variants, Single Amino Acid Variants), indels/insertion tables at intra-MAGs and inter-MAGs level and prepare interactive view of SNVs.
* [6-prepare-additional-data-to-snvs.R](files/script/snv_metapan/6-prepare-additional-data-to-snvs.R): Prepare additionnal data (MLST and wsp assignments) to add in interactive view of SNVs. 
* [7-select-SNV-SCV-SAAVs.R](files/script/snv_metapan/7-select-SNV-SCV-SAAVs.R): Select SCV and SAAV generated from only the filtered SNVs.
* [8-merge_gene_cluster_SNV.R](files/script/snv_metapan/8-merge_gene_cluster_SNV.R): Identify shared SNVs, Shared Single copy Genes (SSG) from the Gene Cluster tables (pangenomics), link the SSG and SNVs informations to export it in SNV interactive databases.
* [9-add-SSG-to-interactive-SNV.sh](files/script/snv_metapan/9-add-SSG-to-interactive-SNV.sh): Add SSG infos in snv interactive databases.
* [10-SNV-figure-for-selected-splits.Rmd](files/script/snv_metapan/10-SNV-figure-for-selected-splits.Rmd): Generate SNV and gene figures for selected splits from *Wolbachia* MAG O11. We then used Inskape to create the final SNV figures for these splits. 


