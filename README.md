#METAGENOMIC WORKFLOW
#FIRSTDAY

#Installing conda
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --set auto_activate_base false

#Add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

##Never install in base environment
conda create -n nf-core-dsl1 -c bioconda nexflow=22.10.6
conda activate nf-core-dsl1

#Find hpc info
sinfo -o “%P %c %m %t %l”
P = Partition
c = cpu
m = memory
t = state
l = time limit

#symlink : shared document without copy

ln -s $PATH/FILE > $PATH/FILE

#Nextflow need singularity, loading module is needed in HPC
spack load apptainer
or
module load singularity

#Running Eager for QC step
screen -R eager #new screen
screen -r eager #resume screen

#Eager Script for adapter removal
nextflow run nf-core/eager \
-r 2.5.0 \
-c envadna.conf \
-profile singularity \
--input '/shared/data/samplesheets/eager_input.tsv' \
--outdir 'eager-screening/' \
--mapper 'bowtie2' \
--fasta '/shared/data/genomes/GCF_000819615.1_ViralProj14015_genomic.fna' \
--run_bam_filtering \
--bam_unmapped_type fastq \
--metagenomic_complexity_filter \
--run_metagenomic_screening \
--metagenomic_tool 'malt' \
--database '/shared/data/databases/malt_db/' \
--run_maltextract \
--maltextract_taxon_list '/shared/data/databases/maixner_mini_targetlist.txt' \
--maltextract_ncbifiles '/shared/data/databases/' \
--skip_preseq \
--skip_deduplication \
--skip_damage_calculation \
--skip_qualimap


#Day2

https://github.com/miwipe/acad_workshop2024
#Additional Information

Filter data:
- Adapter trimming
  AdapterRemoval2
  leeHom (calculate quality phred score > 41 which not supported by downstream metagenomic)
  Cutadapt
  SeqPrep
  Fastp
- Removal of human sequence
  Some package need removal of human DNA (KrakenUniq, MALT)
- low complexity filtering
- Quality trimming
- Duplicate removal can be done before and after Mapping, since metagenomic has strong workflow implying on mapping. Duplicate removal before Mapping is preferable.
- Read length filtering

Mapping:
- Bowtie (Best alignment for metagenomic data, highly taxonomoy data )
- bwa aln
- HOPS

Taxonomyic profiling
- Euka
- Bowtie2+metaDMG(ngsLCA)
- Bowtie2+bam-filter/reassignment+metaDMG(ngsLCA)


#Create Conda env using yaml file
#Make yaml1
name: acad-euks_1
channels:
  - conda-forge
  - bioconda
  - defaults
  - genomewalker
dependencies:
  - python>=3.8,<=3.9
  - Cython>=0.29.24
  - pip
  - pip:
      - pyrle>=0.0.31
      - pyranges>=0.0.112
      - kneed>=0.8.1
      - matplotlib>=3.6.0
  - samtools
  - bowtie2
  - r-base
  - r-ggplot2
  - fastp
  - vsearch
  - sga
  - seqtk
  - fasta-splitter
  - datamash
  - seqkit
#yaml2
name: acad-euks_2
  channels:
    - conda-forge
    - bioconda
    - defaults
  dependencies:
  - beast2=2.6.3=hf1b8bbb_0
  - mafft=7.475=h516909a_0
  - raxml-ng=1.0.1=h7447c1b_0
  - tabview
  - angsd
  - samtools
  - epa-ng
  - gappa
  - seqkit

conda env create -f acad-euks_1.yaml
conda env create -f acad-euks_2.yaml

##BAM FILTER is not installed
conda create -n bamfilter
conda install -c conda-forge -c bioconda -c genomewalker bam-filter



#Analyse data
Quality control of trimmed and merged sequences
When handling large data and mapping against large reference genome collections,it can be important to remove duplicates, to save cpu and run time. For this, I use vsearch https://github.com/torognes/vsearchfast tool that screens for 100% identical sequences (most likely caused by PCR duplication).You can use vsearch --help to familiarize yourself with its options.
vsearch -- Remove the duplicates
Depends on the how big is the fasta in this demo around 10 millions Read
srun -n 1 --mem 1GB  vsearch --fastx_uniques ERR10493277_small-FINAL.fq.gz \
--fastqout ./ERR10493277_small-FINAL.vs.fq --minseqlength 30 --strand both


###Clean out low complex sequences
###Another important aspect is to clean out low complex sequences again several tools can do this,
###I use sga / bbduk in which dust ranges from 1-4 where 1 is the most stringent.
srun -n 1 --mem 5GB sga preprocess -m 30 --dust-threshold=1 ERR10493277_small-FINAL.vs.fq  -o ERR10493277_small-FINAL.vs.d1.fq

# prints the readIDs of each sequence that parsed filters
grep 'M_A00706' ERR10493277_small-FINAL.vs.d1.fq

# we can make this into a text file with all readIDs that parsed
for file in ERR10493277_small-FINAL.vs.d1.fq; do grep 'M_A00706' $file | cut -f2 -d@ > $file.readID.tmp ; done

# we can then inverse grep these readIDs in the original fastq file, to get the readIDs of the filtered sequences
grep 'M_A00706' ERR10493277_small-FINAL.vs.fq | grep -f ERR10493277_small-FINAL.vs.d1.fq.readID.tmp -v | cut -f2 -d@ > ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp

# how many reads did we filter out, and does it correspond to the stdout sga printed?
wc -l ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp

# we can use seqtk to grep out all sequences that where categorized as low complex, what do you see?
seqtk subseq ERR10493277_small-FINAL.vs.fq ERR10493277_small-FINAL.vs.d1.fq.readID_lowcom.tmp


##Extract length distribution
cat ERR10493277_small-FINAL.vs.d1.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.read_length.txt &
