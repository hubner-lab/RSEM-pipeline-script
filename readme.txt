Introduction:

rsem-script-5.0.sh is a bash script that automates the RSEM pipeline.
it can be used for estimating gene and isoform expression levels,
given RNA-Seq data obtained from organisms with a reference genome and annotation.


Software Required:

1. RSEM
2. STAR
3. fastp


Usage Instructions:
(if you feel lost in the Linux command line environment, please check out this friendly tutorial:
https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)

First,
prepare the following dirs:
1. Genome data dir
   contains:
   1.1. a reference genome FASTA file
        # file extension: *.fa or *.fasta
   1.2. an annotation GTF file
        # file extension: *.gtf
2. RNA-Seq data dir
   contains:
   2.1. RNA-Seq FASTQ files
        # file extension: *.fastq or *.fastq.gz
	# for paired-end sequencing data - filenames should be marked with 'R1' or 'R2' in one of the following formats:
	  1. sample_name_R?.fastq ('R?' is between '_' and '.')
	  2. sample_name.R?.fastq ('R?' is between '.' and '.')
          # as '?' can be '1' or '2'

for example:

human
└── LSLNGS2015
    ├── GENOME_data
    │   ├── Homo_sapiens.GRCh38.86.gtf
    │   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
    └── RNASEQ_data
        ├── GM12878.rep1.R1.fastq.gz
        ├── GM12878.rep1.R2.fastq.gz
        ├── GM12878.rep2.R1.fastq.gz
        ├── GM12878.rep2.R2.fastq.gz
        ├── K562.rep1.R1.fastq.gz
        ├── K562.rep1.R2.fastq.gz
        ├── K562.rep2.R1.fastq.gz
        ├── K562.rep2.R2.fastq.gz

Second,
run the execution command specified below.


Execution Command:

nohup bash [path to rsem-script-5.0.sh] [path to genome data dir] [path to RNA-Seq data dir] [no. of threads] > [output filename] &

# 'nohup'  means: the process won't be killed in case the connection to the server fails.
# 'bash'   means: run the script specified next.
# '>'      means: print outputs to the file specified next (through this file you'll be able to follow the pipeline's activity).
# '&'      means: run the script in the background (allows you to keep working with the command line while the pipeline is running).

for example:

nohup \
bash /home/roy/scripts/rsem-script-5.0.sh \
/home/roy/data/human/LSLNGS2015/GENOME_data \
/home/roy/data/human/LSLNGS2015/RNASEQ_data \
10 \
> output.txt \
&

# please notice that the script is expecting 3 arguments to run properly:
  1. a path to the reference genome data dir 
  2. a path to the RNA-Seq data dir
  3. no. of threads for multithread proccesing (if you're not sure - enter '10')



Results:

you'll find the results in your RNA-Seq data dir, under the prefix: '/rsem-pipeline-output-data/rsem_' for each sample.

