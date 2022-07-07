#! /bin/bash

# rsem-script-5.0.sh

Readme()
{
cat << EOF


Introduction:

rsem-script-5.0.sh is a bash script that automates the RSEM pipeline, and meant to simplify things as much as possible.
it can be used for estimating gene and isoform expression levels, given RNA-Seq data obtained from organisms with a reference genome and annotation.


Software Required:

1. RSEM   (https://github.com/deweylab/RSEM)
2. STAR   (https://github.com/alexdobin/STAR)
3. fastp  (https://github.com/OpenGene/fastp)


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


EOF
}


# arguments first check:
if [ $# -ne 3 ]
then
	echo
	echo "$# arguments received, but the pipeline's script needs 3 to operate properly. printing readme file and exiting..."
	Readme
	exit 0
fi


# read arguments:
genome_data_path=$1
RNASeq_data_path=$2
n_threads=$3


# getting reference genome fasta file name:
reference_file_name=$(find $genome_data_path -maxdepth 1 -name *.fa -printf "%f\n")
length=${#reference_file_name}
if [ $length -eq 0 ]
then
	reference_file_name=$(find $genome_data_path -maxdepth 1 -name *.fasta -printf "%f\n")
	length=${#reference_file_name}
	if [ $length -eq 0 ]
	then
		echo "couldn't find fasta file in $genome_data_path. exiting..."
		exit 0
	fi
fi
reference_name=${reference_file_name:0:7}


# getting reference annotation gtf/gff file name:
annotation_file_name=$(find $genome_data_path -maxdepth 1 -name *.gtf -printf "%f\n")
annotation_type='gtf'
length=${#annotation_file_name}
if [ $length -eq 0 ]
then
	annotation_file_name=$(find $genome_data_path -maxdepth 1 -name *.gff -printf "%f\n")
	annotation_type='gff3'
	length=${#annotation_file_name}
	if [ $length -eq 0 ]
	then
		echo "couldn't find gtf/gff file in $genome_data_path. exiting..."
		exit 0
	fi
	echo
	echo "-------------------------------------------------------------------------------------------"
	echo "please note:"
	echo
	echo "couldn't find a standard GTF annotation file. the pipeline will try to convert from gff to gtf, but errors may occur..."
	echo "if RSEM fails because the BAM file declares more reference sequences than RSEM knows (specified in the rsem.log file),
	echo "please try to use a standard GTF annotation file.""
	echo "-------------------------------------------------------------------------------------------"
	echo
fi


# getting RNA-Seq fastq file names:
fastq_file_names=( $(ls $RNASeq_data_path | grep .fastq) )
length=${#fastq_file_names[@]}
if [ $length -eq 0 ]
then
	echo "couldn't find fastq files in $RNASeq_data_path. exiting..."
	exit 0
fi


# n_threads input checks:
if ((n_threads < 1 || n_threads > 1000))
then
	echo "no. of threads received: $n_threads. allowed values: between '1' and '1000'. if you're not sure, please enter '10'. exiting..."
	exit 0
fi


# extracting paired/single end option, and sample names:
is_paired_end=0
sample_names=()
for (( i=0; i<$length; i++ ));
do
	name=${fastq_file_names[i]}
	if [[ "$name" == *".R2."* ]] || [[ "$name" == *"_R2."* ]]
	then
		is_paired_end=1
	else
		# remove suffix and save sample name:
		name=${name%".fastq"*}
		name=${name%"_R1"}
		name=${name%".R1"}
		sample_names+=($name)
	fi
done

n_samples=${#sample_names[@]}


# paired/single end conclusion check:
if [ $is_paired_end -eq 1 ]
then
	let double_samples=2*${#sample_names[@]}
	if [ ${#fastq_file_names[@]} -ne $double_samples ]
	then
		echo "no. of fastq filenames do not match no. of sample names for paired-end seq data (1:2 ratio expected). exiting..."
		exit 0
	fi
else
	if [ ${#fastq_file_names[@]} -ne ${#sample_names[@]} ]
	then
		echo "no. of fastq filenames do not match no. of sample names for single-end seq data (1:1 ratio expected). exiting..."
		exit 0
	fi
fi



echo rsem-script-5.0.sh
date +%d-%m-%y
date +%T
echo

cat << EOF

Configurations:

genome data path:
    $genome_data_path
RNASeq data path:
    $RNASeq_data_path
reference file name:
    $reference_file_name
reference name:
    $reference_name
reference annotation file name:
    $annotation_file_name
annotation type:
    $annotation_type
EOF

echo "fastq file names:"
printf '    %s\n' ${fastq_file_names[@]}

if [ $is_paired_end -eq 1 ]
then
	end="paired-end"
else
	end="single-read"
fi

echo "sample names ($n_samples, $end):"
printf '    %s\n' ${sample_names[@]}

echo "aligner:"
echo "    STAR"
echo "no. of threads: $n_threads"
echo



# Create mapping indices:

# mapping indices for STAR:
: '
usage:
STAR --runMode genomeGenerate \
     --genomeDir path_to_genomedir \
     --genomeFastaFiles reference_fasta_file(s)
'
echo "Creating mapping indices for STAR in $genome_data_path/star..."
mkdir $genome_data_path/star
STAR --runThreadN $n_threads \
     --runMode genomeGenerate \
     --genomeDir $genome_data_path/star \
     --genomeFastaFiles $genome_data_path/$reference_file_name \
     --sjdbGTFfile $genome_data_path/$annotation_file_name
echo "Creating mapping indices for STAR is done."
echo


# mapping indices for RSEM:
: '
usage:
rsem-prepare-reference [options] reference_fasta_file(s) reference_name
'
echo "Creating mapping indices for RSEM in $genome_data_path/rsem..."
mkdir $genome_data_path/rsem
rsem-prepare-reference \
	--$annotation_type \
	$genome_data_path/$annotation_file_name \
	$genome_data_path/$reference_file_name \
	$genome_data_path/rsem/$reference_name
echo "Creating mapping indices for RSEM is done."
echo


# RNA-Seq analysis starts:
echo "RNA-Seq data analysis starts..."
echo
mkdir $RNASeq_data_path/rsem-pipeline-output-data


# Trim with fastp:

echo "Trimming and producing quality reports with fastp..."
date +"%T"
mkdir $RNASeq_data_path/rsem-pipeline-output-data/trimmed
mkdir $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports

if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		r1_fastq_file=${fastq_file_names[2*i]}
		r2_fastq_file=${fastq_file_names[2*i+1]}
		echo
		echo "trimming sample $sample_name"
		date +"%T"
		echo "input files: $r1_fastq_file and $r2_fastq_file"
		echo
		fastp -i $RNASeq_data_path/$r1_fastq_file \
		      -I $RNASeq_data_path/$r2_fastq_file \
		      -o $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.R1.P.fq.gz \
		      -O $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.R2.P.fq.gz \
		      -w 16

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports/$sample_name
		mv fastp.html fastp.json $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports/$sample_name
	done
else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		fastq_file=${fastq_file_names[i]}
		echo
		echo "trimming sample $sample_name"
		date +"%T"
		echo "input file: $fastq_file"
		echo
		fastp -i $RNASeq_data_path/$fastq_file \
		      -o $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.fq.gz \
		      -w 16

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports/$sample_name
		mv fastp.html fastp.json $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports/$sample_name
	done
fi
echo
echo "Trimming with fastp is done."
date +"%T"
echo "trimmed files are here: $RNASeq_data_path/rsem-pipeline-output-data/trimmed"
echo "quality reports are here: $RNASeq_data_path/rsem-pipeline-output-data/trimmed/quality-reports"
echo


# Alignment:

# Mapping with STAR:

echo "Mapping with STAR..."
if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo
		echo "mapping sample $sample_name"
		date +"%T"

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name
		STAR --genomeDir $genome_data_path/star \
		     --readFilesCommand zcat \
		     --readFilesIn $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.R1.P.fq.gz \
		                   $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.R2.P.fq.gz \
		     --outSAMtype BAM SortedByCoordinate \
		     --limitBAMsortRAM 16000000000 \
		     --outSAMunmapped Within \
		     --twopassMode Basic \
		     --outFilterMultimapNmax 1 \
		     --quantMode TranscriptomeSAM \
		     --runThreadN $n_threads \
		     --outFileNamePrefix "$RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name/"
	done
else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo
		echo "mapping sample $sample_name"
		date +"%T"

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name
		STAR --genomeDir $genome_data_path/star \
		     --readFilesCommand zcat \
		     --readFilesIn $RNASeq_data_path/rsem-pipeline-output-data/trimmed/$sample_name.fq.gz \
		     --outSAMtype BAM SortedByCoordinate \
		     --limitBAMsortRAM 16000000000 \
		     --outSAMunmapped Within \
		     --twopassMode Basic \
		     --outFilterMultimapNmax 1 \
		     --quantMode TranscriptomeSAM \
		     --runThreadN $n_threads \
		     --outFileNamePrefix "$RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name/"
	done
fi
echo
echo "Mapping with STAR is done."
echo


# Quantification with RSEM:

: '
usage:
rsem-calculate-expression [optioins] upstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --sam/--bam [--paired-end] input reference_name sample_name
'

echo "Quantifying gene expression with RSEM..."
echo
if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo "quantificating sample $sample_name"
		date +"%T"

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name
		rsem-calculate-expression \
			--bam \
			--no-bam-output \
			-p $n_threads \
			--paired-end \
			--forward-prob 0 \
			$RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name/Aligned.toTranscriptome.out.bam \
			$genome_data_path/rsem/$reference_name \
			$RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name/rsem >& \
			$RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name/rsem.log

		echo "done."
		echo "the resaults are here: $RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name"
		echo
	done

else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo "quantificating sample $sample_name"
		date +"%T"

		mkdir $RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name

		rsem-calculate-expression \
			--bam \
			--no-bam-output \
			-p $n_threads \
			--forward-prob 0 \
			$RNASeq_data_path/rsem-pipeline-output-data/star_$sample_name/Aligned.toTranscriptome.out.bam \
			$genome_data_path/rsem/$reference_name \
			$RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name/rsem >& \
			$RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name/rsem.log

		echo "done."
		echo "the resaults are here: $RNASeq_data_path/rsem-pipeline-output-data/rsem_$sample_name"
		echo ""
	done
fi
echo
echo "quantification with RSEM is done."
date +"%T"
echo
echo "you'll find the results in your RNA-Seq data dir, under the prefix: '/rsem-pipeline-output-data/rsem_' for each sample."



