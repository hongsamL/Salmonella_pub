
#!/bin/bash
# PIPELINE.sh [folder_location] [threads]

INPUT_LOC=$1
NUM_THREADS=$2

checkrequirements (){
	# this function checks if the required packages are installed, if not exits the terminal
	command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "bowtie2 not installed.  Aborting."; exit 1; }
	command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools not installed.  Aborting."; exit 1; }
	command -v varscan >/dev/null 2>&1 || { echo >&2 "varscan not installed.  Aborting."; exit 1; }
	command -v parallel >/dev/null 2>&1 || { echo >&2 "parallel not installed.  Aborting."; exit 1; }
	command -v snp-sites >/dev/null 2>&1 || { echo >&2 "snp-sites not installed.  Aborting."; exit 1; }
	python -c 'import pkgutil; print(1 if pkgutil.find_loader("pandas") else 0)'| grep 1 >/dev/null 2>&1 || { echo >&2 "python package pandas not installed.  Aborting."; exit 1;}
	python -c 'import pkgutil; print(1 if pkgutil.find_loader("Bio.SeqIO") else 0)'| grep 1 >/dev/null 2>&1 || { echo >&2 "biopython not installed.  Aborting."; exit 1;}
	python -c 'import pkgutil; print(1 if pkgutil.find_loader("pyfaidx") else 0)'| grep 1 >/dev/null 2>&1 || { echo >&2 "python package pyfaidx not installed.  Aborting."; exit 1;}
	if [ ! -d "$INPUT_LOC"/inputs ]; then
  		echo "No 'inputs' folder "
		exit 1;
	fi
	if [ ! -d "$INPUT_LOC"/ref ]; then
  		echo "NO 'ref' FOLDER"
		exit 1;
	fi
}

cd $INPUT_LOC
INPUT_LOC=$(pwd)

#----------------Check that reference exists--------
#!/bin/bash
checkRef(){
	numfasta="$(ls *.fasta|wc -l)"
	if [ "$numfasta" -gt 1 ]; then
		echo "multiple fasta files detected\n reference must be in a single .fasta file"
		exit 1
	fi
	if [ "$numfasta" -eq 0 ]; then
		echo "NO fasta file detected\n Reference must be in a single .fasta file"
		exit 1
	fi
	
}
# ---------------Prepping the reference-------------
prepRef (){
	# index reference fasta
	cd ref
	checkRef
	bowtie2-build -f *.fa* reference
	cd ..
}

# ------------Check FASTQ reads exist in ----------
checkFASTQ(){
	#check 2 fastqs in folder
	numfastq="$(ls *.fastq*|wc -l)"
	if [ "$numfastq" -gt 2 ]; then
		echo "multiple fastq files detected\n There can only be 2 .fastq or fastq.gz files in this folder"
		exit 1
	fi
	if [ "$numfastq" -eq 0 ]; then
		echo "NO fasta fastq files found"
		exit 1
	fi
	#check naming of fastq
	if ls *_R1* 1> /dev/null 2>&1; then
    		echo "Read 1 found"
	else
    		echo "Read 1 not found. File name must contain '_R1'"
		exit 1
	fi
	if ls *_R2* 1> /dev/null 2>&1; then
    		echo "Read 2 found"
	else 
    		echo "Read 2 not found. File name must contain '_R2'"
		exit 1
	fi
}

# ------------Aligning reads to reference----------

alignRef (){
	# align each fastq to reference fasta
	REF=$INPUT_LOC/"ref/reference"
	mkdir BAM
	ls inputs | while read dir;do
		mkdir BAM/$dir
		cd inputs/$dir
		checkFASTQ
		read1="$(ls *_R1*)"
		read2="$(ls *_R2*)"
		sampleNAME="$(ls|cut -f1 -d '_'|sort -u)"
		bowtie2 -p $NUM_THREADS -x $REF -1 $read1 -2 $read2 |samtools view -bS - | samtools sort - > ../../BAM/$dir/$sampleNAME.bam
		cd ../..
	done
}

# --------------Alignment metrics-----------------
alignMetrics (){
	# compute alignment metrics for each bam file
	mkdir $INPUT_LOC/align-metrics
	echo > stats_parallel.txt
	ls BAM | while read dir;do
		bam="$(ls BAM/$dir | cut -f1 -d '.')"
		echo "samtools stats BAM/$dir/$bam*.bam | grep ^SN | cut -f 2- >> align-metrics/$bam.stats.tsv" >> stats_parallel.txt
		echo "samtools depth BAM/$dir/$bam*.bam |  awk '{sum+=\$3} END { print \"Average Coverage:"\\t"\" , sum/NR}' >> align-metrics/$bam.stats.tsv" >> stats_parallel.txt
	done
	cat stats_parallel.txt | parallel -j $NUM_THREADS
	# Reformat alignment metrics
	ls align-metrics| while read l;do
		for i in `seq 1 32`;do
			s=$(echo $l | cut -f1 -d '.')
			echo $s >> align-metrics/$l.temp;
		done
		paste align-metrics/$l.temp align-metrics/$l >> align-metrics/$l.new
	done
	rm $INPUT_LOC/align-metrics/*.temp
	cat align-metrics/*.new > alignment_metrics
	python alignment_stats.py -f alignment_metrics -o metrics
}
# ---------------SNP CALLING-----------------------
snpCall() {
	# snp calling with Varscan
	mkdir $INPUT_LOC/snps
	ls BAM | while read dir;do
		sampleNAME="$(ls BAM/$dir|cut -f1 -d '.'|sort -u)"
		echo "samtools mpileup -f $INPUT_LOC/ref/ref.fasta BAM/$dir/$sampleNAME.sorted.bam | \
	varscan mpileup2snp - --min-coverage 8 --min-avg-qual 30 --min-var-freq 0.90 --output-vcf > snps/$sampleNAME.vcf" >> varscan_commands
	done
	cat varscan_commands | parallel -j $NUM_THREADS
}
# ---------------SNP FILTERING----------------------
snpFilter() {
	mkdir filtered_snps
	ls snps | while read file;do
		sampleName="$(echo $file | cut -f1 -d '.')"
		echo $sampleName
		head -n25 snps/$file > snps/$sampleName.vcfhead
		tail -n+25 snps/$file > snps/$sampleName.nohead
		varscan limit snps/$sampleName.nohead --regions-file $INPUT_LOC/filter_regions/*.txt --not-file snps/$sampleName.snps.filtered
		cat snps/$sampleName.vcfhead snps/$sampleName.snps.filtered > snps/$sampleName.snps.filtered.vcf
		rm snps/$sampleName.vcfhead
		rm snps/$sampleName.nohead
		rm snps/$sampleName.snps.filtered
		mv snps/$sampleName.snps.filtered.vcf filtered_snps
	done
}
# --------------Consensus and SNP matrix--------------------
consensus() {
	mkdir consensus
	ls filtered_snps | while read f; do
		python makeconsensus.py -r ref/ref.fasta -f filtered_snps/$f -o consensus ;
	done
	mkdir snpmatrix
	cat consensus/* > consensus/aligned.fasta
	snp-sites -m -o snps.fasta consensus/aligned.fasta
	mv snps.fasta snpmatrix
}

#----------------Pairwise SNPs--------------------------
snpCounts() {
	python pairwise_snps.py -f $INPUT_LOC/snpmatrix/snps.fasta
}

#--------------MAIN -------------------------------------
main() {
	checkrequirements
	echo "preparin thre reference \n"
	prepRef
	echo "algning to reference \n"
	alignRef
	echo "calculating alignment metrics \n"
	#alignMetrics
	echo "initiating SNP calling \n"
	#snpCall
	echo "filtering SNPs \n"
	snpFilter
	echo "creating consesnsus sequences \n"
	consensus
	echo "calculating SNP differences"
	snpCounts
}

main
