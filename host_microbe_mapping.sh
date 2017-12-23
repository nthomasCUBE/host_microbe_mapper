#!/usr/bin/bash

# --------------------------------------------------------------------------

echo "-------------------------------------------------------------------"
echo
echo "               HOST_MICROBE_MAPPING v0.0.6                        "
echo "                  24 MAY 2017                                      "
echo
echo "-------------------------------------------------------------------"

module load htseq
module load cutadapt
module load trimmomatic
module load sortmerna

function contains() {
    local n=$#
    local value=${!n}
    for ((i=1;i < $#;i++)) {
        if [ "${!i}" == "${value}" ]; then
            echo "y"
            return 0
        fi
    }
    echo "n"
    return 1
}

HOST_MAPPER="tophat"

while getopts "F:R:P:H:I:O:A:B:C:X:J:K:L:U" opt; do
  case $opt in
    F)
      FWD_ORIG=$OPTARG
      ;;
    R)
      RVS_ORIG=$OPTARG
      ;;
    P)
      GENOME=$OPTARG
      ;;
    H)
      GENOME2=$OPTARG
      ;;
    I)
      GENOME3=$OPTARG
      ;;
    O)
      HOST_MAPPER=$OPTARG
      ;;
    A)
      EUKARYOTE_INDEX_1=$OPTARG
      ;;
    B)
      EUKARYOTE_INDEX_2=$OPTARG
      ;;
    C)
      PROCS=$OPTARG
      ;;
    X)
      OUTPUT=$OPTARG
      ;;
    J)
      GFF=$OPTARG
      ;;
    K)
      GFF2=$OPTARG
      ;;
    L)
      GFF3=$OPTARG
      ;;
    U)
      PREPROCESSING=TRUE
      ;;
  esac
done

if [ -z ${PROCS} ]
then
        echo "Amount of CPUs (option -C) not defined"
        exit
else
        echo "NUMBER OF CPUs: ${PROCS}"
fi

if [ -z ${FWD_ORIG} ]
then
	echo "Forward reads (option -F) not defined"
	exit
else
        echo "FORWARD READS: ${FWD_ORIG}"
fi

if [ -z ${RVS_ORIG} ]
then
        echo "Reverse reads (option -R) not defined"
        #exit
else
        echo "REVERSE READS: ${RVS_ORIG}"
fi

if [ -z ${OUTPUT} ]
then
        echo "Output dir (option -X) not defined"
        exit
else
        echo "Output dir: ${OUTPUT}"
fi


if [ -z $GENOME ]
then
	echo "No Pathogen genome defined"
else
        echo "Pathogen genome: $GENOME"

fi

if [ -z $GENOME2 ]
then
	echo "No Eukaryotic genome defined"
else
        echo "Eukaryotic genome: $GENOME2"
fi

if [ -z $GENOME3 ]
then
	echo "No second eukaryotic genome defined"
else
        echo "Eukaryotic genome: $GENOME3"
fi

OPTIONS=("tophat" "star")
if [ $(contains "${OPTIONS[@]}" $HOST_MAPPER) == "y" ];
then
    echo "The mapping tools '$HOST_MAPPER' was selected for mapping RNA-seq in the host genomes"
else
    echo "$HOST_MAPPER (option -O) must be either tophat or start"
fi

echo
echo
echo
# --------------------------------------------------------------------------

DIR=$(mktemp -d)

echo "WORKING directory: $DIR"

FWD=$DIR/FWD.fastq
cp ${FWD_ORIG} ${FWD}
if [ -z ${RVS_ORIG} ]
then
	echo ""
else
	RVS=$DIR/RVS.fastq
	cp ${RVS_ORIG} ${RVS}
fi

if [ -z "$GENOME" ]
then
	echo ""
else
	cp $GENOME $DIR/GENOME
	GENOME=$DIR/GENOME
fi

if [ -z "$GENOME2" ]
then
	echo ""
else
        cp $GENOME2 $DIR/GENOME2
        GENOME2=$DIR/GENOME2
fi

if [ -z "$GENOME3" ]
then
	echo ""
else
        cp $GENOME3 $DIR/GENOME3
        GENOME3=$DIR/GENOME3
fi

# --------------------------------------------------------------------------

rm -f $DIR/make_test_runs.log
touch $DIR/make_test_runs.log

mkdir $DIR/results
#mkdir $DIR/orig
# --------------------------------------------------------------------------

echo "step-0: fastqc"

fastqc $FWD -o $DIR &>> $DIR/make_test_runs.log
#cp $FWD $DIR/orig/FWD_orig.fastq
if [ -z ${RVS_ORIG} ]
then
	echo "No reverse reads"
else
	fastqc $RVS -o $DIR &>> $DIR/make_test_runs.log
#	cp $RVS $DIR/orig/RVS_orig.fastq
fi

mv $DIR/*fastqc.html $DIR/results
rm $DIR/*fastqc.zip
fastq-stats $FWD | grep "^reads" | sed 's/reads/INFO::fastq\treads/g' &> $DIR/stats.log
# --------------------------------------------------------------------------


if [ -z ${PREPROCESSING} ]
then
	echo "step-1: cutadapt"
	if [ -z ${RVS_ORIG} ]
	then
	        cutadapt \
	                    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	                    -o ${FWD}.trimmed.fastq ${FWD} &>> $DIR/make_test_runs.log
        	mv ${FWD}.trimmed.fastq ${FWD}
        	cutadapt \
        	            -a AAAAAAAAAAAAAAA \
        	            -o ${FWD}.trimmed.fastq ${FWD} &>> $DIR/make_test_runs.log
        	mv ${FWD}.trimmed.fastq ${FWD}
        	cutadapt \
        	            -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        	            -o ${FWD}.trimmed.fastq ${FWD} &>> $DIR/make_test_runs.log
        	mv ${FWD}.trimmed.fastq ${FWD}
        	cutadapt \
        	            -a TTTTTTTTTTTTTTT \
        	            -o ${FWD}.trimmed.fastq ${FWD} &>> $DIR/make_test_runs.log
        	mv ${FWD}.trimmed.fastq ${FWD}
	else
		echo "=paired=",${FWD},${RVS}
		cutadapt \
		            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
		            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		            -o ${FWD}.trimmed.fastq -p ${RVS}.trimmed.fastq ${FWD} ${RVS} &>> $DIR/make_test_runs.log
		mv ${FWD}.trimmed.fastq ${FWD}
		mv ${RVS}.trimmed.fastq ${RVS}
		cutadapt \
		            -a AAAAAAAAAAAAAAA \
		            -A TTTTTTTTTTTTTTT \
		            -o ${FWD}.trimmed.fastq -p ${RVS}.trimmed.fastq ${FWD} ${RVS} &>> $DIR/make_test_runs.log
		mv ${FWD}.trimmed.fastq ${FWD}
		mv ${RVS}.trimmed.fastq ${RVS}
	fi
	fastq-stats $FWD | grep "^reads" | sed 's/reads/INFO::cutadapt\treads/g' &>> $DIR/stats.log

	# --------------------------------------------------------------------------

	echo "step-2: trimmomatic"

	if [ -z ${RVS_ORIG} ]
	then
	        trimmomatic SE -threads ${PROCS} -phred33 ${FWD} $DIR/${ID}_forward_paired.fq.gz LEADING:8 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50 &>> $DIR/make_test_runs.log
	        zcat $DIR/${ID}_forward_paired.fq.gz > ${FWD}
	else
		ID="trimmomatic"
		trimmomatic PE -threads ${PROCS} -phred33 ${FWD} ${RVS} $DIR/${ID}_forward_paired.fq.gz $DIR/${ID}_forward_unpaired.fq.gz $DIR/${ID}_reverse_paired.fq.gz $DIR/${ID}_reverse_unpaired.fq.gz LEADING:8 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50 &>> $DIR/make_test_runs.log
		zcat $DIR/${ID}_forward_paired.fq.gz > ${FWD}
		zcat $DIR/${ID}_reverse_paired.fq.gz > ${RVS}
		rm $DIR/trimmomatic_*gz
	fi
	fastq-stats $FWD | grep "^reads" | sed 's/reads/INFO::trimming\treads/g' &>> $DIR/stats.log
	# --------------------------------------------------------------------------

	echo "step-3: sortMeRNA"

	FASTQ_FWD=$FWD
	FASTQ_REVERSE=$RVS
	FASTQ_MERGED=${DIR}/merged.fastq
	FASTQ_FQ_OUT1=${FWD}.no_rRNA.fastq
	FASTQ_FQ_OUT2=${RVS}.no_rRNA.fastq

	smr=$(which sortmerna)
	smr=$(dirname $smr)/../

	bash merge-paired-reads.sh ${FASTQ_FWD} ${FASTQ_REVERSE} ${FASTQ_MERGED} &>> $DIR/make_test_runs.log
	if [ -z ${RVS_ORIG} ]
	then
	        sortmerna -a $PROCS --ref ${smr}/rRNA_databases/silva-bac-16s-id90.fasta,${smr}/index/silva-bac-16s-id90.db:\
${smr}/rRNA_databases/silva-bac-23s-id98.fasta,${smr}/index/silva-bac-23s-id98.db:\
${smr}/rRNA_databases/silva-arc-16s-id95.fasta,${smr}/index/silva-arc-16s-id95.db:\
${smr}/rRNA_databases/silva-arc-23s-id98.fasta,${smr}/index/silva-arc-23s-id98.db:\
${smr}/rRNA_databases/silva-euk-18s-id95.fasta,${smr}/index/silva-euk-18s-id95.db:\
${smr}/rRNA_databases/silva-euk-28s-id98.fasta,${smr}/index/silva-euk-28s-id98.db:\
${smr}/rRNA_databases/rfam-5s-database-id98.fasta,${smr}/index/rfam-5s-database-id98.db:\
${smr}/rRNA_databases/rfam-5.8s-database-id98.fasta,${smr}/index/rfam-5.8s-database-id98.db --reads ${FASTQ_FWD} --aligned ${FASTQ_MERGED}_accepted --log --sam --fastx --other ${FASTQ_MERGED}_rejected &>> $DIR/make_test_runs.log
	        less ${FASTQ_MERGED}_rejected* | grep -v "^$" > ${FASTQ_MERGED}_rejected.fastq.fixed
	        cp ${FASTQ_MERGED}_rejected.fastq.fixed ${FASTQ_FQ_OUT1}
	else
		sortmerna -a $PROCS --ref ${smr}/rRNA_databases/silva-bac-16s-id90.fasta,${smr}/index/silva-bac-16s-id90.db:\
${smr}/rRNA_databases/silva-bac-23s-id98.fasta,${smr}/index/silva-bac-23s-id98.db:\
${smr}/rRNA_databases/silva-arc-16s-id95.fasta,${smr}/index/silva-arc-16s-id95.db:\
${smr}/rRNA_databases/silva-arc-23s-id98.fasta,${smr}/index/silva-arc-23s-id98.db:\
${smr}/rRNA_databases/silva-euk-18s-id95.fasta,${smr}/index/silva-euk-18s-id95.db:\
${smr}/rRNA_databases/silva-euk-28s-id98.fasta,${smr}/index/silva-euk-28s-id98.db:\
${smr}/rRNA_databases/rfam-5s-database-id98.fasta,${smr}/index/rfam-5s-database-id98.db:\
${smr}/rRNA_databases/rfam-5.8s-database-id98.fasta,${smr}/index/rfam-5.8s-database-id98.db --reads ${FASTQ_MERGED} --paired_in --num_alignments 1 --aligned ${FASTQ_MERGED}_accepted --log --sam --fastx --other ${FASTQ_MERGED}_rejected &>> $DIR/make_test_runs.log
	        less ${FASTQ_MERGED}_rejected* | grep -v "^$" > ${FASTQ_MERGED}_rejected.fastq.fixed
		bash unmerge-paired-reads.sh ${FASTQ_MERGED}_rejected.fastq.fixed ${FASTQ_FQ_OUT1} ${FASTQ_FQ_OUT2} &>> $DIR/make_test_runs.log
	fi
	cp ${FASTQ_MERGED}_accepted.sam ${DIR}/results/rRNA_matches.sam
	cp ${FASTQ_MERGED}_accepted*log ${DIR}/results/SortMeRNA_log.txt
	mv ${FWD}.no_rRNA.fastq ${FWD}
	fastq-stats $FWD | grep "^reads" | sed 's/reads/INFO::sortmerna\treads/g' &>> $DIR/stats.log

	if [ -z ${RVS_ORIG} ]
	then
		echo "no reverse reads"
	else
		mv ${RVS}.no_rRNA.fastq ${RVS}
	fi
	rm ${FASTQ_MERGED}*accepted*
	rm ${FASTQ_MERGED}*rejected*
	#rm ${FASTQ_MERGED}
	cd ${cur_path}
else
	echo "no preprocessing will be performed..."
fi

# --------------------------------------------------------------------------

echo "step-4: bwa (mapping to bacterium)"

if [ -z "$GENOME" ]
then
	echo "INFO: No mapping on the pathogen"
else
	bwa index $GENOME &>> $DIR/make_test_runs.log

	if [ -z ${RVS_ORIG} ]
	then
	        bwa aln -t ${PROCS} $GENOME $FWD > ${DIR}/fwd.sai 2>> $DIR/make_test_runs.log
	        bwa samse $GENOME $DIR/fwd.sai $FWD > $DIR/both.sam 2>> $DIR/make_test_runs.log
	else
		bwa aln -t ${PROCS} $GENOME $FWD > $DIR/fwd.sai 2>> $DIR/make_test_runs.log
		bwa aln -t ${PROCS} $GENOME $RVS > $DIR/rev.sai  2>> $DIR/make_test_runs.log
		bwa sampe $GENOME $DIR/fwd.sai $DIR/rev.sai $FWD $RVS > $DIR/both.sam  2>> $DIR/make_test_runs.log
	fi
	samtools view --threads ${PROCS} -b -h $DIR/both.sam > $DIR/both_sort.bam
	samtools sort --threads ${PROCS} $DIR/both_sort.bam > $DIR/both.bam  2>> $DIR/make_test_runs.log
	rm $DIR/both_sort.bam
	samtools flagstat $DIR/both.bam
	cp $DIR/both.bam $DIR/results/pathogen_mapping.bam
	samtools view --threads ${PROCS} -h -f 0x4 ${DIR}/both.bam > ${DIR}/unmapped1.bam
	samtools sort --threads ${PROCS} -n ${DIR}/unmapped1.bam > ${DIR}/unmapped1.bam.sorted
	bedtools bamtofastq -i ${DIR}/unmapped1.bam.sorted -fq ${FWD} -fq2 ${RVS}
	echo "GENOME=${GENOME}"
	echo "DIR=$DIR"
	rm ${GENOME}.*
	rm $DIR/fwd.sai $DIR/rev.sai $DIR/both.sam $DIR/both.bam ${DIR}/unmapped1.*
	samtools flagstat ${DIR}/results/pathogen_mapping.bam | grep -w "properly" | cut -f 1 -d " " |  awk '{ print $1/2 }' | sed 's/^/INFO::pathogen\treads\t/g' &>> $DIR/stats.log

fi

# --------------------------------------------------------------------------

echo "step-5: $HOST_MAPPER (mapping to Host)"

if [ $HOST_MAPPER == 'tophat' ]
then
	if [ -z ${GENOME2} ];
	then
		echo "step-5: not run, no host genome"
	else
		ID=tophat
                if [ -z ${EUKARYOTE_INDEX_1} ]
                then
 			bowtie2-build $GENOME2 ${GENOME2} &>> ${DIR}/make_test_runs.log
		else
			bname=`basename ${EUKARYOTE_INDEX_1}`
			bdir=`dirname ${EUKARYOTE_INDEX_1}`
			echo "bname=${bname}"
			for i in ${bname}*bt2
			do
				new_name=`echo $i | sed "s/${bname}//"`
				new_name=${DIR}/GENOME2${new_name}
			        ln -s "$bdir/$i" $new_name
			done
		fi
		tophat -p $PROCS -o $DIR/${ID} ${GENOME2} ${FWD} ${RVS} &>> $DIR/make_test_runs.log
		tail -n 50 $DIR/${ID}/*align*
		zip ${DIR}/eukaryotic_host1_mapped.bam.zip ${DIR}/${ID}/accepted*bam
		samtools sort --threads ${PROCS} -n ${DIR}/${ID}/unmapped.bam > ${DIR}/${ID}/unmapped.bam.sorted
		rm ${GENOME2}.*
		samtools flagstat ${DIR}/${ID}/accepted*bam | grep -w "properly" | cut -f 1 -d " " |  awk '{ print $1/2 }' | sed 's/^/INFO::host1\treads\t/g' &>> $DIR/stats.log
	fi
else
        if [ -z ${GENOME2} ];
	then
                echo "step-5: not run, no host genome"
	else
		if [ -z "$EUKARYOTE_INDEX_1" ]
		then
			echo "INFO: Making new index: ${EUKARYOTE_INDEX_1}"
		        mkdir ${DIR}/Genome2_Index
			STAR --runThreadN ${PROCS} --runMode genomeGenerate --genomeDir ${DIR}/Genome2_Index --genomeFastaFiles ${GENOME2}
		else
			ln -s ${EUKARYOTE_INDEX_1} ${DIR}/Genome2_Index
			echo "INFO: Using existing index: ${EUKARYOTE_INDEX_1}"
		fi
		if [ -z ${RVS_ORIG} ];
		then
                        STAR --genomeDir ${DIR}/Genome2_Index --outSAMstrandField intronMotif --readFilesIn ${FWD} --outFileNamePrefix ${DIR}/star_ --runThreadN ${PROCS}
		else
			STAR --genomeDir ${DIR}/Genome2_Index --outSAMstrandField intronMotif --readFilesIn ${FWD} ${RVS} --outFileNamePrefix ${DIR}/star_ --runThreadN ${PROCS}
		fi
		samtools view --threads ${PROCS} -h -F 524 ${DIR}/star_Aligned.out.sam > ${DIR}/star_Aligned.out.sam.ProperPairs
		samtools sort --threads ${PROCS} -n ${DIR}/star_Aligned.out.sam.ProperPairs > ${DIR}/star_Aligned.out.sam.ProperPairs.bam
		cp $DIR/star_Aligned.out.sam.ProperPairs.bam $DIR/results/host1_mapping.bam
		rm ${DIR}/star_Aligned*
		rm -r ${DIR}/Genome2_Index
	        rm ${DIR}/star_Log*out ${DIR}/star_SJ.out.tab
                samtools flagstat ${DIR}/results/host1_mapping.bam | grep -w "properly" | cut -f 1 -d " " | awk '{ print $1/2 }' |  sed 's/^/INFO::host1\treads\t/g' &>> $DIR/stats.log

	fi
fi

# --------------------------------------------------------------------------

echo "step-6: $HOST_MAPPER (mapping to Host-2)"

if [ $HOST_MAPPER == 'tophat' ]
then
	if [ -z ${GENOME3} ];
	then
	        echo "step-6: not run, no host genome"
	else
		ID=tophat2
		if [ -z ${EUKARYOTE_INDEX_2} ]
		then
			bowtie2-build $GENOME3 ${GENOME3} &>> make_test_runs.log
		else
                        bname=`basename ${EUKARYOTE_INDEX_2}`
                        bdir=`dirname ${EUKARYOTE_INDEX_2}`
                        echo "bname=${bname}"
                        for i in ${bname}*bt2
                        do
                                new_name=`echo $i | sed "s/${bname}//"`
                                new_name=${DIR}/GENOME3${new_name}
                                ln -s "$bdir/$i" $new_name
                        done
		fi
		tophat -p $PROCS -o $DIR/${ID} ${GENOME3} ${FWD} ${RVS} &>> $DIR/make_test_runs.log
                zip ${DIR}/eukaryotic_host2_mapped.bam.zip ${DIR}/${ID}/accepted*bam
		tail -n 50 $DIR/${ID}/*align*
		samtools sort --threads ${PROCS} -n ${DIR}/${ID}/unmapped.bam > ${DIR}/${ID}/unmapped.bam.sorted
		rm ${GENOME3}.*
                samtools flagstat ${DIR}/${ID}/accepted*bam | grep -w "properly" | cut -f 1 -d " " |  awk '{ print $1/2 }' | sed 's/^/INFO::host2\treads\t/g' &>> $DIR/stats.log
	fi
else
	if [ -z ${GENOME3} ];
	then
                echo "step-6: not run, no host genome"
	else
		if [ -z ${EUKARYOTE_INDEX_2} ]
		then
		        mkdir ${DIR}/Genome3_Index
	                echo "INFO: Making new index: ${EUKARYOTE_INDEX_2}"
		        STAR --runThreadN ${PROCS} --runMode genomeGenerate --genomeDir ${DIR}/Genome3_Index --genomeFastaFiles ${GENOME3}
		else
        	        ln -s ${EUKARYOTE_INDEX_2} ${DIR}/Genome3_Index
			echo "INFO: Using existing index ${EUKARYOTE_INDEX_2}"
		fi
                if [ -z ${RVS_ORIG} ];
                        STAR --genomeDir ${DIR}/Genome3_Index --outSAMstrandField intronMotif --readFilesIn ${FWD} --outFileNamePrefix ${DIR}/star_ --runThreadN ${PROCS}
                then
		        STAR --genomeDir ${DIR}/Genome3_Index --outSAMstrandField intronMotif --readFilesIn ${FWD} ${RVS} --outFileNamePrefix ${DIR}/star_ --runThreadN ${PROCS}
		fi
	        samtools view --threads ${PROCS} -h -F 524 ${DIR}/star_Aligned.out.sam > ${DIR}/star_Aligned.out.sam.ProperPairs
	        samtools sort --threads ${PROCS} -n ${DIR}/star_Aligned.out.sam.ProperPairs > ${DIR}/star_Aligned.out.sam.ProperPairs.bam
                cp $DIR/star_Aligned.out.sam.ProperPairs.bam $DIR/results/host2_mapping.bam
		rm -r ${DIR}/Genome3_Index
		rm ${DIR}/star_Log*out ${DIR}/star_SJ.out.tab ${DIR}/star_Aligned*
                samtools flagstat ${DIR}/results/host2_mapping.bam | grep -w "properly" | cut -f 1 -d " " |  awk '{ print $1/2 }' | sed 's/^/INFO::host2\treads\t/g' &>> $DIR/stats.log
	fi
fi

# --------------------------------------------------------------------------

if [ -z $GFF ]
then
        echo "No count extraction for - pathogen"
else
	samtools view --threads ${PROCS} -hS ${DIR}/results/pathogen_mapping.bam > ${DIR}/results/pathogen_mapping.sam
	htseq-count ${DIR}/results/pathogen_mapping.sam ${GFF} -t gene > ${DIR}/results/pathogen_mapping_counts.txt
	rm ${DIR}/results/pathogen_mapping.sam
fi

if [ -z $GFF2 ]
then
        echo "No count extraction for - host1"
else
        samtools view --threads ${PROCS} -hS ${DIR}/results/host1_mapping.bam > ${DIR}/results/host1_mapping.sam
        htseq-count ${DIR}/results/host1_mapping.sam ${GFF2} -t gene > ${DIR}/results/host1_mapping_counts.txt
	for type in mRNA processed_transcript NMD_transcript_variant aberrant_processed_transcript exon
	do
	        htseq-count ${DIR}/results/host1_mapping.sam ${GFF2} -t ${type} > ${DIR}/results/host1_mapping_counts_${type}.txt
	done
        rm ${DIR}/results/host1_mapping.sam
fi

if [ -z $GFF3 ]
then
        echo "No count extraction for - host2"
else
        samtools view --threads ${PROCS} -hS ${DIR}/results/host2_mapping.bam > ${DIR}/results/host2_mapping.sam
        htseq-count ${DIR}/results/host2_mapping.sam ${GFF3} -t gene > ${DIR}/results/host2_mapping_counts.txt
        for type in mRNA processed_transcript NMD_transcript_variant aberrant_processed_transcript exon
        do
                htseq-count ${DIR}/results/host2_mapping.sam ${GFF3} -t ${type} > ${DIR}/results/host2_mapping_counts_${type}.txt
        done
        rm ${DIR}/results/host2_mapping.sam
fi

# --------------------------------------------------------------------------

echo "INFO: Calculations are stored under: $DIR"

rm ${DIR}/GENOME*
rm ${FASTQ_MERGED}
rm ${DIR}/results/*bam
rm ${DIR}/results/*sam
rm ${DIR}/FWD.fastq
rm ${DIR}/RVS.fastq

pwd
echo "OUTPUT=",${OUTPUT}
mv ${DIR} ${OUTPUT}

