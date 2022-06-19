#!/bin/bash  
#  

cd $PBS_O_WORKDIR

# all module versions correct for artemis. kallisto was changed, RSeQc included in python
module load star/2.5.2a
module load picard/2.7.1
module load samtools/1.6
module load gcc/4.9.3
module load python/3.5.1
module load R/3.3.2-intel
module load rsem/1.3.0
module load ucsc-userapps/348
module load trimmomatic/0.36
module load kallisto/0.43.1
module load fastqc/0.11.3
module load stringtie/1.3.3b 
module load bedtools/2.26.0
module load perl/5.24.0
module load trinity/2.1.1


if [ -z $1 ]; then
	echo "Need to submit name for output files/folders" && exit
fi

if [ -z $2 ]; then
	echo "Need to provide number of processors to be used e.g. 6" && exit
fi

## Define project folder
outDir="/pipeline_output/${1}"
echo $outDir && mkdir -p $outDir

# fastq files stored inside a directory named based on sample IDs
inDir="/pipeline_input/${1}"
echo "inDir= $inDir"

## Specify number of cores
numcores=${2}
echo "numcores= $numcores"


## Annotations
indexDir="/annotations/hg38"
echo "indexDir is $indexDir"

genomeDir=$indexDir"/GRCh38.p10.genome.fa"
echo "genomeDir is $genomeDir"

gtfDir=$indexDir"/gencode.v27.annotation.gtf"
echo "gtf is $gtfDir"

rseqcIndexDir=$indexDir"/hg38_RefSeq.bed"
echo "RSEQC index is $rseqcIndexDir"

genSize=$indexDir"/star/chrNameLength.txt"
echo "GenomeSize is $genSize"

kalIndexDir=$indexDir"/kallisto/kallisto_index_gencode_27.idx"
echo "Kallisto index is $kalIndexDir"

rsemIndexDir=$indexDir"/rsem"
echo "RSEMIndex is $rsemIndexDir"

##Global Internal file structure

fastQCDir=${outDir}"/fastQC/"${1}
echo "FastQCDir is $fastQCDir" && mkdir -p $fastQCDir

fastQC_trim_Dir=${outDir}"/fastQC_trim/"${1}
echo "FastQCDir_trim is $fastQC_trim_Dir" && mkdir -p $fastQC_trim_Dir

trim_Dir=${outDir}"/uniform_trim/${1}"

echo "TrimmedDir is $trim_Dir" && mkdir -p ${trim_Dir}

starDir=${outDir}"/star/"${1}
echo "STARDir is $starDir " && mkdir -p ${starDir}

kalDir=${outDir}"/kallisto/"${1}
echo "KallistoDir is $kalDir" && mkdir -p ${kalDir}

rsemDir=${outDir}"/rsem/"${1}
echo "RSEMDir is $rsemDir" && mkdir -p ${rsemDir}

logDir=${outDir}"/logs"
echo "LogDir is $logDir"

resDir=${outDir}"/splicing/"${1}
echo "SplicingRseDir is $resDir" && mkdir -p ${resDir}

juncDir=${resDir}"/rseqc"
echo "JuncDir is $juncDir" && mkdir -p ${juncDir}

read_len="101"
echo "Readlength is $read_len"

gownDir=${outDir}"/stringtie/"${1} 		
echo "BallgownDir is $gownDir" && mkdir -p $gownDir

TMPDIR=${outDir}"/tmp"
echo "TMPDIR is $TMPDIR" && mkdir -p $TMPDIR
resDir_tmp=$TMPDIR"/${1}_splicing"
echo "resDir tmp is $resDir_tmp" && mkdir -p $resDir_tmp

juncDir_tmp=$resDir_tmp"/rseqc"
echo "juncDir tmp is $juncDir_tmp" && mkdir -p $juncDir_tmp



###Check if files are there

# minFileSize="1M"

minFileSize="1M"
inFile1=${inDir}"/*R1.fq.gz"
echo "File1 is $inFile1"

inFile2=${inDir}"/*R2.fq.gz"
echo "File2 is $inFile2"

find $inFile1 -type f -size $minFileSize -delete
find $inFile2 -type f -size $minFileSize -delete

if [ ! -f $inFile1 ];
then
	echo "Can't input inFile1"
else
	echo "Found" $inFile1
fi

if [ ! -f $inFile2 ];
then
	echo "Can't input inFile2"
else
	echo "Found" $inFile2
fi



##########################################################################################################################
##########################################################################################################################
###FastQC preTrim
##########################################################################################################################
##########################################################################################################################

# # can skip for total/capseq rerun 2017
# if [ ! -f $fastQCDir/${1}_R1.fq_fastqc.zip ];
# then
# 	echo "Running fastQC preTrim"
# 	# could alternativly nest the scp in an elif fragment
# 	if [ ! -f $TMPDIR/${1}*fq.gz && ! -f $TMPDIR/${1}*R2.fq.gz ];
# 		then
# 			scp -r $inDir/*fq.gz $TMPDIR
# 		else 
# 			echo "Already in tmpdir"
# 	fi
 
#  	inFile1=${TMPDIR}"/*R1.fq.gz"
# 	echo "Infile1"$inFile1

# 	inFile2=${TMPDIR}"/*R2.fq.gz"
# 	echo "Infile2"$inFile2
	
# 	fastqc_pretrim_tmp_Dir=$TMPDIR"/${1}_QC"
# 	echo $fastqc_pretrim_tmp_Dir && mkdir -p $fastqc_pretrim_tmp_Dir

# 	echo "fastqc -t $numcores --outdir $fastQCDir $inFile1  $inFile2"
# 	time fastqc \
# 	-t $numcores \
# 	--outdir $fastQCDir \
# 	$inFile1  $inFile2 
# else
# 	echo "Found" $fastQCDir/${1}_R1.fq_fastqc.zip
# fi

# mv $fastqc_pretrim_tmp_Dir $fastQCDir


##########################################################################################################################
##########################################################################################################################
###STAR 
##########################################################################################################################
##########################################################################################################################


star_tmp_Dir=$TMPDIR"/${1}_star"
echo $star_tmp_Dir && mkdir -p $star_tmp_Dir
scp ${starDir}/${1}_star/${1}Aligned.sortedByCoord.out.bam $starDir
scp ${starDir}/${1}_star/${1}Aligned.toTranscriptome.out.bam $starDir
find $starDir/${1}Aligned.sortedByCoord.out.bam -type f -size $minFileSize -delete
find ${starDir}/${1}Aligned.toTranscriptome.out.bam -type f -size $minFileSize -delete

if [ ! -f $starDir/${1}Aligned.sortedByCoord.out.bam ];
then
	echo "Running STAR"

	#if [[ ! -f $TMPDIR/${1}*uniform_R1.fq.gz && ! -f $TMPDIR/${1}*uniform_R2.fq.gz ]];
	if [[ ! -f $TMPDIR/*R1.fq.gz && ! -f $TMPDIR/*R2.fq.gz ]];
		then
			#scp -r $trim_Dir/* $TMPDIR
			scp $inDir/*.fq.gz $TMPDIR
		else 
			echo "Already in tmpdir"
	fi
	
	# inFile1=${TMPDIR}"/*R1_trim.fq.gz"
	#inFile1=${TMPDIR}"/*uniform_R1.fq.gz"
	inFile1=${TMPDIR}/*R1.fq.gz
	echo $inFile1

	# inFile2=${TMPDIR}"/*R2_trim.fq.gz"
	#inFile2=${TMPDIR}"/*uniform_R2.fq.gz"
	inFile2=${TMPDIR}/*R2.fq.gz
	echo $inFile2

	echo "STAR --runMode alignReads \
		--readFilesIn $inFile1 $inFile2 \
		--readFilesCommand zcat \
		--outFileNamePrefix ${star_tmp_Dir}/${1} \
		--genomeDir ${indexDir} \
		--sjdbGTFfile ${gtfDir} \
		--outSJfilterReads Unique \
		--sjdbOverhang 100 \
		--twopassMode Basic \
		--runThreadN ${numcores:=6} \
		--genomeLoad NoSharedMemory \
		--outFilterType BySJout \
		--outFilterMultimapNmax 100 \
		--outFilterMismatchNmax 33 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMismatchNoverLmax 0.3 \
		--outFilterScoreMinOverLread 0.3 \
		--limitOutSJcollapsed 1000000 \
		--limitSjdbInsertNsj 1000000 \
		--alignEndsType EndToEnd \
		--alignSJDBoverhangMin 3  \
		--alignSJoverhangMin 8 \
		--alignIntronMin 20 \
		--winAnchorMultimapNmax 50 \
		--seedSearchStartLmax 12 \
		--chimSegmentMin 20 \
		--outSAMattributes All \
		--outSAMstrandField intronMotif \
		--quantMode TranscriptomeSAM \
		--outSAMattrIHstart 0 \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate"
	time STAR --runMode alignReads \
	--readFilesIn $inFile1 $inFile2 \
	--readFilesCommand zcat \
	--outFileNamePrefix ${star_tmp_Dir}/${1} \
	--genomeDir ${indexDir} \
	--sjdbGTFfile ${gtfDir} \
	--outSJfilterReads Unique \
	--sjdbOverhang 100 \
	--twopassMode Basic \
	--runThreadN ${numcores} \
	--genomeLoad NoSharedMemory \
	--outFilterType BySJout \
	--outFilterMultimapNmax 100 \
	--outFilterMismatchNmax 33 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMismatchNoverLmax 0.3 \
	--outFilterScoreMinOverLread 0.3 \
	--limitOutSJcollapsed 1000000 \
	--limitSjdbInsertNsj 1000000 \
	--alignEndsType EndToEnd \
	--alignSJDBoverhangMin 3  \
	--alignSJoverhangMin 8 \
	--alignIntronMin 20 \
	--winAnchorMultimapNmax 50 \
	--seedSearchStartLmax 12 \
	--chimSegmentMin 20 \
	--outSAMattributes All \
	--outSAMstrandField intronMotif \
	--quantMode TranscriptomeSAM \
	--outSAMattrIHstart 0 \
	--outSAMunmapped Within \
	--outSAMtype BAM SortedByCoordinate

	
		# to protect from overwriting previous runs results with data copied for current run
		if [ ! -f ${starDir}/${1}.sorted.bam ];
		then
			echo "copying $star_tmp_Dir to $starDir"
			scp -r $star_tmp_Dir $starDir
		else 
			echo "copied in previous run"
		fi

else
	echo "Found" $starDir/${1}Aligned.sortedByCoord.out.bam
fi

find $star_tmp_Dir/${1}Aligned.sortedByCoord.out.bam -type f -size $minFileSize -delete
find $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam -type f -size $minFileSize -delete

if [[ ! -f $star_tmp_Dir/${1}Aligned.sortedByCoord.out.bam && $starDir/${1}Aligned.sortedByCoord.out.bam ]];
	then
		scp $starDir"/"${1}"Aligned.sortedByCoord.out.bam" $star_tmp_Dir
		echo "copying star bam"
		in_Bam=$star_tmp_Dir"/"${1}"Aligned.sortedByCoord.out.bam"
		echo $in_Bam
	else 
		echo "Already in tmpdir"
		in_Bam=$star_tmp_Dir"/"${1}"Aligned.sortedByCoord.out.bam"
		echo $in_Bam
fi


##########################################################################################################################
##########################################################################################################################
### Create introns BED file
##########################################################################################################################
##########################################################################################################################

# all happens off local scratch but bam in local
#resDir_tmp=$TMPDIR"/${1}_splicing"
#mkdir $resDir_tmp
# mv everything out of $resDir_tmp to $resDir at the end

if 	[ ! -f $resDir_tmp"/${1}.spliced.bed" ]; then
	bamToBed -bed12 -i $in_Bam |  awk ' $10 >1 {print $0 }' > $resDir_tmp/${1}.spliced.bed
else
	echo "Found $resDir_tmp/${1}.spliced.bed"
fi

if [ ! -f $resDir_tmp"/${1}.8names.txt" ]; then
	awk 'OFS="\t" {print $4, $11}' $resDir_tmp/${1}.spliced.bed  | sed 's/,/\t/g' | awk '$3 >=8' | awk '$2 >= 8 {print $1}' > $resDir_tmp/${1}.8names.txt
else
	echo "Found $resDir_tmp/${1}.8names.txt"
fi

if [ ! -f $resDir_tmp"/${1}.filtered_8.bed" ]; then
	LC_ALL=C grep -wF -f $resDir_tmp/${1}.8names.txt $resDir_tmp/${1}.spliced.bed > $resDir_tmp/${1}.filtered_8.bed
else
	echo "Found $resDir_tmp/${1}.filtered_8.bed"
fi

if [ ! -f $resDir_tmp"/${1}.introns.bed" ]; then	
	perl /project/RDS-SMS-brainomics-RW/rosmap/scripts/bed2introns.pl $resDir_tmp/${1}.filtered_8.bed $resDir_tmp/${1}.introns.bed
else
	echo "Found $resDir_tmp/${1}.introns.bed"
fi

##RSeQC Analysis

#juncDir_tmp=$resDir_tmp"/rseqc"
#echo $juncDir_tmp && mkdir -p $juncDir_tmp

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}filtered_8.bam" ]; then
	bedToBam -i $resDir_tmp/${1}.filtered_8.bed -g ${genSize} -bed12 > $juncDir_tmp/${1}filtered_8.bam
else
	echo "Found $resDir_tmp/${1}filtered_8.bam"
fi

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}bam_stat.txt" ]; then
        echo "Running bam_stat.py -i $juncDir_tmp/${1}filtered_8.bam > $juncDir_tmp/${1}bam_stat.txt 2>&1"
        bam_stat.py -i $juncDir_tmp/${1}filtered_8.bam > $juncDir_tmp/${1}bam_stat.txt 2>&1
else
        echo "Found" $juncDir/${1}bam_stat.txt
fi

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}read_dis.txt" ]; then
        echo "Running read_distribution.py"
        read_distribution.py -r $rseqcIndexDir -i $juncDir_tmp/${1}filtered_8.bam > $juncDir_tmp/${1}read_dis.txt 2>&1
else
        echo "Found" $juncDir/${1}read_dis.txt
fi

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}.splice_junction.pdf" ]; then
        echo "Running junction_annotation.py -r $rseqcIndexDir -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1} 2>&1"
        junction_annotation.py -r $rseqcIndexDir -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1} 2>&1
else
        echo "Found" $juncDir/${1}.splice_junction.pdf
fi

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}.junctionSaturation_plot.pdf" ]; then
        echo "Running junction_saturation.py -r $rseqcIndexDir -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1}"
        junction_saturation.py -r $rseqcIndexDir -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1}
else
        echo "Found" $juncDir/${1}.junctionSaturation_plot.pdf
fi

if [ ! -f $resDir"/${1}_splicing/rseqc/${1}filtered_8.bam.bai" ]; 
then
	echo "Running samtools index $juncDir_tmp/${1}filtered_8.bam for bam2wig.py"
	time samtools index $juncDir_tmp/${1}filtered_8.bam
else
	echo "Found" $juncDir_tmp/${1}filtered_8.bam.bai
fi
# edited because of existing output file structure - can be rewritten to something more sensible 
if [ ! -f $resDir"/${1}_splicing/rseqc"/${1}.bw ]; then
        echo "Running bam2wig.py --skip-multi-hits -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1} --chromSize /share/Temp/borgue/harvard/rseqc/hg19.chrom.sizes"
        bam2wig.py --skip-multi-hits -i $juncDir_tmp/${1}filtered_8.bam --out-prefix $juncDir_tmp/${1} --chromSize $genSize #-d "1+-,1-+,2++,2--" --wigsum=TOTAL_WIGSUM
        rm $juncDir_tmp/${1}.wig
else
        echo "Found" $juncDir/${1}.bw
fi

# necessary to allow for files from new runs to be added
# to resDir as previously mv couldn't overwrite the existing directory
if [[ ! -d $resDir"/${1}_splicing/rseqc" && ! -d $resDir"/${1}_splicing" ]]; then
		mkdir $resDir"/${1}_splicing/rseqc"
		mkdir $resDir"/${1}_splicing"
fi
mv $resDir_tmp/* $resDir"/${1}_splicing"
mv $juncDir_tmp/* $resDir"/${1}_splicing/rseqc"
##########################################################################################################################
##########################################################################################################################
### Stringtie
##########################################################################################################################
##########################################################################################################################

if [ ! -f $gownDir/${1}.gtf ]; then
        echo "Running stringtie"
        if [ ! -f $TMPDIR/${1}Aligned.sortedByCoord.out.bam ]; 	# re-work around $TMPDIR for new SGE tmp management
        then
           	scp $starDir/${1}_star/${1}Aligned.sortedByCoord.out.bam $TMPDIR 						# where in_bam is aligned to genome $starDir"/"${1}"Aligned.sortedByCoord.out.bam"
           	echo "copying bam for stringtie run 1"
        else
        	echo "Found" 
        fi

        inFile=$TMPDIR"/${1}Aligned.sortedByCoord.out.bam"
        echo "New infile" $inFile

        echo "stringtie ${inFile} -G ${gtfDir} -o ${gownDir}/${1}.gtf -a 8 -p ${numcores} -v -b ${gownDir}"
        stringtie ${inFile} -G ${gtfDir} -o ${gownDir}/${1}.gtf -a 8 -p ${numcores} -v -b ${gownDir} 				# semi-denovo over the gtf annotation
else
        echo "Found"  $gownDir/${1}.gtf
fi



##########################################################################################################################
##########################################################################################################################
###Indexing using SamTools 
##########################################################################################################################
##########################################################################################################################

if [ ! -f $in_Bam".bai" ]; 
then
	echo "Running samtools index $in_Bam"
	time samtools index $in_Bam
else
	echo "Found" $in_Bam".bai"
fi

##########################################################################################################################
##########################################################################################################################
###Sorting the transcriptome bam for RSEM
##########################################################################################################################
##########################################################################################################################

find ${starDir}/${1}.sorted.bam -type f -size $minFileSize -delete

if [ ! -f ${starDir}/${1}.sorted.bam ]; then
	echo "Sorting the transcriptome bam staroutput for RSEM"
	
	if [ ! -f $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam ];
		then
			scp $starDir"/"${1}"Aligned.sortedByCoord.out.bam" $star_tmp_Dir
			echo "copying star bam"
		else 
			echo "Already in tmpdir"
	fi

	echo "samtools view -@ ${numcores:=6} $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam -f 3 -b > ${star_tmp_Dir}/${1}.out.bam"
	time samtools view -@ ${numcores} -f 3 -b $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam > $star_tmp_Dir/${1}.out.bam

	time convert-sam-for-rsem -p ${2} --memory-per-thread 8G $star_tmp_Dir/${1}.out.bam ${star_tmp_Dir}/${1}.sorted
	echo "scp ${star_tmp_Dir}/${1}.sorted.bam ${starDir}"
	scp ${star_tmp_Dir}/${1}.sorted.bam ${starDir}
	#rm ${star_tmp_Dir}/${1}.out.bam
else
	echo "Found" ${starDir}/${1}.sorted.bam
fi

find ${starDir}/${1}.sorted.bam -type f -size $minFileSize -delete

##########################################################################################################################
##########################################################################################################################
###Counting the transcriptome bam in RSEM
##########################################################################################################################
##########################################################################################################################

find $rsemDir/${1}.genes.results -type f -size $minFileSize -delete

if [ ! -f ${rsemDir}/${1}.genes.results ]; then
	if [ ! -f ${star_tmp_Dir}/${1}.sorted.bam ];
		then
			echo "copying ${starDir}/${1}.sorted.bam ${star_tmp_Dir}"
			scp ${starDir}/${1}.sorted.bam ${star_tmp_Dir}
		else 
			echo "Already in tmpdir"
	fi

	echo "time rsem-calculate-expression --paired-end --bam --forward-prob 0 --no-bam-output -p $numcores ${star_tmp_Dir}/${1}.sorted.bam ${rsemIndexDir} ${rsemDir}/${1}"
	time rsem-calculate-expression \
		--paired-end \
		--bam \
		--forward-prob 0 \
		--no-bam-output \
		-p ${numcores} \
		${star_tmp_Dir}/${1}.sorted.bam \
		${rsemIndexDir} \
		${rsemDir}/${1}
		(echo ${1}; awk '{print $5}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_count
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $5}' $rsemDir/${1}.genes.results) > $rsemDir/${1}_rsem_gene_count
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.genes.results) > $rsemDir/${1}_rsem_gene_fpkm
else
	echo "${rsemDir}/${1}.genes.results is already there."
fi

find $rsemDir/${1}.genes.results -type f -size $minFileSize -delete




##########################################################################################################################
##########################################################################################################################
## Running RSeQC
##########################################################################################################################
##########################################################################################################################


echo "Running RSeQC"

mkdir -p ${starDir}/rseqc

if [ ! -f $star_tmp_Dir"/"${1}"Aligned.sortedByCoord.out.bam" ];
	then
		echo "copying inbam to TMP"
		scp $in_Bam $star_tmp_Dir
		scp $in_Bam".bai" $star_tmp_Dir
		in_Bam=$star_tmp_Dir"/"${1}"Aligned.sortedByCoord.out.bam"
	else 
		echo "Already in tmpdir"
fi

if [ ! -f ${starDir}/rseqc/${1}bam_stat.txt.gz ]; then
	echo "Running bam_stat.py -i ${in_Bam} > ${starDir}/rseqc/${1}bam_stat.txt 2>&1"
	time bam_stat.py -i ${in_Bam} > ${starDir}/rseqc/${1}bam_stat.txt 2>&1
else
	echo "Found" ${starDir}/rseqc/${1}bam_stat.txt.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.clipping_profile.pdf.gz ]; then
	echo "Running clipping_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time clipping_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.clipping_profile.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.deletion_profile.pdf.gz ]; then
	echo "Running deletion_profile.py -i ${in_Bam} -l $read_len --out-prefix ${starDir}/rseqc/${1}"
	time deletion_profile.py -i ${in_Bam} -l $read_len --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.deletion_profile.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.geneBodyCoverage.curves.pdf.gz ]; then
	echo "Running geneBody_coverage.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time geneBody_coverage.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/${1}.geneBodyCoverage.curves.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}infer_exp.txt.gz ]; then
	echo "Running infer_experiment.py"
	time infer_experiment.py -i ${in_Bam} -r ${rseqcIndexDir} > ${starDir}/rseqc/${1}infer_exp.txt
else
	echo "Found" ${starDir}/rseqc/${1}infer_exp.txt.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.inner_distance_plot.pdf.gz ]; then
	echo "Running inner_distance.py"
	time inner_distance.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.inner_distance_plot.pdf.gz
fi

if [ ! -f ${starDir}/${1}read_dis.txt.gz ]; then
	echo "Running read_distribution.py"
	time read_distribution.py -r ${rseqcIndexDir} -i ${in_Bam} > ${starDir}/${1}read_dis.txt 2>&1
else
	echo "Found" ${starDir}/${1}read_dis.txt.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.insertion_profile.txt.gz ]; then
	echo "Running insertion_profile.py"
	time insertion_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} 
else
	echo "Found" ${starDir}/rseqc/${1}.insertion_profile.txt.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.insertion_profile.pdf.gz ]; then
	echo "Running mismatch_profile.py"
	time mismatch_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -l $read_len
else
	echo "Found" ${starDir}/rseqc/${1}.insertion_profile.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.DupRate_plot.pdf.gz ]; then
	echo "Running read_duplication.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_duplication.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.DupRate_plot.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.GC_plot.pdf.gz ]; then
	echo "Running read_GC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_GC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.GC_plot.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.NVC_plot.pdf.gz ]; then
	echo "Running read_NVC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_NVC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.NVC_plot.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.qual.heatmap.pdf.gz ]; then
	echo "Running read_quality.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_quality.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.qual.heatmap.pdf.gz
fi


if [ ! -f ${starDir}/rseqc/${1}_read_count.xls.gz ]; then
	echo "Running RPKM_count.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--""
	time RPKM_count.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--"
else
	echo "Found" ${starDir}/${1}_read_count.xls.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.eRPKM.xls.gz ]; then
	echo "Running RPKM_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--""
	time RPKM_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--"
else
	echo "Found" ${starDir}/rseqc/${1}.eRPKM.xls.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.splice_junction.pdf.gz ]; then
	echo "Running junction_annotation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} 2>&1"
	time junction_annotation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} 2>&1
else
	echo "Found" ${starDir}/rseqc/${1}.splice_junction.pdf.gz
fi

if [ ! -f ${starDir}/rseqc/${1}.junctionSaturation_plot.pdf.gz ]; then
	echo "Running junction_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time junction_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.junctionSaturation_plot.pdf.gz
fi

cd ${starDir}

if [ ! -f ${starDir}/rseqc/${1}Aligned.sortedByCoord.out.tin.xls.gz ]; then
	echo "Running tin.py -r ${rseqcIndexDir} -i ${in_Bam}"
	time tin.py -r ${rseqcIndexDir} -i ${in_Bam}
	mv ${starDir}/${1}*.tin.xls ${starDir}/rseqc/
else
	echo "Found" ${starDir}/rseqc/${1}Aligned.sortedByCoord.out.tin.xls.gz
fi


if [ ! -f ${starDir}/rseqc/${1}.Forward.bw.gz ]; then
	echo "Running bam2wig.py --skip-multi-hits -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} --chromSize ${indexDir}/genome.chrom.sizes"
	time bam2wig.py --skip-multi-hits -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} --chromSize ${indexDir}/chrNameLength.txt --wigsum=1000000000 --strand='1++,1--,2+-,2-+'
	rm ${starDir}/rseqc/${1}*.wig
else
	echo "Found" ${starDir}/rseqc/${1}.Forward.bw.gz
fi

gzip ${starDir}/rseqc/*



##########################################################################################################################
##########################################################################################################################
###Cleanup
##########################################################################################################################
##########################################################################################################################
if [ ! -z $TMPDIR ]; then
	>&2 echo "Removing contents of directory $TMPDIR"
 	rm -rf $TMPDIR/* 
else
 	>&2 echo “[ERROR] Attempting to remove an empty variable”
fi

#mv ${outDir}/*.e* ${logDir}
#mv ${outDir}/*.o* ${logDir}


find ./*.po* -type f -size 1M -delete
find ./*.pe* -type f -size 1M -delete


asd() {
cat <<"EOT"



                      _           _        _____                      _      _           _
    /\               | |         (_)      / ____|                    | |    | |         | |
   /  \   _ __   __ _| |_   _ ___ _ ___  | |     ___  _ __ ___  _ __ | | ___| |_ ___  __| |
  / /\ \ | '_ \ / _` | | | | / __| / __| | |    / _ \| '_ ` _ \| '_ \| |/ _ \ __/ _ \/ _` |
 / ____ \| | | | (_| | | |_| \__ \ \__ \ | |___| (_) | | | | | | |_) | |  __/ ||  __/ (_| |
/_/    \_\_| |_|\__,_|_|\__, |___/_|___/  \_____\___/|_| |_| |_| .__/|_|\___|\__\___|\__,_|
                         __/ |                                 | |
                        |___/                                  |_|
EOT
}

asd
echo "Analysis finished. Have a good day"



