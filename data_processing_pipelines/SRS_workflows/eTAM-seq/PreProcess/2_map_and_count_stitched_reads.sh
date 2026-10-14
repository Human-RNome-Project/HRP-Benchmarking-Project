#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## main setting
# All the reference files live under $ref_dir. Override it from the environment
# (REF_DIR=/path/to/b38 ./run) or as the second argument of this script.
# See README.md ("Reference data: where to get it") for the expected content.
ref_dir=${2:-${REF_DIR:-/home/h/OneDrive/Unil/data/b38}}
# sample fastq file path
fastq_dir=${FASTQ_DIR:-.}
# hisat-3n executable
hisat3n_bin=${HISAT3N:-$HOME/src/github/hisat-3n/hisat-3n}
# hisat-3n rRNA index name: e.g., rRNA
hisat3n_rep_index_name=${REP_INDEX_NAME:-hs_rrna}
# hisat-3n rRNA index path
hisat3n_rep_index=${REP_INDEX:-$ref_dir/HS_rep/$hisat3n_rep_index_name}
# hisat-3n genome index name: e.g., GRCh38_tran
hisat3n_index_name=${GENOME_INDEX_NAME:-b38}
# hisat-3n genome index path
hisat3n_index=${GENOME_INDEX:-$ref_dir/HS_rep/$hisat3n_index_name}
# hisat-3n alignment result output path
hisat3n_out_dir=${OUT_DIR:-./outstitch}
# rRNA sequneces in fasta format
rep_fa=${REP_FA:-$ref_dir/hs_rrna.fa}
# genome sequneces in fasta format
genome_fa=${GENOME_FA:-$ref_dir/b38.fa}
# threads to run the bash script
ncpus=${NCPUS:-24}
# minimum percentage of converted As in a read to be considered as a valid read
A2G_percent=${A2G_PERCENT:-0.5}
# minimum read coverage in a position to be considered for output
cov=${COV:-1}
# sample name: e.g., hela.polya.wt.ftom.rep1
sample=$1
# library type: F for forward, R for reverse
strandness=${STRANDNESS:-R}
# here R2 reads were used.
sample_r2_fq_file=$fastq_dir/$sample.stitched.fq.gz

if [[ -z $sample ]]; then
	echo "usage: $0 <sample> [ref_dir]   (ref_dir defaults to \$REF_DIR)" >&2
	exit 1
fi

# fail early, with an explicit message, if the reference folder is incomplete
for f in "$rep_fa" "$genome_fa"; do
	if [[ ! -s $f ]]; then
		echo "Missing reference file: $f" >&2
		echo "Check REF_DIR (currently '$ref_dir') -- see README.md, 'Reference data: where to get it'." >&2
		exit 1
	fi
done
for idx in "$hisat3n_rep_index" "$hisat3n_index"; do
	# .ht2 for a normal index, .ht2l for a large one (human genome)
	if ! compgen -G "$idx*.ht2*" > /dev/null; then
		echo "Missing hisat-3n index: $idx*.ht2[l]" >&2
		echo "Check REF_DIR (currently '$ref_dir') -- see README.md, 'Reference data: where to get it'." >&2
		exit 1
	fi
done

echo -e "\n$sample\n"
echo "reference folder: $ref_dir"

if [[ -d $hisat3n_out_dir/$sample ]];then rm -rf $hisat3n_out_dir/$sample;fi
mkdir -p $hisat3n_out_dir/$sample

echo -e "\nhisat3n align -- rRNA mapping\n"
$hisat3n_bin -p $ncpus --time --base-change A,G --no-spliced-alignment --no-softclip --norc --no-unal --rna-strandness $strandness -x $hisat3n_rep_index -U $sample_r2_fq_file --un-gz $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
	samtools sort -@ $ncpus -T $hisat3n_out_dir/$sample -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam -

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam

samtools stats -r $rep_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam > $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.stats

# filter reads by the percentage of converted As
samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -
# Count converted As and unconverted As. Here hisat-3n-table is used as an example. Other similar tools can also be used to complete the task. 
#samtools view -@ $ncpus $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam | hisat-3n-table -u -p $ncpus --alignments - --ref $rep_fa --output-name $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.conversion.flt.txt --base-change A,G
# Here, pileup2var is used.
pileup2var -f 524 -s $strandness -g $rep_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.pileup2var.flt.txt

echo -e "\nhisat3n align -- genome mapping\n"
$hisat3n_bin -p $ncpus --time --base-change A,G --repeat --repeat-limit 1000 --bowtie2-dp 0 --no-unal --rna-strandness $strandness -x $hisat3n_index -U $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
	samtools view -@ $ncpus -Shb - -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam

# filter out multiple-loci mapped reads
samtools view -@ $ncpus -q 60 -hb $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam

# deduplication
umi_tools dedup --random-seed=123 --method=unique --spliced-is-unique -I $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -S $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam \
	--output-stats=$hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup -L $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.log

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam

samtools stats -r $genome_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam > $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.stats

# filter reads by the percentage of converted As
samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -
# Count converted As and unconverted As. Here hisat-3n-table is used as an example. Other similar tools can also be used to complete the task. 
#samtools view -@ $ncpus $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam | hisat-3n-table -p $ncpus --alignments - --ref $genome_fa --output-name $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.conversion.flt.txt --base-change A,G
# Here, pileup2var is used.
pileup2var -t $ncpus -f 524 -a A -c $cov -s $strandness -g $genome_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.pileup2var.flt.txt

rm $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam.bai
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam.bai
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam.bai
wait

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
