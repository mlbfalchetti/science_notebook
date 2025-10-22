#-- 
ssh maf0964@o2.hms.harvard.edu
srun --pty -p interactive -t 0-12:00 --mem=200G -c 20 bash

#-- Prepare
conda create -n rnaseq star rsem samtools fastp fastqc rseqc pigz -c bioconda -c conda-forge -y

#-- Run
conda activate rnaseq
export THREADS=64
mkdir -p ref fastqc trim work
cd ref
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz
gunzip -c GRCm39.primary_assembly.genome.fa.gz > GRCm39.fa
gunzip -c gencode.vM35.annotation.gtf.gz > gencode.vM35.annotation.gtf
rsem-prepare-reference \
  --gtf gencode.vM35.annotation.gtf \
  --star \
  -p $THREADS \
  GRCm39.fa rsem_grcm39_vM35
cd ..
export REFPREFIX="$(pwd)/ref/rsem_grcm39_vM35"
SAMPLES=($(ls *_R1_001.fastq.gz | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$//' | sort -u))
printf '%s\n' "${SAMPLES[@]}"
# SAMPLES=("${SAMPLES[@]:0:4}")
for S in "${SAMPLES[@]}"; do
  fastqc --nogroup -t "$THREADS" -o fastqc ${S}_S*_L*_R1_001.fastq.gz ${S}_S*_L*_R2_001.fastq.gz
done
for S in "${SAMPLES[@]}"; do
  fastp \
    -i ${S}_S*_L*_R1_001.fastq.gz -I ${S}_S*_L*_R2_001.fastq.gz \
    -o trim/${S}_R1.trim.fq.gz   -O trim/${S}_R2.trim.fq.gz \
    -j fastqc/${S}_fastp.json    -h fastqc/${S}_fastp.html \
    -w "$THREADS"
done
for S in "${SAMPLES[@]}"; do
  fastqc --nogroup -t "$THREADS" -o fastqc trim/${S}_R1.trim.fq.gz trim/${S}_R2.trim.fq.gz
done
for S in "${SAMPLES[@]}"; do
  rsem-calculate-expression \
    --star --star-gzipped-read-file \
    --paired-end \
    --strandedness reverse \
    --append-names --estimate-rspd \
    -p "$THREADS" \
    trim/${S}_R1.trim.fq.gz trim/${S}_R2.trim.fq.gz \
    "$REFPREFIX" work/${S}
done
