
git submodule update --init --recursive

mkdir -p data/refs
cd data/refs
wget https://github.com/f0t1h/3M-february-2018/raw/master/3M-february-2018.txt.gz
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz -O homo_sapiens.annot.gtf.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O homo_sapiens.cdna.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -O homo_sapiens.dna.fa.gz
gunzip homo_sapiens.*.gz
awk '$0~/^#/ || $1==21' homo_sapiens.annot.gtf > homo_sapiens.chr21.annot.gtf
awk '$0~/^>/ {FLAG=match($0,"chromosome:GRCh38:21:")} FLAG' homo_sapiens.cdna.fa > homo_sapiens.chr21.cdna.fa
awk '$0~/^>/ {FLAG=$1==">21"} FLAG' homo_sapiens.dna.fa > homo_sapiens.chr21.dna.fa
cd ../..

mkdir -p data/samples
cd data/samples
axel -n 12 http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_A549_directcDNA_replicate1_run3/SGNex_A549_directcDNA_replicate1_run3.fastq.gz
axel -n 12 http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_Hct116_directcDNA_replicate4_run1/SGNex_Hct116_directcDNA_replicate4_run1.fastq.gz
axel -n 12 http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_HepG2_directcDNA_replicate5_run3/SGNex_HepG2_directcDNA_replicate5_run3.fastq.gz
axel -n 12 http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_K562_directcDNA_replicate1_run2/SGNex_K562_directcDNA_replicate1_run2.fastq.gz
axel -n 12 http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_MCF7_cDNA_replicate1_run3/SGNex_MCF7_cDNA_replicate1_run3.fastq.gz
cd ../..

mamba env create -f env.yaml
mamba activate sc_bench
snakemake -j96 --use-conda --conda-create-envs-only