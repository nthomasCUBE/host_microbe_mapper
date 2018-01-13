#!/usr/bin/bash
#
#SBATCH --job-name=test
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --nice=3000
#SBATCH --partition=basic

ulimit -v 30000000
tmp=$(mktemp)
tmp2=$(mktemp)

echo "tmp::${tmp}"
echo "tmp2::${tmp2}"

zcat BLA_R1_001.fastq.gz | head -n 10000 > $tmp
zcat BLA_R2_001.fastq.gz | head -n 10000 > $tmp2

wc -l $tmp
wc -l $tmp2

FWD=$tmp
RVS=$tmp2

X=/proj/Neisseria_meningitidis/5.mapping_pipeline/demo_data

HS_INDEX=HOMO_SAPIENS_INDEX
MS_INDEX=MUS_MUSCULUS_INDEX

NM_REF=references/GCF_000026965.1_ASM2696v1_genomic.fna
HS_REF=references/Homo_sapiens.GRCh38.dna.primary_assembly.fa
MS_REF=references/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa

NM_GFF=gff/GCF_000026965.1_ASM2696v1_genomic_V2.gff
HS_GFF=gff/Homo_sapiens.GRCh38.87_modified.gff3
MS_GFF=gff/Mus_musculus.GRCm38.87_modified.gff3

C=1

echo "FWD::$FWD"
echo "RVS::$RVS"
echo "X::$X"
echo "P::${NM_REF}"
echo "O::star"
echo "A::${HS_INDEX}"
echo "B::${MS_INDEX}"
echo "H::${HS_REF}"
echo "I::${MS_REF}"
echo "J::${NM_GFF}"
echo "K::${HS_GFF}"
echo "L::${MS_GFF}"

bash host_microbe_mapping.sh -F ${FWD} -R ${RVS} -X ${X} -P ${NM_REF} -O star -C ${C} -A ${HS_INDEX} -B ${MS_INDEX} -H ${HS_REF} -I ${MS_REF} -J ${NM_GFF} -K ${HS_GFF} -L ${MS_GFF}

