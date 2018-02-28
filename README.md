
# host_pathogen_mapping

## Requirements

- FastQC (https://github.com/s-andrews/FastQC)
- cutadapt (https://cutadapt.readthedocs.io/en/stable/)
- trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- sortMeRNA (http://bioinfo.lifl.fr/RNA/sortmerna/)
- bwa (https://sourceforge.net/projects/bio-bwa/)
- tophat/STAR (https://ccb.jhu.edu/software/tophat/index.shtml; https://github.com/alexdobin/STAR)
- samtools (http://samtools.sourceforge.net/)
- htseq-count (https://htseq.readthedocs.io/en/release_0.9.1/)

## General description of the script
The host_microbe_mapper was developed to analyse dual RNA-seq expression data with a host and a pathogen.
The entire script with all commands assigned are summarized in "host_pathogen_mapping.sh" which 
includes the commands for running Tophat and Bowtie. All the required input parameters can be specified 
in a seperate Graphical User Interface (GUI) which is provided within this data package.
For any questions please use the Github Tracking option.
![alt text](https://github.com/nthomasCUBE/host_pathogen_mapping/blob/master/misc/pix.png)
![alt text](https://github.com/nthomasCUBE/host_pathogen_mapping/blob/master/misc/GUI_v0.png)

## Assembling of the 16S genes by the output of SortMeRNA
Furthermore, we can extract the output of the sortMeRNA output form mapped reads to the 16S rRNA
gene and used that as input for the REAGO, in case that paired end sequencing information was used.
This is done by using the script 'data/extract_16.py' and by pointing to the output directory from
the host_microbe_mapper 'python data/extract_16.py <OUTPUT_DIR>'

## How to run the script
'python GUI_host_pathogen_mapping.py' open a GUI where the FastQ files can be integrated and where
the reference sequences can be defined.
It creates a 'pipeline.sh' file, that can be started and does the mapping with help of the
'data/host_microbe_mapping.sh' script.

## Demo data
We used human, mouse and the pathogen Neisseria for demonstrating the pipeline using chromosomes chr19
and by generating reads with help of ArtificialFastqGenerator.jar (Frampton et al. 2012).
Human (NC_000019.10), mouse (NC_000085.6) and N. meningitidis (NC_003112.2).

## To run demo data without GUI

It is possible to run the host_microbe_mapper also without the GUI, which requires the following.

```
FWD=Nm_cds_reads.1.fastq
RVS=Nm_cds_reads.2.fastq
REF=../read_simulator/references/Neisseria_meningitidis_genome
HUMAN=../read_simulator/references/Homo_sapiens_chr19.fasta
MOUSE=../read_simulator/references/Mus_musculus_chr19.fasta

./data/host_pathogen_mapping.sh -F ${FWD} -R ${RVS} -P ${REF} -C 4 -X /naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/output2 -H $HUMAN -I $MOUSE
```
where F represents the forwards reads, F represents the reverse, P represents the microbial genome,
whereas C represents the processor amount, X represents the output directory, H represents the first host genome
(here: human) reference and I represents the second host reference (here: mouse).

Parameter | Meaning
------------ | -------------
F | forward reads in the FASTQ format
R | reverse reads in the FASTQ format
P | microbe reference in FASTA format
C | processor amount 
X | output directory
H | human mapping reference
I | mouse mapping reference









