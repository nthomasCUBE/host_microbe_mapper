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
`
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

Species | Gene Bank Identifier
------------ | -------------
Homo sapiens | GCF_000001405.37 - 38.p11
Mus musculus | GCF_000001635.26 - 38.p6
Neisseria meningitidis | GCF_000008805.1 - ASM880v1

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

## Example output of the hostmicrobemapper

```
-------------------------------------------------------------------

               HOST_PATHOGEN_MAPPING v0.0.7
                  14 DECEMBER 2017

-------------------------------------------------------------------
FWD.fastq
RVS.fastq
../read_simulator/references/Neisseria_meningitidis_genome
NUMBER OF CPUs: 4
FORWARD READS: FWD.fastq
REVERSE READS: RVS.fastq
Output dir: /naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/output2_28feb18
Pathogen genome: ../read_simulator/references/Neisseria_meningitidis_genome
Eukaryotic genome: ../read_simulator/references/Homo_sapiens_chr19.fasta
Eukaryotic genome: ../read_simulator/references/Mus_musculus_chr19.fasta
The mapping tools 'tophat' was selected for mapping RNA-seq in the host genomes



WORKING directory: /gpfs/proj/abc/tmp.qtfwssmZHr
no subsets - complete set to use
step-0: fastqc
/gpfs/proj/abc/tmp.qtfwssmZHr/FWD.fastq
/gpfs/proj/abc/tmp.qtfwssmZHr
Started analysis of FWD.fastq
Approx 10% complete for FWD.fastq
Approx 25% complete for FWD.fastq
Approx 40% complete for FWD.fastq
Approx 50% complete for FWD.fastq
Approx 65% complete for FWD.fastq
Approx 80% complete for FWD.fastq
Approx 90% complete for FWD.fastq
Started analysis of RVS.fastq
Approx 10% complete for RVS.fastq
Approx 25% complete for RVS.fastq
Approx 40% complete for RVS.fastq
Approx 50% complete for RVS.fastq
Approx 65% complete for RVS.fastq
Approx 80% complete for RVS.fastq
Approx 90% complete for RVS.fastq
step-1: cutadapt
=paired=,/gpfs/proj/abc/tmp.qtfwssmZHr/FWD.fastq,/gpfs/proj/abc/tmp.qtfwssmZHr/RVS.fastq
step-2: trimmomatic
TrimmomaticPE: Started with arguments:
 -threads 4 -phred33 /gpfs/proj/abc/tmp.qtfwssmZHr/FWD.fastq /gpfs/proj/abc/tmp.qtfwssmZHr/RVS.fastq /gpfs/proj/abc/tmp.qtfwssmZHr/trimmomatic_forward_paired.fq.gz /gpfs/proj/abc/tmp.qtfwssmZHr/trimmomatic_forward_unpaired.fq.gz /gpfs/proj/abc/tmp.qtfwssmZHr/trimmomatic_reverse_paired.fq.gz /gpfs/proj/abc/tmp.qtfwssmZHr/trimmomatic_reverse_unpaired.fq.gz LEADING:8 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
Input Read Pairs: 7500 Both Surviving: 7500 (100.00%) Forward Only Surviving: 0 (0.00%) Reverse Only Surviving: 0 (0.00%) Dropped: 0 (0.00%)
TrimmomaticPE: Completed successfully
step-3: sortMeRNA
/naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/sortmerna-2.1b
step-4: bwa (mapping to bacterium)
[bwa_index] Pack FASTA... 0.02 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.50 seconds elapse.
[bwa_index] Update BWT... 0.02 sec
[bwa_index] Pack forward-only FASTA... 0.01 sec
[bwa_index] Construct SA from BWT and Occ... 0.28 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME
[main] Real time: 1.301 sec; CPU: 0.840 sec
14970 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
3858 + 0 mapped (25.77% : N/A)
14970 + 0 paired in sequencing
7485 + 0 read1
7485 + 0 read2
2674 + 0 properly paired (17.86% : N/A)
2848 + 0 with itself and mate mapped
1010 + 0 singletons (6.75% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
[bam_sort_core] merging from 0 files and 4 in-memory blocks...

GENOME=/gpfs/proj/abc/tmp.qtfwssmZHr/GENOME
DIR=/gpfs/proj/abc/tmp.qtfwssmZHr
step-5: tophat (mapping to Host)
Building a SMALL index

[2018-02-28 17:29:36] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-02-28 17:29:36] Checking for Bowtie
                  Bowtie version:        2.3.3.1
[2018-02-28 17:29:36] Checking for Bowtie index files (genome)..
[2018-02-28 17:29:36] Checking for reference FASTA file
        Warning: Could not find FASTA file /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME2.fa
[2018-02-28 17:29:36] Reconstituting reference FASTA file from Bowtie index
  Executing: /naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/bowtie2-2.3.3.1-linux-x86_64/bowtie2-inspect /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME2 > /gpfs/proj/abc/tmp.qtfwssmZHr/tophat/tmp/GENOME2.fa
[2018-02-28 17:29:39] Generating SAM header for /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME2
[2018-02-28 17:29:39] Preparing reads
         left reads: min. length=186, max. length=200, 5052 kept reads (0 discarded)
        right reads: min. length=194, max. length=200, 5052 kept reads (0 discarded)
[2018-02-28 17:29:39] Mapping left_kept_reads to genome GENOME2 with Bowtie2
[2018-02-28 17:29:41] Mapping left_kept_reads_seg1 to genome GENOME2 with Bowtie2 (1/8)
[2018-02-28 17:29:41] Mapping left_kept_reads_seg2 to genome GENOME2 with Bowtie2 (2/8)
[2018-02-28 17:29:41] Mapping left_kept_reads_seg3 to genome GENOME2 with Bowtie2 (3/8)
[2018-02-28 17:29:42] Mapping left_kept_reads_seg4 to genome GENOME2 with Bowtie2 (4/8)
[2018-02-28 17:29:42] Mapping left_kept_reads_seg5 to genome GENOME2 with Bowtie2 (5/8)
[2018-02-28 17:29:42] Mapping left_kept_reads_seg6 to genome GENOME2 with Bowtie2 (6/8)
[2018-02-28 17:29:43] Mapping left_kept_reads_seg7 to genome GENOME2 with Bowtie2 (7/8)
[2018-02-28 17:29:43] Mapping left_kept_reads_seg8 to genome GENOME2 with Bowtie2 (8/8)
[2018-02-28 17:29:44] Mapping right_kept_reads to genome GENOME2 with Bowtie2
[2018-02-28 17:29:45] Mapping right_kept_reads_seg1 to genome GENOME2 with Bowtie2 (1/8)
[2018-02-28 17:29:46] Mapping right_kept_reads_seg2 to genome GENOME2 with Bowtie2 (2/8)
[2018-02-28 17:29:46] Mapping right_kept_reads_seg3 to genome GENOME2 with Bowtie2 (3/8)
[2018-02-28 17:29:47] Mapping right_kept_reads_seg4 to genome GENOME2 with Bowtie2 (4/8)
[2018-02-28 17:29:47] Mapping right_kept_reads_seg5 to genome GENOME2 with Bowtie2 (5/8)
[2018-02-28 17:29:47] Mapping right_kept_reads_seg6 to genome GENOME2 with Bowtie2 (6/8)
[2018-02-28 17:29:48] Mapping right_kept_reads_seg7 to genome GENOME2 with Bowtie2 (7/8)
[2018-02-28 17:29:48] Mapping right_kept_reads_seg8 to genome GENOME2 with Bowtie2 (8/8)
[2018-02-28 17:29:49] Searching for junctions via segment mapping
[2018-02-28 17:29:51] Retrieving sequences for splices
[2018-02-28 17:29:52] Indexing splices
Building a SMALL index
[2018-02-28 17:29:53] Mapping left_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 17:29:53] Mapping left_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 17:29:54] Mapping left_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 17:29:54] Mapping left_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 17:29:54] Mapping left_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 17:29:54] Mapping left_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 17:29:55] Mapping left_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 17:29:55] Mapping left_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 17:29:55] Joining segment hits
[2018-02-28 17:29:57] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 17:29:57] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 17:29:58] Mapping right_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 17:29:58] Mapping right_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 17:29:58] Mapping right_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 17:29:58] Mapping right_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 17:29:59] Mapping right_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 17:29:59] Mapping right_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 17:29:59] Joining segment hits
[2018-02-28 17:30:01] Reporting output tracks
-----------------------------------------------
[2018-02-28 17:30:03] A summary of the alignment counts can be found in /gpfs/proj/abc/tmp.qtfwssmZHr/tophat/align_summary.txt
[2018-02-28 17:30:03] Run complete: 00:00:27 elapsed
Left reads:
          Input     :      5052
           Mapped   :      2120 (42.0% of input)
            of these:        80 ( 3.8%) have multiple alignments (0 have >20)
Right reads:
          Input     :      5052
           Mapped   :      2083 (41.2% of input)
            of these:        81 ( 3.9%) have multiple alignments (0 have >20)
41.6% overall read mapping rate.

Aligned pairs:      1745
     of these:        67 ( 3.8%) have multiple alignments
34.5% concordant pair alignment rate.
  adding: gpfs/proj/abc/tmp.qtfwssmZHr/tophat/accepted_hits.bam (deflated 9%)
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
step-6: tophat (mapping to Host-2)
Building a SMALL index

[2018-02-28 17:30:36] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-02-28 17:30:36] Checking for Bowtie
                  Bowtie version:        2.3.3.1
[2018-02-28 17:30:36] Checking for Bowtie index files (genome)..
[2018-02-28 17:30:36] Checking for reference FASTA file
        Warning: Could not find FASTA file /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME3.fa
[2018-02-28 17:30:36] Reconstituting reference FASTA file from Bowtie index
  Executing: /naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/bowtie2-2.3.3.1-linux-x86_64/bowtie2-inspect /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME3 > /gpfs/proj/abc/tmp.qtfwssmZHr/tophat2/tmp/GENOME3.fa
[2018-02-28 17:30:39] Generating SAM header for /gpfs/proj/abc/tmp.qtfwssmZHr/GENOME3
[2018-02-28 17:30:40] Preparing reads
         left reads: min. length=186, max. length=200, 5052 kept reads (0 discarded)
        right reads: min. length=194, max. length=200, 5052 kept reads (0 discarded)
[2018-02-28 17:30:40] Mapping left_kept_reads to genome GENOME3 with Bowtie2
[2018-02-28 17:30:41] Mapping left_kept_reads_seg1 to genome GENOME3 with Bowtie2 (1/8)
[2018-02-28 17:30:41] Mapping left_kept_reads_seg2 to genome GENOME3 with Bowtie2 (2/8)
[2018-02-28 17:30:42] Mapping left_kept_reads_seg3 to genome GENOME3 with Bowtie2 (3/8)
[2018-02-28 17:30:42] Mapping left_kept_reads_seg4 to genome GENOME3 with Bowtie2 (4/8)
[2018-02-28 17:30:43] Mapping left_kept_reads_seg5 to genome GENOME3 with Bowtie2 (5/8)
[2018-02-28 17:30:43] Mapping left_kept_reads_seg6 to genome GENOME3 with Bowtie2 (6/8)
[2018-02-28 17:30:43] Mapping left_kept_reads_seg7 to genome GENOME3 with Bowtie2 (7/8)
[2018-02-28 17:30:44] Mapping left_kept_reads_seg8 to genome GENOME3 with Bowtie2 (8/8)
[2018-02-28 17:30:44] Mapping right_kept_reads to genome GENOME3 with Bowtie2
[2018-02-28 17:30:46] Mapping right_kept_reads_seg1 to genome GENOME3 with Bowtie2 (1/8)
[2018-02-28 17:30:46] Mapping right_kept_reads_seg2 to genome GENOME3 with Bowtie2 (2/8)
[2018-02-28 17:30:46] Mapping right_kept_reads_seg3 to genome GENOME3 with Bowtie2 (3/8)
[2018-02-28 17:30:47] Mapping right_kept_reads_seg4 to genome GENOME3 with Bowtie2 (4/8)
[2018-02-28 17:30:47] Mapping right_kept_reads_seg5 to genome GENOME3 with Bowtie2 (5/8)
[2018-02-28 17:30:47] Mapping right_kept_reads_seg6 to genome GENOME3 with Bowtie2 (6/8)
[2018-02-28 17:30:48] Mapping right_kept_reads_seg7 to genome GENOME3 with Bowtie2 (7/8)
[2018-02-28 17:30:48] Mapping right_kept_reads_seg8 to genome GENOME3 with Bowtie2 (8/8)
[2018-02-28 17:30:49] Searching for junctions via segment mapping
[2018-02-28 17:30:50] Retrieving sequences for splices
[2018-02-28 17:30:52] Indexing splices
Building a SMALL index
[2018-02-28 17:30:53] Mapping left_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 17:30:53] Mapping left_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 17:30:53] Mapping left_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 17:30:54] Mapping left_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 17:30:54] Mapping left_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 17:30:54] Mapping left_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 17:30:55] Mapping left_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 17:30:55] Mapping left_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 17:30:55] Joining segment hits
[2018-02-28 17:30:57] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 17:30:57] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 17:30:57] Mapping right_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 17:30:58] Mapping right_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 17:30:58] Mapping right_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 17:30:58] Mapping right_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 17:30:59] Mapping right_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 17:30:59] Mapping right_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 17:30:59] Joining segment hits
[2018-02-28 17:31:01] Reporting output tracks
-----------------------------------------------
[2018-02-28 17:31:03] A summary of the alignment counts can be found in /gpfs/proj/abc/tmp.qtfwssmZHr/tophat2/align_summary.txt
[2018-02-28 17:31:04] Run complete: 00:00:27 elapsed
  adding: gpfs/proj/abc/tmp.qtfwssmZHr/tophat2/accepted_hits.bam (deflated 8%)
Left reads:
          Input     :      5052
           Mapped   :      2196 (43.5% of input)
            of these:        31 ( 1.4%) have multiple alignments (0 have >20)
Right reads:
          Input     :      5052
           Mapped   :      2266 (44.9% of input)
            of these:        34 ( 1.5%) have multiple alignments (0 have >20)
44.2% overall read mapping rate.

Aligned pairs:      1978
     of these:        29 ( 1.5%) have multiple alignments
39.2% concordant pair alignment rate.
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
No count extraction for - pathogen
No count extraction for - host1
No count extraction for - host2
INFO: Calculations are stored under: tmp.qtfwssmZHr
OUTPUT=,output2_28feb18
```







