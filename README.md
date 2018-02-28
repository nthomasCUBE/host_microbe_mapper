
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
Output dir: /naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/output2
Pathogen genome: ../read_simulator/references/Neisseria_meningitidis_genome
Eukaryotic genome: ../read_simulator/references/Homo_sapiens_chr19.fasta
Eukaryotic genome: ../read_simulator/references/Mus_musculus_chr19.fasta
The mapping tools 'tophat' was selected for mapping RNA-seq in the host genomes



WORKING directory: /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx
step-0: fastqc
/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/FWD.fastq
/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx
Started analysis of FWD.fastq
Approx 5% complete for FWD.fastq
Approx 10% complete for FWD.fastq
Approx 15% complete for FWD.fastq
Approx 20% complete for FWD.fastq
Approx 25% complete for FWD.fastq
Approx 30% complete for FWD.fastq
Approx 35% complete for FWD.fastq
Approx 40% complete for FWD.fastq
Approx 45% complete for FWD.fastq
Approx 50% complete for FWD.fastq
Approx 55% complete for FWD.fastq
Approx 60% complete for FWD.fastq
Approx 65% complete for FWD.fastq
Approx 70% complete for FWD.fastq
Approx 75% complete for FWD.fastq
Approx 80% complete for FWD.fastq
Approx 85% complete for FWD.fastq
Approx 90% complete for FWD.fastq
Approx 95% complete for FWD.fastq
Started analysis of RVS.fastq
Approx 5% complete for RVS.fastq
Approx 10% complete for RVS.fastq
Approx 15% complete for RVS.fastq
Approx 20% complete for RVS.fastq
Approx 25% complete for RVS.fastq
Approx 30% complete for RVS.fastq
Approx 35% complete for RVS.fastq
Approx 40% complete for RVS.fastq
Approx 45% complete for RVS.fastq
Approx 50% complete for RVS.fastq
Approx 55% complete for RVS.fastq
Approx 60% complete for RVS.fastq
Approx 65% complete for RVS.fastq
Approx 70% complete for RVS.fastq
Approx 75% complete for RVS.fastq
Approx 80% complete for RVS.fastq
Approx 85% complete for RVS.fastq
Approx 90% complete for RVS.fastq
Approx 95% complete for RVS.fastq
step-1: cutadapt
=paired=,/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/FWD.fastq,/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/RVS.fastq
cutadapt: error: Reads are improperly paired. Read name 'HWI-ST745_0097:7:1101:1001:1000#0/1' in file 1 does not match 'HWI-ST745_0097:7:1101:9439:1017#0/2' in file 2.
step-2: trimmomatic
TrimmomaticPE: Started with arguments:
 -threads 4 -phred33 /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/FWD.fastq /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/RVS.fastq /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/trimmomatic_forward_paired.fq.gz /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/trimmomatic_forward_unpaired.fq.gz /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/trimmomatic_reverse_paired.fq.gz /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/trimmomatic_reverse_unpaired.fq.gz LEADING:8 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
Input Read Pairs: 467260 Both Surviving: 467259 (100.00%) Forward Only Surviving: 0 (0.00%) Reverse Only Surviving: 1 (0.00%) Dropped: 0 (0.00%)
TrimmomaticPE: Completed successfully
step-3: sortMeRNA
/naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/test/sortmerna-2.1b
step-4: bwa (mapping to bacterium)
[bwa_index] Pack FASTA... 0.01 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 0.52 seconds elapse.
[bwa_index] Update BWT... 0.01 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.30 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME
[main] Real time: 2.179 sec; CPU: 0.872 sec
1846764 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
433226 + 0 mapped (23.46% : N/A)
1846764 + 0 paired in sequencing
923382 + 0 read1
923382 + 0 read2
329968 + 0 properly paired (17.87% : N/A)
342012 + 0 with itself and mate mapped
91214 + 0 singletons (4.94% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
GENOME=/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME
DIR=/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx
step-5: tophat (mapping to Host)
Building a SMALL index

[2018-02-28 15:16:58] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-02-28 15:16:58] Checking for Bowtie
                  Bowtie version:        2.3.3.1
[2018-02-28 15:16:58] Checking for Bowtie index files (genome)..
[2018-02-28 15:16:58] Checking for reference FASTA file
        Warning: Could not find FASTA file /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME2.fa
[2018-02-28 15:16:58] Reconstituting reference FASTA file from Bowtie index
  Executing: /naslx/projects_mpiio/pr74ma/ge34juq2/apps/bowtie2-2.3.3.1-linux-x86_64/bowtie2-inspect /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME2 > /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/tmp/GENOME2.fa
[2018-02-28 15:17:00] Generating SAM header for /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME2
[2018-02-28 15:17:01] Preparing reads
         left reads: min. length=186, max. length=200, 706769 kept reads (0 discarded)
        right reads: min. length=186, max. length=200, 706769 kept reads (0 discarded)
[2018-02-28 15:17:20] Mapping left_kept_reads to genome GENOME2 with Bowtie2
[2018-02-28 15:18:47] Mapping left_kept_reads_seg1 to genome GENOME2 with Bowtie2 (1/8)
[2018-02-28 15:18:54] Mapping left_kept_reads_seg2 to genome GENOME2 with Bowtie2 (2/8)
[2018-02-28 15:19:00] Mapping left_kept_reads_seg3 to genome GENOME2 with Bowtie2 (3/8)
[2018-02-28 15:19:07] Mapping left_kept_reads_seg4 to genome GENOME2 with Bowtie2 (4/8)
[2018-02-28 15:19:13] Mapping left_kept_reads_seg5 to genome GENOME2 with Bowtie2 (5/8)
[2018-02-28 15:19:19] Mapping left_kept_reads_seg6 to genome GENOME2 with Bowtie2 (6/8)
[2018-02-28 15:19:25] Mapping left_kept_reads_seg7 to genome GENOME2 with Bowtie2 (7/8)
[2018-02-28 15:19:31] Mapping left_kept_reads_seg8 to genome GENOME2 with Bowtie2 (8/8)
[2018-02-28 15:19:37] Mapping right_kept_reads to genome GENOME2 with Bowtie2
[2018-02-28 15:21:10] Mapping right_kept_reads_seg1 to genome GENOME2 with Bowtie2 (1/8)
[2018-02-28 15:21:16] Mapping right_kept_reads_seg2 to genome GENOME2 with Bowtie2 (2/8)
[2018-02-28 15:21:22] Mapping right_kept_reads_seg3 to genome GENOME2 with Bowtie2 (3/8)
[2018-02-28 15:21:28] Mapping right_kept_reads_seg4 to genome GENOME2 with Bowtie2 (4/8)
[2018-02-28 15:21:34] Mapping right_kept_reads_seg5 to genome GENOME2 with Bowtie2 (5/8)
[2018-02-28 15:21:41] Mapping right_kept_reads_seg6 to genome GENOME2 with Bowtie2 (6/8)
[2018-02-28 15:21:48] Mapping right_kept_reads_seg7 to genome GENOME2 with Bowtie2 (7/8)
[2018-02-28 15:21:54] Mapping right_kept_reads_seg8 to genome GENOME2 with Bowtie2 (8/8)
[2018-02-28 15:22:02] Searching for junctions via segment mapping
[2018-02-28 15:22:24] Retrieving sequences for splices
[2018-02-28 15:22:26] Indexing splices
Building a SMALL index
[2018-02-28 15:22:27] Mapping left_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 15:22:30] Mapping left_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 15:22:32] Mapping left_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 15:22:34] Mapping left_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 15:22:37] Mapping left_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 15:22:39] Mapping left_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 15:22:41] Mapping left_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 15:22:44] Mapping left_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 15:22:46] Joining segment hits
[2018-02-28 15:22:50] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 15:22:52] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 15:22:54] Mapping right_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 15:22:56] Mapping right_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 15:22:59] Mapping right_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 15:23:01] Mapping right_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 15:23:03] Mapping right_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 15:23:06] Mapping right_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 15:23:08] Joining segment hits
[2018-02-28 15:23:12] Reporting output tracks
        [FAILED]
Error running /naslx/projects_mpiio/pr74ma/ge34juq2/apps/tophat-2.1.1.Linux_x86_64/tophat_reports --min-anchor 8 --splice-mismatches 0 --min-report-intron 50 --max-report-intron 500000 --min-isoform-fraction 0.15 --output-dir /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/ --max-multihits 20 --max-seg-multihits 40 --segment-length 25 --segment-mismatches 2 --min-closure-exon 100 --min-closure-intron 50 --max-closure-intron 5000 --min-coverage-intron 50 --max-coverage-intron 20000 --min-segment-intron 50 --max-segment-intron 500000 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --read-realign-edit-dist 3 --max-insertion-length 3 --max-deletion-length 3 -z gzip -p4 --inner-dist-mean 50 --inner-dist-std-dev 20 --no-closure-search --no-coverage-search --no-microexon-search --sam-header /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/tmp/GENOME2_genome.bwt.samheader.sam --report-discordant-pair-alignments --report-mixed-alignments --samtools=/naslx/projects_mpiio/pr74ma/ge34juq2/apps/tophat-2.1.1.Linux_x86_64/samtools_0.1.18 --bowtie2-max-penalty 6 --bowtie2-min-penalty 2 --bowtie2-penalty-for-N 1 --bowtie2-read-gap-open 5 --bowtie2-read-gap-cont 3 --bowtie2-ref-gap-open 5 --bowtie2-ref-gap-cont 3 /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/tmp/GENOME2.fa /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/junctions.bed /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/insertions.bed /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/deletions.bed /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/fusions.out /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/tmp/accepted_hits  /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/tmp/left_kept_reads.bam
        Loading ...done

tail: cannot open `/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/*align*' for reading: No such file or directory
        zip warning: name not matched: /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/accepted*bam

zip error: Nothing to do! (/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/eukaryotic_host1_mapped.bam.zip)
[E::hts_open_format] Failed to open file /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/unmapped.bam
samtools sort: can't open "/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/unmapped.bam": No such file or directory
[E::hts_open_format] Failed to open file /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/accepted*bam
samtools flagstat: Cannot open input file "/gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat/accepted*bam": No such file or directory
step-6: tophat (mapping to Host-2)
Building a SMALL index

[2018-02-28 15:23:49] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-02-28 15:23:49] Checking for Bowtie
                  Bowtie version:        2.3.3.1
[2018-02-28 15:23:49] Checking for Bowtie index files (genome)..
[2018-02-28 15:23:49] Checking for reference FASTA file
        Warning: Could not find FASTA file /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME3.fa
[2018-02-28 15:23:49] Reconstituting reference FASTA file from Bowtie index
  Executing: /naslx/projects_mpiio/pr74ma/ge34juq2/apps/bowtie2-2.3.3.1-linux-x86_64/bowtie2-inspect /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME3 > /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat2/tmp/GENOME3.fa
[2018-02-28 15:23:52] Generating SAM header for /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/GENOME3
[2018-02-28 15:23:53] Preparing reads
         left reads: min. length=186, max. length=200, 706769 kept reads (0 discarded)
        right reads: min. length=186, max. length=200, 706769 kept reads (0 discarded)
[2018-02-28 15:24:12] Mapping left_kept_reads to genome GENOME3 with Bowtie2
[2018-02-28 15:26:11] Mapping left_kept_reads_seg1 to genome GENOME3 with Bowtie2 (1/8)
[2018-02-28 15:26:18] Mapping left_kept_reads_seg2 to genome GENOME3 with Bowtie2 (2/8)
[2018-02-28 15:26:27] Mapping left_kept_reads_seg3 to genome GENOME3 with Bowtie2 (3/8)
[2018-02-28 15:26:34] Mapping left_kept_reads_seg4 to genome GENOME3 with Bowtie2 (4/8)
[2018-02-28 15:26:41] Mapping left_kept_reads_seg5 to genome GENOME3 with Bowtie2 (5/8)
[2018-02-28 15:26:48] Mapping left_kept_reads_seg6 to genome GENOME3 with Bowtie2 (6/8)
[2018-02-28 15:26:54] Mapping left_kept_reads_seg7 to genome GENOME3 with Bowtie2 (7/8)
[2018-02-28 15:27:01] Mapping left_kept_reads_seg8 to genome GENOME3 with Bowtie2 (8/8)
[2018-02-28 15:27:08] Mapping right_kept_reads to genome GENOME3 with Bowtie2
[2018-02-28 15:28:41] Mapping right_kept_reads_seg1 to genome GENOME3 with Bowtie2 (1/8)
[2018-02-28 15:28:49] Mapping right_kept_reads_seg2 to genome GENOME3 with Bowtie2 (2/8)
[2018-02-28 15:28:56] Mapping right_kept_reads_seg3 to genome GENOME3 with Bowtie2 (3/8)
[2018-02-28 15:29:04] Mapping right_kept_reads_seg4 to genome GENOME3 with Bowtie2 (4/8)
[2018-02-28 15:29:10] Mapping right_kept_reads_seg5 to genome GENOME3 with Bowtie2 (5/8)
[2018-02-28 15:29:17] Mapping right_kept_reads_seg6 to genome GENOME3 with Bowtie2 (6/8)
[2018-02-28 15:29:26] Mapping right_kept_reads_seg7 to genome GENOME3 with Bowtie2 (7/8)
[2018-02-28 15:29:32] Mapping right_kept_reads_seg8 to genome GENOME3 with Bowtie2 (8/8)
[2018-02-28 15:29:39] Searching for junctions via segment mapping
[2018-02-28 15:33:10] Retrieving sequences for splices
[2018-02-28 15:33:12] Indexing splices
Building a SMALL index
[2018-02-28 15:33:14] Mapping left_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 15:33:17] Mapping left_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 15:33:21] Mapping left_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 15:33:24] Mapping left_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 15:33:27] Mapping left_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 15:33:30] Mapping left_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 15:33:34] Mapping left_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 15:33:37] Mapping left_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 15:33:40] Joining segment hits
[2018-02-28 15:34:13] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/8)
[2018-02-28 15:34:16] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/8)
[2018-02-28 15:34:19] Mapping right_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/8)
[2018-02-28 15:34:23] Mapping right_kept_reads_seg4 to genome segment_juncs with Bowtie2 (4/8)
[2018-02-28 15:34:26] Mapping right_kept_reads_seg5 to genome segment_juncs with Bowtie2 (5/8)
[2018-02-28 15:34:30] Mapping right_kept_reads_seg6 to genome segment_juncs with Bowtie2 (6/8)
[2018-02-28 15:34:33] Mapping right_kept_reads_seg7 to genome segment_juncs with Bowtie2 (7/8)
[2018-02-28 15:34:36] Mapping right_kept_reads_seg8 to genome segment_juncs with Bowtie2 (8/8)
[2018-02-28 15:34:40] Joining segment hits
[2018-02-28 15:35:13] Reporting output tracks
-----------------------------------------------
[2018-02-28 15:35:56] A summary of the alignment counts can be found in /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat2/align_summary.txt
[2018-02-28 15:35:56] Run complete: 00:12:07 elapsed
  adding: gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx/tophat2/accepted_hits.bam (deflated 17%)
Left reads:
          Input     :    706769
           Mapped   :    580709 (82.2% of input)
            of these:     13823 ( 2.4%) have multiple alignments (829 have >20)
Right reads:
          Input     :    706769
           Mapped   :    590875 (83.6% of input)
            of these:     14086 ( 2.4%) have multiple alignments (829 have >20)
82.9% overall read mapping rate.

Aligned pairs:    579632
     of these:     13812 ( 2.4%) have multiple alignments
                  573628 (99.0%) are discordant alignments
 0.8% concordant pair alignment rate.
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
No count extraction for - pathogen
No count extraction for - host1
No count extraction for - host2
INFO: Calculations are stored under: /gpfs/scratch/pr74ma/ge34juq2/ge34juq2/tmp.Hdf7VwbZYx

```







