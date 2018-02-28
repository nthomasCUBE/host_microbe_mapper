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
![alt text](https://github.com/nthomasCUBE/host_pathogen_mapping/blob/master/misc/pix.png)

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



