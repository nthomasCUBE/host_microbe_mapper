import sys

# ----------------------------------------------------------------------------

DIR=sys.argv[1]
fh=file("bin/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta")
allowed={}
for line in fh.readlines():
	line=line.strip()
	vals=line.split()
	if(len(line)>0 and line[0]==">"):
		allowed[vals[0][1:]]=1

# ----------------------------------------------------------------------------

fh=file(DIR+"/merged.fastq_accepted.sam")
fh=file(sys.argv[1])
MAP={}
for line in fh.readlines()[1:]:
	line=line.strip()
	vals=line.split()
	if(MAP.get(vals[2])==None):
		MAP[vals[2]]={}
	if(MAP[vals[2]].get(vals[0])==None):
		MAP[vals[2]][vals[0]]={}
	MAP[vals[2]][vals[0]][vals[9]]=1

# ----------------------------------------------------------------------------

fw=file(DIR+"/extract_16S.fa","w")
for MAP_ in MAP:
	if(allowed.get(MAP_)!=None):
		gg=MAP[MAP_]
		for gg_ in gg:
			ee=MAP[MAP_][gg_]
			for ee_ in ee:
				fw.write(">"+MAP_+"_"+gg_+"\n")
				fw.write(ee_+"\n")
fw.close()

# ----------------------------------------------------------------------------

