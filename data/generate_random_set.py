import sys

N=0
fh=file(sys.argv[1])
for line in fh.readlines():
	line=line.strip()
	N=N+1
N=(N/4)*int(sys.argv[2])/100.0
N=N*4
M=0
fh=file(sys.argv[1])
for line in fh.readlines():
	line=line.strip()
	if(M < N):
		print line
	else:
		sys.exit()
	M=M+1
