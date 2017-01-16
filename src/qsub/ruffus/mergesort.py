import os
import sys
import subprocess
import fnmatch

mergedBamfn = sys.argv[1]
finalbamdir = sys.argv[2]

command = ""
os.chdir(finalbamdir)
matches = []

for root,dirnames, filenames in os.walk(finalbamdir):
    for filename in fnmatch.filter(filenames, '*.bam'):
        
        path = os.path.join(root, filename)
        if os.path.islink(path):
            path = os.path.realpath(path)
            
        if (not matches.__contains__(path)):
            matches.append(path)
            command = " ".join([path, command])
        
command = " ".join(["/mnt/work1/software/sambamba/0.5.4/sambamba merge", mergedBamfn, command,"-t 4" ])
subprocess.check_output(command, shell = True)       