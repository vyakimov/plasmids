import sys
import os
cmd="bamtools merge"
for arg in sys.argv[1:len(sys.argv)-1]:
    cmd = cmd + " -in " + arg
cmd = cmd + " -out " + sys.argv[len(sys.argv)-1]
os.system(cmd)
