from ROOT import *
from array import *
from extLibFull import convertTreeToTable,debug,createPUweight
from input import *
import datetime
from multiprocessing import Pool,cpu_count

print "Starting execution:",str(datetime.datetime.now())
ncpu=cpu_count()
createPUweight()

# create processing pool with available CPU's
cpool = Pool(ncpu)
base='inputData-v3/'
# get the list of files we have to process
flist=[]
for s in samplesConvert:
    for f in groups[s]['sets']:
        flist.append(f)
flist.append("data")
for f in ppN:
    flist.append(f)

for f in apN:
    flist.append(f)

import sys

if len(sys.argv)>1:
    print "Only doing:",sys.argv[1]
    convertTreeToTable(sys.argv[1])
else:
    cpool.map(convertTreeToTable,flist)
