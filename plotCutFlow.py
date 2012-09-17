from input import *
flist=[]
#samples.remove('stop')
for s in samples:
    for f in groups[s]['sets']:
        flist.append(f)
flist.append('data')

from extLibFull import drawCutFlow
from pylab import show
from cutsHDF import *
from copy import deepcopy
import sys
from string import Template
from pylab import draw

models=dict()
models['prompt']=["eeee", "emem", "mmmm" ]
models['ltau']=["etet","mtmt","BP1","BP2","BP3","BP4"]

model=sys.argv[1]
mass=350
flist.append('ap-'+str(mass))
flist.append('pp-'+str(mass))
mmin=0.5*mass
mmax=1.1*mass
clist=[]
if model in models['prompt']: 
    clist=deepcopy(cuts['prompt'])
    mmin=0.9*mass
    clist[1].append(Template('(dz>0)'))
else: 
    clist=deepcopy(cuts['ltau'])

if model in extras.keys(): clist.append(extras[model])
clist[0].append(Template('(mp>'+str(mmin)+') & (mp<'+str(mmax)+')'))
clist[1].append(Template('(mp>'+str(mmin)+') & (mp<'+str(mmax)+') & (mn>'+str(mmin)+') & (mn<'+str(mmax)+')'))
drawCutFlow(flist,samples,clist,mass,model,model+'-cutflow-3l',[0])
draw()
samples.remove('stop')
drawCutFlow(flist,samples,clist,mass,model,model+'-cutflow-4l',[1])
