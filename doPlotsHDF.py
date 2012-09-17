from input import *
import datetime
from extLibFull import getModel,filterByCutHDF,plotHDF
from pylab import *
start=datetime.datetime.now()
modMass={'eeee':350, 'emem':350, 'etet':350, 'mmmm':350, 'mtmt':350, 'BP1':350, 'BP2':350, 'BP3':350, 'BP4':350}
flist=[]
for s in samples:
    for f in groups[s]['sets']:
        flist.append(f)
flist.append('data')

from cutsHDF import *
models=dict()
models['prompt']=["eeee", "emem", "mmmm" ]
models['ltau']=["etet","mtmt","BP1","BP2","BP3","BP4"]

xl3=r'$m(\ell^{\pm}\ell^{\pm})$ [GeV]'
xl4=[r'$m(\ell^{+}\ell^{+})$ [GeV]',r'$m(\ell^{-}\ell^{-})$ [GeV]']
yl=r'N events / 20 GeV'
yl4=[yl,yl]
var3="mp"
var4=["mp","mn"]
nb3=25
nb4=[nb3,nb3]
xr3=[10,510]
xr4=[xr3,xr3]

modList=modMass.keys()
import sys
if len(sys.argv)>1: 
        modList=[sys.argv[1]]
for mod in modList:
    m=modMass[mod]
    mmin=0.5*m
    mmax=1.1*m
    cMod='ltau'
    if mod in ['eeee','emem','mmmm']: 
        mmin=0.9*m
        cMod='prompt'
    for il in [0,1]:
        nl=il+3
        cut="(nlep>2)"
        if mod in extras.keys():
            cut+=" & "+extras[mod][il]
        for c in range(len(cuts[cMod][il])):
            cut+=" & "+cuts[cMod][il][c].substitute(mass=m)
            vlist=var3
            x=xl3
            y=yl
            nb=nb3
            xr=xr3
            ll=1
            lloc=1
            if c: lloc=3
            sc='linear'
            vlistFilt=["mp"]
            if nl==4:
                vlist=var4
                vlistFilt=var4
                x=xl4
                y=yl4
                nb=nb4
                xr=xr4
                ll=[1,1]
                lloc=[1,1]
                if c: lloc=[3,3]
                sc=['linear','linear']
            print il,c,cut
            filt=filterByCutHDF(flist+['ap-'+str(m),'pp-'+str(m)],cut,vlistFilt,getModel(mod))
            fname=mod+'-m'+str(m)+'-'+str(nl)+'l-cut'+str(c)
            plotHDF(filt,samples,'appp',m,getModel(mod),vlist,x,y,nb,xr,fname+'.pdf',lloc)
            savefig(fname+'.png')
            draw()
            fname=mod+'-m'+str(m)+'-'+str(nl)+'l-cut'+str(c)+'-lin'
            plotHDF(filt,samples,'appp',m,getModel(mod),vlist,x,y,nb,xr,fname+'.pdf',ll,sc)
            savefig(fname+'.png')
