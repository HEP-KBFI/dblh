import datetime
start=datetime.datetime.now()

from input import *
flist=[]
#samples=samplespowheg
#samples.remove('stop')
for s in samples:
    for f in groups[s]['sets']:
        flist.append(f)

flist.append("data")
flist.append("ap-350")
flist.append("pp-350")

from extLibFull import filterByCutHDF,plotHDF,getModel
model=getModel('BP4')

from pylab import *
print "Starting processing...",
fnam="BP4-4l"
#fnam="BP4-3l"
varList=['mp']
nbin=[25]
xl=[r'$m(\ell^\pm\ell^\pm)$ [GeV]']
yl=[r'Events / 20 GeV']
xr=[[10,510]]
#xl=[r'$E^{miss}_{T}$ [GeV]']
#yl=[r'N events / 5 GeV']
#xr=[[0,100]]
ll=[1]

#cut1='~lowMassRes & (nlep>3) & (isoworst<0.35) & (sipworst<4) & (tWorstId>1)'
cut1='~lowMassRes & (nlep>3) & (ntau<3) & (isoworst<0.35) & (sipworst<4) & (tWorstId>1)' # 4-lep
#cut1='~lowMassRes & (nlep==3) & (ntau<2) & (isoworst<0.35) & (sipworst<4) & (tWorstId>1)' # 3-lep
#cut1='(nlep==3) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*275+125) & (dz>80) & (phi1<0.005*275+1.15) & (met>20)'
#cut1='(nlep==3) & (ntau==0) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*275+125) & (dz>80) & (phi1<0.005*275+1.15) & (met>20) &  ~lowMassRes & (tWorstId>1)'
#cut2='(nlep==3) & (ntau==1) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*275+125) & (dz>80) & (phi1<0.005*275+1.15) & (met>20) &  ~lowMassRes & (tWorstId>1)'
#cut1='recoPass & (nlep==3) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*130+125) & (dz>80) & (phi1<0.005*130+1.15) & (met>20) & ~lowMassRes'
#cut1='recoPass & (nlep==3) & (ntau==0) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*130+125) & (dz>80) & (phi1<0.005*130+1.15) & (met>20) &  ~lowMassRes'
#cut1='recoPass & (nlep==3) & (ntau==0) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*130+125) & (dz>80) & (phi1<0.005*130+1.15) & (met>20) & ~lowMassRes'
#cut2='recoPass & (nlep==3) & (ntau==1) & (isoworst<0.35) & (sipworst<4) & (sumpt>0.85*130+125) & (dz>80) & (phi1<0.005*130+1.15) & (met>20) &  ~lowMassRes'
#cut1='(nlep==3) & (ntau<2) & (isoworst<0.35) & (sipworst<4) & ~lowMassRes & recoPass & (tWorstId>1) & (sumpt>0.85*250+125) & (dz>80) & (phi1<0.005*250+1.15) & (met>20)'
#cut1='recoPass & (nlep==3) & (ntau==0) & (isoworst<0.35) & (sipworst<4)'
#cut1='recoPass & ~lowMassRes & (nlep==3) & (isoworst<0.35) & (sipworst<4) & (tWorstId>1) & (dz<30) & (sumpt>150)'
#cut1='recoPass & ~lowMassRes & (nlep>2) & (isoworst<0.35) & (sipworst<4) & (tWorstId>1)'
#cut1='recoPass & ~lowMassRes & (nlep==3) & (isoworst<0.35) & (sipworst<4) & (dz<80) & (sumpt>0.85*275+125) & (phi1<0.005*275+1.15) & (met>20)'
filt1=filterByCutHDF(flist,cut1,varList,model)
#filt2=filterByCutHDF(flist,cut2,varList,model)
#fnam="BP4-4l"
plotHDF(filt1,samples,"appp",350,getModel("BP4"),varList,xl,yl,nbin,xr,fnam+'.pdf',ll,showRatio=False)
savefig(fnam+'.png')
draw()
#fnam="BP4-3l-1tau-medTau"
#plotHDF(filt2,samples,"appp",275,getModel("BP4"),varList,xl,yl,nbin,xr,fnam+'.png')
#draw()
