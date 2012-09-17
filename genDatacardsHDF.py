from input import *
from string import Template
import datetime
from math import *
from extLibFull import getModel,getYieldHDF,filterByCutHDF
start=datetime.datetime.now()
flist=[]
dir="datacards/Fall11-lumi4.9-noDY/"
samples=samplesDC
for s in samples:
    for f in groups[s]['sets']:
        flist.append(f)
flist.append('data')
texLabel={
        'eeee':"100\% decay to electrons",
        'emem':"100\% decay to e-$\mu$ pairs",
        'mmmm':"100\% decay to muons",
        'etet':"100\% decay to e-$\\tau$ pairs",
        'mtmt':"100\% decay to $\mu\\tau$ pairs",
        'BP1':"benchmark point 1",
        'BP2':"benchmark point 2",
        'BP3':"benchmark point 3",
        'BP4':"benchmark point 4"
}

#unLumi=1.045
unLumi=1.022
unLep=0.02
unTau=0.06
unXsec=1.1

dcTemplate=Template("""Doubly charged higgs. Model=$model. Mass=$mass GeV
----------------------------------------------------------------------------------
imax 1
jmax 1
kmax 7
----------------------------------------------------------------------------------
Observation $obs
----------------------------------------------------------------------------------
bin 1 1
process 0 1
----------------------------------------------------------------------------------
rate $sig $bg
----------------------------------------------------------------------------------
lumi lnN $unLumi 1.00                           ! Luminosity uncertainty
lepton lnN $unLep 1.00                          ! lepton ID + isolation uncertainty
trigger lnN 1.015 1.015                         ! trigger+PV uncertainty
xs lnN $unXsec 1.00                             ! cross section uncertainty
cr$name lnN 1.00 $unCR                          ! stat uncert in CR
sigstat lnN $unSigStat 1.00                     ! Statistical uncertainty of the signal
r$name lnN 1.00 $unRatio                        ! Systematics on ratio (PDF, QCD scale, lepton energy scale)
""")

from cutsHDF import *
models=dict()
models['prompt']=["eeee", "emem", "mmmm" ]
models['ltau']=["etet","mtmt","BP1","BP2","BP3","BP4"]
modList=texLabel.keys()
import sys
if len(sys.argv)>1: 
    modList=sys.argv
    modList.remove(sys.argv[0])
print modList
for mod in modList:
    modstart=datetime.datetime.now()
    textout="""
    \\begin{table*}[p]
    \\begin{center}
    \\scriptsize{
    \\caption{Background estimation for """+texLabel[mod]+""".}
    \\label{tab:bgest-"""+mod+"""}
    \\begin{tabular}{|c|c|clc|c|c|c|}
    \\hline
    Mass & Final state & MC estimate & Estimate from data & Observation & Pair-production & Asso. production \\\\ \n"""
    for mass in masses:
        massStart=datetime.datetime.now()
        print "Starting:",mod,mass
        textout+="\\hline\n"
        tmax=[2,3]

        mmin=0.5*mass
        mmax=1.1*mass
        cMod="ltau"
        if mod in ['eeee','emem','mmmm']: 
            mmin=0.9*mass
            tmax=[1,1]
            cMod="prompt"

        nDataCR=[[0,0],[0,0,0]]
        nDataSR=[[0,0],[0,0,0]]
        nBGCR=[[0,0],[0,0,0]]
        nBGSR=[[0,0],[0,0,0]]
        nSigPP=[[0,0],[0,0,0]]
        nSigAP=[[0,0],[0,0,0]]
        flist2=flist+['ap-'+str(mass),'pp-'+str(mass)]
        
        # Let's start with 3 leptons
        for il in range(2):
            nl=il+3
            for nt in range(tmax[il]):
                print "il:",il,"nt:",nt
                mx=len(cuts[cMod][il])
                cutCR="(ntau=="+str(nt)+") & ((mp<"+str(mmin)+") | (mp>"+str(mmax)+"))"
                cutSR="(ntau=="+str(nt)+") & (mp>"+str(mmin)+") & (mp<"+str(mmax)+")"
                if il:
                    cutCR+=" & ((mn<"+str(mmin)+") | (mn>"+str(mmax)+"))"
                    cutSR+=" & (mn>"+str(mmin)+") & (mn<"+str(mmax)+")"
                
                if mod in extras.keys():
                    cutCR+=" & "+extras[mod][il]
                    cutSR+=" & "+extras[mod][il]
                
                cutCR+=" & "+cuts[cMod][il][0].substitute(mass=mass)
                for c in range(mx):
                    cutSR+=" & "+cuts[cMod][il][c].substitute(mass=mass)

                print "Debug:",il,nt,cutCR,cutSR
                crdata=filterByCutHDF(flist2,cutCR,["mp"],getModel(mod))
                srdata=filterByCutHDF(flist2,cutSR,["mp"],getModel(mod))
                yldCR=getYieldHDF(crdata,samples,getModel(mod),addErr=True)
                yldSR=getYieldHDF(srdata,samples,getModel(mod),addErr=True)

                nDataCR[il][nt]=yldCR['data']
                nDataSR[il][nt]=yldSR['data']

                nBGCR[il][nt]=yldCR['bg']
                nBGSR[il][nt]=yldSR['bg']
                #for s in samples:
                #   nBGCR[il][nt]+=yldCR[s]
                #   nBGSR[il][nt]+=yldSR[s]
           
                nSigPP[il][nt]+=yldSR['pp-'+str(mass)]
                nSigAP[il][nt]+=yldSR['ap-'+str(mass)]

                alpha=nBGSR[il][nt]*1./nBGCR[il][nt]
                dest=round(alpha*(nDataCR[il][nt]+1),2)
                unAlpha=yldSR['bg-err']/yldSR['bg']+0.05
                derr1=round(dest*unAlpha,3)
                derr2=round(1./sqrt(nDataCR[il][nt]+1),3)
                name='l'*(il+3-nt)+'t'*nt
                nameTeX='$'+'\ell'*(il+3-nt)+'\\tau'*nt+'$'
                lepUn=sqrt((nl-nt)*unLep**2+nt*unTau**2)+1
                unSigPP=yldSR['pp-'+str(mass)+'-err']/yldSR['pp-'+str(mass)]+1
                unSigAP=yldSR['pp-'+str(mass)+'-err']/yldSR['pp-'+str(mass)]+1
                unSigComb=sqrt((unSigPP-1)**2+(unSigAP-1)**2)+1
                dcPP=dcTemplate.substitute(model=mod,mass=mass,obs=nDataSR[il][nt],name=name,sig=nSigPP[il][nt],bg=dest,unLumi=unLumi,unLep=lepUn,unXsec=unXsec,unCR=derr2+1,unSigStat=unSigPP,unRatio=1+unAlpha)
                fout=open(dir+mod+'-m'+str(mass)+'-'+name+'-SP-PP.dat','wb')
                fout.write(dcPP)
                fout.close()
                dcAP=dcTemplate.substitute(model=mod,mass=mass,obs=nDataSR[il][nt],name=name,sig=nSigAP[il][nt],bg=dest,unLumi=unLumi,unLep=lepUn,unXsec=unXsec,unCR=derr2+1,unSigStat=unSigAP,unRatio=1.+unAlpha)
                fout=open(dir+mod+'-m'+str(mass)+'-'+name+'-SP-AP.dat','wb')
                fout.write(dcAP)
                fout.close()
                dcComb=dcTemplate.substitute(model=mod,mass=mass,obs=nDataSR[il][nt],name=name,sig=nSigPP[il][nt]+nSigAP[il][nt],bg=dest,unLumi=unLumi,unLep=lepUn,unXsec=unXsec,unCR=derr2+1,unSigStat=unSigComb,unRatio=1.+unAlpha)
                fout=open(dir+mod+'-m'+str(mass)+'-'+name+'-SP-Comb.dat','wb')
                fout.write(dcComb)
                fout.close()
                textout+=str(mass)+"~\GeV & "+nameTeX+" & $"+str(round(yldSR['bg'],2))+"\pm"+str(round(yldSR['bg-err'],3))+"$ & $"+str(dest)+"\pm"+str(derr1)+"\pm"+str(round(derr2*dest,3))+"$ & $"+str(nDataSR[il][nt])+"\pm"+str(round(sqrt(nDataSR[il][nt]),2))+"$ & $"+str(round(yldSR['pp-'+str(mass)],2))+"\pm"+str(round(yldSR['pp-'+str(mass)+'-err'],3))+"$ & $"+str(round(yldSR['ap-'+str(mass)],2))+"\pm"+str(round(yldSR['ap-'+str(mass)+'-err'],3))+"$ \\\\ \n"
        print "Mass running took:",str(datetime.datetime.now()-massStart)
    textout+="""
    \\hline
    \end{tabular}
    }
    \end{center}
    \end{table*}
    """
    ftex=open(dir+mod+'-bgcheck-SP.tex','wb')
    ftex.write(textout)
    ftex.close()
    print "Model took:",str(datetime.datetime.now()-modstart)

