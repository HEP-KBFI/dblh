import matplotlib as mpl
from types import *
mpl.use('QT4Agg')
params = {
            'axes.labelsize': 20,
            'text.fontsize': 20,
            'legend.fontsize': 14,
            'xtick.labelsize': 20,
            'ytick.labelsize': 20,
            'text.usetex': True
}
mpl.rcParams.update(params)

from multiprocessing import Pool,cpu_count
from pylab import *
from ROOT import *
from input import *
from array import *
from copy import copy,deepcopy
from math import *
import cPickle
from tables import *
import numpy
from pylab import draw,show

intList=["run","lumi","ev","np"]
#intArrList=["pid","status","id"]
intArrList=["pid","status","id","conv"]
doubleList=["npu","vx","vy","vz","metx","mety","rho"]
doubleArrList=["px", "py", "pz", "E", "M", "MC", "SIP", "Q", "TkIso", "EcalIso", "HcalIso" ]
allvars=intList+intArrList+doubleList+doubleArrList
ArrMax=100
EvArrMax=10
weight=TH1D("weight","weight",1000,0,50)
weightint=TH1D("weightint","weight",25,0,25)
debug=false
ngen=dict()
fstates=['data','eeee','eeem','eeet','eemm', 'eemt', 'eett', 'emem','emet','emmm','emmt','emtt','etet','etmm','etmt','ettt','mmmm','mmmt','mmtt','mtmt','mttt','tttt']
unLumi=0.045
unXsec=0.05
unLep=0.02
unTau=0.06

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

cutLabel=[['Pre-selection', r'$\sum p_\mathrm{T}$', 'Z veto', r'$\Delta\varphi(\ell\ell)$','Mass window'], ['Pre-selection',r'$\sum p_\mathrm{T}$','Z veto','Mass window']]



fall11=[
        0.003388501,
        0.010357558,
        0.024724258,
        0.042348605,
        0.058279812,
        0.068851751,
        0.072914824,
        0.071579609,
        0.066811668,
        0.060672356,
        0.054528356,
        0.04919354,
        0.044886042,
        0.041341896,
        0.0384679,
        0.035871463,
        0.03341952,
        0.030915649,
        0.028395374,
        0.025798107,
        0.023237445,
        0.020602754,
        0.0180688,
        0.015559693,
        0.013211063,
        0.010964293,
        0.008920993,
        0.007080504,
        0.005499239,
        0.004187022,
        0.003096474,
        0.002237361,
        0.001566428,
        0.001074149,
        0.000721755,
        0.000470838,
        0.00030268,
        0.000184665,
        0.000112883,
        6.74043E-05,
        3.82178E-05,
        2.22847E-05,
        1.20933E-05,
        6.96173E-06,
        3.4689E-06,
        1.96172E-06,
        8.49283E-07,
        5.02393E-07,
        2.15311E-07,
        9.56938E-08
  ]

class Event(IsDescription):
    recoPass    = BoolCol()
    lowMassRes  = BoolCol()
    candle      = StringCol(2)
    nlep        = UInt8Col()
    nele        = UInt8Col()
    nmu         = UInt8Col()
    ntau        = UInt8Col()
    tWorstId    = UInt8Col()
    fstate      = UInt8Col()
    nbs         = UInt8Col()
    run         = UInt16Col()
    lumi        = UInt16Col()
    ev          = UInt16Col()
    pid         = Int32Col(shape=(EvArrMax,))
    sumpt       = Float64Col()
    dz          = Float64Col()
    mp          = Float64Col()
    mn          = Float64Col()
    hpt1        = Float64Col()
    hpt2        = Float64Col()
    phi1        = Float64Col()
    phi2        = Float64Col()
    met         = Float64Col()
    isoworst    = Float64Col()
    sipworst    = Float64Col()
    puWeight    = Float64Col()
    lepWeight   = Float64Col()
    iso         = Float64Col(shape=(EvArrMax,))
    sip         = Float64Col(shape=(EvArrMax,))

class Candle(IsDescription):
    run         = UInt16Col()
    lumi        = UInt16Col()
    ev          = UInt16Col()
    nlep        = UInt16Col()
    fstate      = UInt16Col()
    pid1        = Int32Col()
    pid2        = Int32Col()
    mp          = Float64Col()
    pt1         = Float64Col()
    pt2         = Float64Col()
    eta1        = Float64Col()
    eta2        = Float64Col()
    iso1        = Float64Col()
    iso2        = Float64Col()
    sip1        = Float64Col()
    sip2        = Float64Col()
    puWeight    = Float64Col()
    lepWeight   = Float64Col()
    phi         = Float64Col()

class Leptons(IsDescription):
    nlep        = UInt16Col()
    fstate      = UInt16Col()
    tid         = UInt16Col()
    pid         = Int32Col()
    pt          = Float64Col()
    eta         = Float64Col()
    puWeight    = Float64Col()
    lepWeight   = Float64Col()
    iso         = Float64Col()
    sip         = Float64Col()
    isoworst    = Float64Col()
    sipworst    = Float64Col()

def AxNum(x,pos):
    #print x,pos,log(x,10),str(int(log(x,10))),str(int(x))
    if x<1 or x>=100:
        return r'$10^{'+str(int(round(log(x,10))))+'}$'
    if x>=1 and x<100:
        return r'$'+str(int(x))+'$'

def initTree(varlist,tree):
    """
    Deactivates all other branches, creates the relevant variables and associates them to the tree
    """
    tree.SetBranchStatus("*",0)
    cacheSize=10000000
    tree.SetCacheSize(cacheSize)
    for v in varlist:
            tree.SetBranchStatus(v,1)

    vList=dict([(v,makeVar(v)) for v in varlist])
    if debug: print "Debug: ",vList
    for v in varlist:
        tree.SetBranchAddress(v,vList[v])
        tree.AddBranchToCache(v,kTRUE)

    tree.StopCacheLearningPhase()
    return vList

def readInput(flist,level=0):
    """
    Reads the input into a dictionary
    """
    res=dict()
    for f in flist:
        inf=open(base+f+'-l'+str(level)+'.dat','rb')
        res[f]=cPickle.load(inf)

    return res

def makeVar(vName):
    """
    Returns the proper variable type
    """
    if vName in intList: return array('i',[0])
    if vName in intArrList: return array('i',[0]*ArrMax)
    if vName in doubleList: return array('d',[0])
    if vName in doubleArrList: return array('d',[0]*ArrMax)

def createPUweight():
    """
    Creates the weight list
    """
    pufile="pileup-2011-int.root"
    f=TFile(pufile)
    pu=f.Get("pileup")
    pu.Scale(1./pu.Integral(0,25))
    f11tot=0
    for i in range(25):
        f11tot+=fall11[i]
    f11scale=1./f11tot
    for i in range(25):
        weightint.Fill(i,pu.GetBinContent(i+1)/(f11scale*fall11[i]))

def createPUweightOld():
    """
    Creates the weight list
    """
    pufile="pileup-2011.root"
    f=TFile(pufile)
    pu=f.Get("pileup")
    pu.Scale(1./pu.Integral(0,1001))
    pu.Rebin(2)
    fall11=TFile("fall2011exp.root")
    fexp=fall11.Get("F2011exp")
    fexp.Scale(1./fexp.Integral(0,500))
    for i in range(502):
        #print "i=",i,"data=",pu.GetBinLowEdge(i),pu.GetBinWidth(i),"F11=",fexp.GetBinLowEdge(i),fexp.GetBinWidth(i),
        if fexp.GetBinContent(i):
            weight.SetBinContent(i,pu.GetBinContent(i)/fexp.GetBinContent(i))
        else:
            weight.SetBinContent(i,0)
        #print weight.GetBinContent(i)

def getPUweight(v,data=false):
    """
    Calculates the specific event weight from pileup
    """
    w=1.
    if not data:
        w=weightint.GetBinContent(weightint.FindBin(v['npu'][0]))
    v['puWeight']=w
    return v

def getEleWeights(v,data=false):
    """
    calculates the weight from electron T&P
    """
    elist=v['esel']
    weights= [
        [0.997, 0.988, 1.036, 1.028], # pT->15-20
        [0.987, 1.009, 1.012, 1.007], # pT->20-30
        [0.995, 0.997, 1.009, 1.004], # pT->30-40
        [0.994, 0.997, 1.003, 1.002], # pT->40-50
        [0.994, 0.996, 1.006, 1.002], # pT->50-60
        [0.997, 0.998, 1.006, 0.999] # pT->60+
    ]
    w=[1.]*v['np'][0]
    # Loop over the electrons
    if not data:
        for i in elist:
            pt=v['lvs'][i].Pt()
            eta=abs(v['lvs'][i].Eta())
            ptb=-1
            etb=-1
            # now determine which bin to use
            if pt<20: ptb=0
            if pt>=20 and pt<30: ptb=1
            if pt>=30 and pt<40: ptb=2
            if pt>=40 and pt<50: ptb=3
            if pt>=50 and pt<60: ptb=4
            if pt>=60: ptb=5
            if eta<0.78: etb=0
            if eta>=0.78 and eta<1.4442: etb=1
            if eta>=1.566 and eta<2: etb=2
            if eta>=2 and eta<2.5: etb=3
            w[i]=weights[ptb][etb]

    v['eWeight']=w
    return v

def getMuWeights(v,data=false):
    """
    calculates weight for muon T&P
    """
    w=[1.]*v['np'][0]
    v['mWeight']=w
    return v 

def getTauWeights(v,data=false):
    """
    calculates weight for tau T&P
    """
    w=[1.]*v['np'][0]
    v['tWeight']=w
    return v 

def addLVs(v):
    """
    This function constructs the TLorentzVectors and adds them to the dictionary for this event
    """
    lvs=[]
    for i in range(v['np'][0]):
        lvs.append(TLorentzVector(v['px'][i],v['py'][i],v['pz'][i],v['E'][i]))

    v["lvs"]=lvs
    return v

def addLVsHDF(v):
    lvs=[]
    for i in range(v['np']):
        lv=v['lv'][i]
        lvs.append(TLorentzVector(lv[0],lv[1],lv[2],lv[3]))
    v['lvs']=lvs
    return v

def selectRecoElectrons(v):
    """
    This function loops over all particles in the event and selects electrons, that pass our requirements
    """
    global debug
    passList=[]
    eidReq=4
    eptMin=15
    for i in range(v['np'][0]):
        if abs(v['pid'][i]) == 11 and v['status'][i]:
            # Well at least it's a reconstructed electron
            id=v['id'][i]
            conv=v['conv'][i] 
            ch=abs(v['Q'][i])
            pt=v['lvs'][i].Pt()
            eta=abs(v['lvs'][i].Eta())
            if debug: print "Found electron, particle nr:",i
            if id>=eidReq and conv<2 and ch==1 and pt>eptMin and eta<2.5:
                passList.append(i)

    v['esel']=passList
    return v

def selectRecoMuons(v):
    """
    This function loops over all particles in the event and selects muons, that pass our requirements
    """
    global debug
    passList=[]
    midReq=1+4+8+16+32+64#+128
    muMinPt=5
    for i in range(v['np'][0]):
        if abs(v['pid'][i]) == 13 and v['status'][i]:
            # It's a reco muon
            id=v['id'][i]
            pt=v['lvs'][i].Pt()
            eta=abs(v['lvs'][i].Eta())
            if debug: print "Found muon, particle nr:",i
            if id & midReq==midReq and pt>muMinPt and eta<2.4:
                passList.append(i)
    
    v['msel']=passList
    return v

def selectRecoTaus(v):
    """
    Thif function loops over all particles in the event and selects tau jets that pass requirements
    """
    global debug
    passList=[]
    tauMinPt=15
    wID=4
    tauIdMin=1+32+128+512 # Loose dBeta iso + Tight mu + medium ele + electron MVA
    for i in range(v['np'][0]):
        if abs(v['pid'][i]) == 15 and v['status'][i]:
            # it's a reco tau jet
            id=v['id'][i]
            pt=v['lvs'][i].Pt()
            eta=abs(v['lvs'][i].Eta())
            if debug: print "Found tau, particle nr:",i
            if id & tauIdMin == tauIdMin and pt>tauMinPt and eta<2.1 and (eta<1.46 or eta>1.558):
                # now check overlaps
                overlap=false
                for l in v['esel']+v['msel']:
                    if v['lvs'][i].DeltaR(v['lvs'][l])<0.1:
                        overlap=true

                if not overlap:
                    passList.append(i)
                    if id & 1 == 1 and id & 2 == 0 and wID>1: wID=1
                    if id & 2 == 2 and id & 4 == 0 and wID>2: wID=2
                    if id & 4 == 4 and wID>3: wID=3

    v['tsel']=passList
    v['tWorstId']=wID
    return v

def tagBs(v):
    """
    This tags any b quarks in the event
    """
    b=0
    for i in range(v['np'][0]):
        if not v['status'][i] and abs(v['pid'][i])==5:
            b+=1
    #if b: print "Found b's:",b
    v['nbs']=b
    return v

def removeQuarkonia(v):
    """
    This function removes low mass resonances (m<12 GeV, no flavor requirement)
    """
    totlep=len(v['esel']+v['msel']+v['tsel'])
    selLep=v['esel']+v['msel']
    nl=len(selLep)
    v['lowMassRes']=False
    for i in range(nl-1):
        for j in range(i+1,nl):
            il=selLep[i]
            jl=selLep[j]
            l1=v['lvs'][il]
            l2=v['lvs'][jl]
            #print il,jl,v['pid'][il],v['pid'][jl],(l1+l2).M()
            #if totlep>3: print i,j,v['pid'][i],v['pid'][j],(l1+l2).M()
            if (l1+l2).M() < 12:
                v['lowMassRes']=True
                # This shit is quarkonia!!!
                #if i in v['esel']: v['esel'].remove(i)
                #if j in v['esel']: v['esel'].remove(j)
                #if i in v['msel']: v['msel'].remove(i)
                #if j in v['msel']: v['msel'].remove(j)
    return v

def doBasicStuff(v,data=false):
    """
    This is a helper function that does the usual stuff per event
    in one go without having to think about it. I.e. selects leptons
    calculates weights etc...
    """
    v=addLVs(v)
    v=selectRecoElectrons(v)
    v=selectRecoMuons(v)
    v=selectRecoTaus(v)
    v=removeQuarkonia(v)
    v=tagBs(v)
    v=getEleWeights(v,data)
    v=getMuWeights(v,data)
    v=getTauWeights(v,data)
    v=getPUweight(v,data)
    v=calcIso(v)
    v=checkIsoSIP(v)
    v=calcAdditionalVariables(v)
    v=getFinalState(v,data)
    return v

def calcIso(v):
    """
    This function calculates the relative isolation of every lepton and
    does the effective area PU removal. The result is stored back
    """
    efEl=[0.078, 0.046, 0.026, 0.072];
    efMu=[0.087, 0.049, 0.042, 0.059];
    rho=v['rho'][0]
    np=v['np'][0]
    isoList=[0]*np
    for p in range(np):
        pid=abs(v['pid'][p])
        iso=0
        if pid==11 or pid==13:
            tk=v['TkIso'][p]
            ecal=v['EcalIso'][p]
            hcal=v['HcalIso'][p]
            pt=v['lvs'][p].Pt()
            eta=abs(v['lvs'][p].Eta())
            eRho=0;
            hRho=0;
            if pid==11 and eta<1.479:
                eRho=efEl[0]*rho
                hRho=efEl[2]*rho

            if pid==11 and eta>1.479:
                eRho=efEl[1]*rho
                hRho=efEl[3]*rho

            if pid==13 and eta<1.479:
                eRho=efMu[0]*rho
                hRho=efMu[2]*rho

            if pid==13 and eta>1.479:
                eRho=efMu[1]*rho
                hRho=efMu[3]*rho
            iso=(tk+ecal+hcal-eRho-hRho)/pt
            isoList[p]=iso

    v['iso']=isoList
    return v

def getModel(mod):
    """
    Returns the 21 branching ratios
    """
    model=[]
    if not isinstance(mod,ListType):
        if mod=="eeee": 
            mod=[1.,0,0,0,0,0]
        elif mod=="emem":
            mod=[0,1.,0,0,0,0]
        elif mod=="etet":
            mod=[0,0,1.,0,0,0]
        elif mod=="mmmm":
            mod=[0,0,0,1.,0,0]
        elif mod=="mtmt":
            mod=[0,0,0,0,1.,0]
        elif mod=="tttt":
            mod=[0,0,0,0,0,1.]
        elif mod=="BP1":
            mod=[0.0032,0.0064,0.0064,0.3017,0.3807,0.3017]
        elif mod=="BP2":
            mod=[0.4972,0,0,0.1256,0.2513,0.1256]
        elif mod=="BP3":
            mod=[0.34,0,0,0.33,0,0.33]
        elif mod=="BP4":
            mod=[1./6]*6
        else:
            mod=[0]*6

    for i in range(6):
        for j in range(i,6):
            f=2
            if i==j: f=1
            model.append(f*mod[i]*mod[j])

    return model

def getFinalState(v,data=false,debug=False):
    """
    This function will define the final state based on MC information
    """
    # First off let's filter out generator level objects and wether or not it's even signal
    if data:
        v['fstate']=0
        return v

    genList=[]
    isSig=false
    isAP=false
    if debug: print "From:",v['run'],v['lumi'],v['ev']
    for p in range(v['np'][0]):
        if debug: print "Inp:",p,v['pid'][p],v['Q'][p],v['status'][p]
        if v['status'][p]==0:
            if abs(v['pid'][p]) > 10 and abs(v['pid'][p])<17: genList.append(p)

        if abs(v['pid'][p])==37: isAP=true
        if abs(v['pid'][p])==38: isAP=true
        if abs(v['pid'][p])==9900041 or isAP: isSig=true

    if not isSig:
        # ok, it's not signal, but it's MC so we can define a few additional final states ... but later :)
        v['fstate']=99
        return v

    posList=[]
    negList=[]
    nuList=[]
    for p in genList:
        if v['Q'][p]>0: posList.append(abs(v['pid'][p]))
        if v['Q'][p]<0: negList.append(abs(v['pid'][p]))
        if not v['Q'][p]: nuList.append(abs(v['pid'][p]))

    pos=['t' if x==15 else x for x in ['m' if x==13 else x for x in ['e' if x==11 else x for x in posList]]]
    neg=['t' if x==15 else x for x in ['m' if x==13 else x for x in ['e' if x==11 else x for x in negList]]]
    nu=['t' if x==16 else x for x in ['m' if x==14 else x for x in ['e' if x==12 else x for x in nuList]]]
    if debug: print "Debug: ",pos,neg,nu
    if len(pos+neg+nu)!=4:
        v['recoPass']=false
        return v

    if isAP:
        if len(pos)==1: pos.append(nu[0])
        if len(neg)==1: neg.append(nu[0])
    pos.sort()
    neg.sort()
    sel=[pos,neg]
    sel.sort()
    sel=[x for a in sel for x in a]
    fs=''
    for x in sel:
        fs+=x

    if len(sel)==4:
        if debug: print "Final:",fs,fstates.index(fs)
        v['fstate']=fstates.index(fs)
    else:
        v['recoPass']=false
    return v

def fillNgen(fn,model=0):
    if 'data' in fn: return
    global ngen
    if fn not in ngen.keys():
        sk=TFile(base+fn+"_skim.root")
        if not model:
            ngen[fn]=sk.Get("trig/eff").GetBinContent(1)
        else:
            br=sk.Get("branch/br")
            ngen[fn]=[0]*21
            for j in range(21):
                ix=0
                iy=0
                if (j < 7):
                    ix=2
                    iy=j+2
                if (j>5 and j<11):
                    ix=3
                    iy=j-6+3
                if (j>10 and j<15):
                    ix=4
                    iy=j-6-5+4
                if (j>14 and j<18):
                    ix=5
                    iy=j-6-5-4+5
                if (j>17 and j<20):
                    ix=6
                    iy=j-6-5-4-3+6
                if (j==20):
                    ix=7
                    iy=7
                fs=j
                ng=br.GetBinContent(ix,iy)+br.GetBinContent(iy,ix)
                if ix==iy:
                    ng=ng/2
                ngen[fn][fs]=ng
    return

def getScaleFactor(v,fn,model=0):
    """
    Get the scale factor sigma*lumi/Ngenerated for backgrounds and for signal it takes into account
    the final state and model
    """
    if 'data' in fn: return 1.
    global ngen
    w=1.
    fillNgen(fn,model)
    # ok, now let's split the math to signal vs background
    if not model:
        w=datasets[fn]['xsec']*lumi/ngen[fn]
    else:
        fs=v['fstate']-1
        tp=fn[0:2]
        mass=int(fn[3:])
        xsec=0
        if tp=="ap":
            xsec=xsAP[masses.index(mass)]
        if tp=="pp":
            xsec=xsPP[masses.index(mass)]

        #print "SF debug:",fs,xsec,lumi,fn,ngen[fn],model[fs]
        w=model[fs]*xsec*lumi/ngen[fn][fs]

    return w

def checkIsoSIP(v,isoC=0.35,sipC=4,label='isosipcheck'):
    """
    Check if the event passes isolation/SIP requirements
    """
    ok=false
    elist=v['esel']
    mlist=v['msel']
    list=elist+mlist
    iso1=0
    iso2=0
    sip1=0
    for p in list:
        iso=v['iso'][p]
        sip=v['SIP'][p]
        if iso>iso2 and iso<iso1:
            iso2=iso

        if iso>iso1:
            iso2=iso1
            iso1=iso
        
        if sip>sip1:
            sip1=sip

    if iso1+iso2<isoC and sip1<sipC:
        ok=true

    v[label]=ok
    v['isoworst']=iso1+iso2
    v['sipworst']=sip1
    return v

def getPhi(f1,f2):
    """
    Calculates the difference in opening angle
    """
    pi=3.1415
    if f1<0: f1+=2*pi
    if f2<0: f2+=2*pi
    dphi=abs(f1-f2)
    if dphi>pi: dphi=2*pi-dphi
    return dphi

def getMt(l1,l2):
    """
    Calculates transverse mass between two lorentz vectors
    """
    a1=(l1.Pt()+l2.Pt())**2
    b1=(l1+l2).Pt()**2
    return sqrt(a1+b1)

def calcAdditionalVariables(v):
    """
    Adds variables for final analysis like H++ mass, pT sum etc
    """
    sel=v['esel']+v['msel']+v['tsel']
    nsel=len(sel)
    npos=0
    nneg=0
    reco=false
    mp=0
    mn=0
    phi1=0
    phi2=0
    hpt1=0
    hpt2=0
    sumpt=0
    dz=9999
    lepUsed=[]
    mx=v['metx'][0]
    my=v['mety'][0]
    metV=TLorentzVector(mx,my,0,0)
    met=metV.Pt()
    if 'candle' in v.keys(): del v['candle']
    if nsel>2:
        for p in sel:
            if v['Q'][p]>0: npos+=1
            if v['Q'][p]<0: nneg+=1
        for i in range(nsel-1):
            for j in range(i+1,nsel):
                ip=sel[i]
                jp=sel[j]
                if v['pid'][ip]==-v['pid'][jp]:
                    m=(v['lvs'][ip]+v['lvs'][jp]).M()
                    if abs(m-91.2)<dz: dz=abs(m-91.2)

    if nsel==2:
        l1=v['lvs'][sel[0]]
        l2=v['lvs'][sel[1]]
        sumpt=l1.Pt()+l2.Pt()
        mp=(l1+l2).M()
        hpt1=(l1+l2).Pt()
        phi1=getPhi(l1.Phi(),l2.Phi())
        fs='xx'
        p1=v['pid'][sel[0]]
        p2=v['pid'][sel[1]]
        if p1==-p2 and v['iso'][sel[0]]<0.2 and v['iso'][sel[1]]<0.2 and l1.Pt()>20 and l2.Pt()>20:
            if abs(p1)==11: fs='ee'
            if abs(p1)==13: fs='mm'
        if fs != 'xx': 
            v['candle']=fs  
        reco=false

    if nsel==3:
        # check charge topology
        if (nneg==1 and npos==2) or (nneg==2 and npos==1):
            # Ok, we actually have the right topology
            for p in sel:
                sumpt+=v['lvs'][p].Pt()
                lepUsed.append(p)

            for i in range(2):
                for j in range(i+1,3):
                    ip=sel[i]
                    jp=sel[j]
                    if v['Q'][ip]*v['Q'][jp] > 0:
                        # found our H++ candidate
                        phi1=getPhi(v['lvs'][ip].Phi(),v['lvs'][jp].Phi())
                        mp=(v['lvs'][ip]+v['lvs'][jp]).M()
                        k=list(sel)
                        k.remove(ip)
                        k.remove(jp)
                        kp=k[0]
                        mn=getMt(v['lvs'][kp],metV)
                        hpt1=(v['lvs'][ip]+v['lvs'][jp]).Pt()
                        hpt2=(v['lvs'][kp]+metV).Pt()
                        reco=true

    if nsel>3:
        # check charge topology
        if nneg>1 and npos>1:
            # Ok, we have at least one right topology
            posCands=[]
            negCands=[]
            # loop over all combinations and find all pairs
            for i in range(nsel-1):
                for j in range(i+1,nsel):
                    ip=sel[i]
                    jp=sel[j]
                    if v['Q'][ip]*v['Q'][jp] > 0:
                        if v['Q'][ip]>0: posCands.append([ip,jp])
                        else: negCands.append([ip,jp])

            # now loop over all candidates to find closest in mass
            mdif=9999
            bestPos=-1
            bestNeg=-1
            for i in range(len(posCands)):
                for j in range(len(negCands)):
                    pmass=(v['lvs'][posCands[i][0]]+v['lvs'][posCands[i][1]]).M()
                    nmass=(v['lvs'][negCands[j][0]]+v['lvs'][negCands[j][1]]).M()
                    if abs(pmass-nmass)<mdif:
                        mdif=abs(pmass-nmass)
                        bestPos=i
                        bestNeg=j

            lepList=posCands[bestPos]+negCands[bestNeg]
            lepUsed=copy(lepList)
            for i in lepList:
                sumpt+=v['lvs'][i].Pt()

            mp=(v['lvs'][posCands[bestPos][0]]+v['lvs'][posCands[bestPos][1]]).M()
            mn=(v['lvs'][negCands[bestNeg][0]]+v['lvs'][negCands[bestNeg][1]]).M()
            hpt1=(v['lvs'][posCands[bestPos][0]]+v['lvs'][posCands[bestPos][1]]).Pt()
            hpt2=(v['lvs'][negCands[bestNeg][0]]+v['lvs'][negCands[bestNeg][1]]).Pt()
            phi1=getPhi(v['lvs'][posCands[bestPos][0]].Phi(),v['lvs'][posCands[bestPos][1]].Phi())
            phi2=getPhi(v['lvs'][negCands[bestNeg][0]].Phi(),v['lvs'][negCands[bestNeg][1]].Phi())
            reco=true
    
    v['recoPass']=reco
    w=1.
    for p in lepUsed:
        if abs(v['pid'][p])==11: w=w*v['eWeight'][p]
        if abs(v['pid'][p])==13: w=w*v['mWeight'][p]
        if abs(v['pid'][p])==15: w=w*v['tWeight'][p]

    v['sumpt']=sumpt
    v['met']=met
    v['mp']=mp
    v['mn']=mn
    v['hpt1']=hpt1
    v['hpt2']=hpt2
    v['phi1']=phi1
    v['phi2']=phi2
    v['dz']=dz
    v['lepWeight']=w

    return v

def convertTree(f,level=0):
    # get input data
    tFile=TFile(base+f+"_tree.root")
    tree=tFile.Get("ntuple/cands")
    sf=1.
    data=true
    if not "data" in f:
        data=false

    print "Starting:",f
    # init Tree
    v=initTree(allvars,tree)
    ntot=tree.GetEntries()
    vList=[]
    #ntot=100
    getList=['sumpt','dz','mp','mn','hpt1','hpt2','phi1','phi2','met', 'isoworst', 'sipworst', 'puWeight', 'lepWeight', 'tWorstId','fstate','nbs', 'lowMassRes']
    for i in range(ntot):
        if i % 50000 == 0: print f,i
        tree.LoadTree(i)
        tree.GetEntry(i)
        v=doBasicStuff(v,data)
        ev=dict()
        if not v['recoPass']: continue
        ev['nele']=len(v['esel'])
        ev['nmu']=len(v['msel'])
        ev['ntau']=len(v['tsel'])
        if level:
            ev['pid']=[]
            ev['pt']=[]
            ev['eta']=[]
            for p in v['esel']+v['msel']+v['tsel']:
                ev['pid'].append(v['pid'][p])
                ev['pt'].append(v['lvs'][p].Pt())
                ev['eta'].append(v['lvs'][p].Eta())

            if level>1:
                ev['iso']=[]
                ev['sip']=[]
                ev['id']=[]
                for p in v['esel']+v['msel']+v['tsel']:
                    ev['iso'].append(v['iso'][p])
                    ev['sip'].append(v['SIP'][p])
                    ev['id'].append(v['id'][p])

        for var in getList:
            if var in v.keys():
                ev[var]=v[var]

        vList.append(copy(ev))

    return [f,vList]

def convertTreeToTable(f):
    # get input data
    tFile=TFile(base+f+"_tree.root")
    tree=tFile.Get("ntuple/cands")
    tabFile=openFile(base+f+".h5",mode='w')
    tab=tabFile.createTable('/','events',Event,"Event listing",expectedrows=10**6)
    cand=tabFile.createTable('/','candle',Candle,"Candle cands",expectedrows=10**7)
    lepton=tabFile.createTable('/','leptons',Leptons,"Leptons",expectedrows=10**7)
    ev=tab.row
    can=cand.row
    lep=lepton.row
    sf=1.
    data=true
    if not "data" in f:
        data=false

    print "Starting:",f
    # init Tree
    v=initTree(allvars,tree)
    ntot=tree.GetEntries()
    vList=[]
    #ntot=10
    getList=['sumpt','dz','mp','mn','hpt1','hpt2','phi1','phi2','met', 'isoworst', 'sipworst', 'puWeight', 'lepWeight', 'tWorstId','fstate','nbs','recoPass','lowMassRes','candle']
    for i in range(ntot):
        if i % 50000 == 0: print f,i
        tree.LoadTree(i)
        tree.GetEntry(i)
        #if v['np'][0] > 49: 
        #    print "Ran out of particle slots"
        #    continue
        v=doBasicStuff(v,data)
        #if not v['recoPass']: continue
        nele=len(v['esel'])
        nmu=len(v['msel'])
        ntau=len(v['tsel'])
        nlep=nele+nmu+ntau
        ev['nele']=nele
        ev['nmu']=nmu
        ev['ntau']=ntau
        ev['nlep']=nlep
        lsel=v['esel']+v['msel']+v['tsel']
        ev['pid']=numpy.array([v['pid'][i] for i in lsel]+[0]*(EvArrMax-nlep))
        ev['iso']=numpy.array([v['iso'][i] for i in lsel]+[0]*(EvArrMax-nlep))
        ev['sip']=numpy.array([v['SIP'][i] for i in lsel]+[0]*(EvArrMax-nlep))

        for var in getList:
            if var in v.keys():
                ev[var]=v[var]

        for var in ['run','lumi','ev']:
            ev[var]=v[var][0]

        ev.append()

        for i in range(nlep-1):
            for j in range(i+1,nlep):
                x1=lsel[i]
                x2=lsel[j]
                l1=v['lvs'][x1]
                l2=v['lvs'][x2]
                can['run']=v['run'][0]
                can['fstate']=v['fstate']
                can['lumi']=v['lumi'][0]
                can['ev']=v['ev'][0]
                can['nlep']=nlep
                can['pid1']=v['pid'][x1]
                can['pid2']=v['pid'][x2]
                can['mp']=(l1+l2).M()
                can['phi']=0
                can['pt1']=l1.Pt()
                can['pt2']=l2.Pt()
                can['eta1']=l1.Eta()
                can['eta2']=l2.Eta()
                can['iso1']=v['iso'][x1]
                can['iso2']=v['iso'][x2]
                can['sip1']=v['SIP'][x1]
                can['sip2']=v['SIP'][x2]
                can['puWeight']=v['puWeight']
                can['lepWeight']=v['lepWeight']
                can.append()

        for i in lsel:
            if abs(v['pid'][i])==15: lep['tid']=v['tWorstId']
            else: lep['tid']=0
            lep['pid']=v['pid'][i]
            lep['pt']=v['lvs'][i].Pt()
            lep['eta']=v['lvs'][i].Eta()
            lep['nlep']=len(lsel)
            lep['fstate']=v['fstate']
            lep['puWeight']=v['puWeight']
            lep['lepWeight']=v['lepWeight']
            lep['iso']=v['iso'][i]
            lep['sip']=v['SIP'][i]
            lep['isoworst']=v['isoworst']
            lep['sipworst']=v['sipworst']
            lep.append()

    tab.flush()
    tab.cols.nlep.createIndex()
    tab.cols.ntau.createIndex()
    tab.cols.isoworst.createIndex()
    tab.cols.sipworst.createIndex()
    tab.cols.sumpt.createIndex()
    tab.cols.dz.createIndex()
    tab.cols.mp.createIndex()
    cand.flush()
    cand.cols.mp.createIndex()
    cand.cols.nlep.createIndex()
    cand.cols.pid1.createIndex()
    cand.cols.pid2.createIndex()
    lepton.flush()
    lepton.cols.nlep.createIndex()
    lepton.cols.pid.createIndex()
    tabFile.close()
    return

def getLeptonPtEtaDistributions(f,histos,pureweigh=true):
    """
    Creates the lepton pT and Eta distributions
    """
    # Get input data
    tFile=TFile(base+f+"_tree.root")
    tree=tFile.Get("ntuple/cands")
    sf=1.
    data=true
    if not "data" in f:
        sFile=TFile(base+f+"_skim.root")
        ngen=sFile.Get("trig/eff").GetBinContent(1)
        xs=datasets[f]['xsec']
        # calculate scale factor
        sf=lumi*xs/ngen
        data=false

    # init Tree
    v=initTree(allvars,tree)
    print f,sf
    puW=0
    nopuW=0
    for i in range(tree.GetEntries()):
        tree.LoadTree(i)
        tree.GetEntry(i)
        v=doBasicStuff(v,data)
        if not v['isosipcheck']: continue
        w=1.
        if pureweigh:
            w=v['puWeight']
        nopuW+=1
        puW+=w
        if debug: print f,w,pureweigh
        if len(v['esel']+v['msel'])==2:
            for p in v['esel']:
                histos[0].Fill(v['lvs'][p].Pt(),v['eWeight'][p]*w)
                histos[3].Fill(v['lvs'][p].Eta(),v['eWeight'][p]*w)

            for p in v['msel']:
                histos[1].Fill(v['lvs'][p].Pt(),v['mWeight'][p]*w)
                histos[4].Fill(v['lvs'][p].Eta(),v['mWeight'][p]*w)

        if len(['tsel'])==1:
            for p in v['tsel']:
                histos[2].Fill(v['lvs'][p].Pt(),v['tWeight'][p]*w)
                histos[5].Fill(v['lvs'][p].Eta(),v['tWeight'][p]*w)

    for i in histos:
        i.Scale(sf*nopuW/puW)

    print "Done with ",f,"pt totals:",histos[0].Integral(),histos[1].Integral(),histos[2].Integral()
    print "Reading",tFile.GetBytesRead(),"bytes in",tFile.GetReadCalls(),"transactions"
    return [f,histos]

def getStandardCandles(f,histos,pureweigh=true):
    """
    Draw the standard candles
    """
    # get input data
    tFile=TFile(base+f+"_tree.root")
    tree=tFile.Get("ntuple/cands")
    sf=1.
    data=true
    nopuW=0
    puW=0
    if not "data" in f:
        sFile=TFile(base+f+"_skim.root")
        ngen=sFile.Get("trig/eff").GetBinContent(1)
        xs=datasets[f]['xsec']
        sf=lumi*xs/ngen
        data=false

    # init Tree
    v=initTree(allvars,tree)
    for i in range(tree.GetEntries()):
        tree.LoadTree(i)
        tree.GetEntry(i)
        v=doBasicStuff(v,data)
        if not v['isosipcheck']: continue
        w=1.
        if pureweigh:
            w=v['puWeight']
        puW+=w
        nopuW+=1
        es=v['esel']
        ms=v['msel']
        if len(es+ms)==2:
            # Possible 2l stuff
            if len(es)==2 and v['pid'][es[0]]==-v['pid'][es[1]]:
                # E+E- final state
                histos[0].Fill((v['lvs'][es[0]]+v['lvs'][es[1]]).M(),w)
                if v['iso'][es[0]]+v['iso'][es[1]]<0.15 and v['SIP'][es[0]]<2 and v['SIP'][es[1]]<2:
                    histos[2].Fill((v['lvs'][es[0]]+v['lvs'][es[1]]).M(),w)

            if len(ms)==2 and v['pid'][ms[0]]==-v['pid'][ms[1]]:
                # Mu+Mu- final state
                histos[1].Fill((v['lvs'][ms[0]]+v['lvs'][ms[1]]).M(),w)
                if v['iso'][ms[0]]+v['iso'][ms[1]]<0.15 and v['SIP'][ms[0]]<2 and v['SIP'][ms[1]]<2:
                    histos[3].Fill((v['lvs'][ms[0]]+v['lvs'][ms[1]]).M(),w)

    
    for i in histos:
        i.Scale(sf*nopuW/puW)

    print "Residual PU difference:",1.*nopuW/puW
    print "Done with ",f,"masses totals:",histos[0].Integral(),histos[1].Integral(),
    print " Reading",tFile.GetBytesRead(),"bytes in",tFile.GetReadCalls(),"transactions"
    return [f,histos]

def plotStackAndData(inD,samples,xleg,yleg,nbin,xr,ll=1):
    var=inD[0]
    weight=inD[1]
    flist=[]
    md=[]
    wd=[]
    cl=[]
    for s in samples:
        m=[]
        w=[]
        for f in groups[s]['sets']:
            flist.append(f)
            skim=TFile(base+f+"_skim.root")
            ngen=skim.Get("trig/eff").GetBinContent(1)
            xs=datasets[f]['xsec']
            sf=xs*lumi/ngen
            m+=var[f]
            w+=weight[f]

        md.append(m)
        wd.append(w)
        cl.append(groups[s]['mpcol'])

    fig=figure(figsize=[8,8],facecolor='w',dpi=75)
    fig.add_axes([0.15,0.4,0.8,0.5],xlim=xr,yscale='log',xlabel=xleg,ylabel=yleg)
    b=[0.02]*nbin
    for i in range(len(md)):
        a=hist(md[i],bins=nbin,range=xr,histtype='bar',lw=0,rwidth=1,weights=wd[i],color=cl[i],bottom=b,label=groups[samples[i]]['mplabel'])
        b=[b[i]+a[0][i] for i in range(nbin)]

    dh=hist(var['data'],bins=nbin,range=xr,histtype='step',lw=0)
    bw=(xr[1]-xr[0])/nbin
    xran=[xr[0]+(i+0.5)*bw for i in range(nbin)]
    errorbar(xran,dh[0],xerr=bw/2,yerr=[sqrt(i)-0.021 for i in dh[0]],fmt='ko',label='Data')
    mx=3*max(dh[0])
    legend(loc=ll)
    gca().set_ylim(0.02,mx)
    gca().set_xlim(xr[0],xr[1])
    fig.add_axes([0.15,0.05,0.8,0.25],xlim=xr)
    r=[dh[0][i]/b[i] for i in range(nbin)]
    re=[1./b[i]*(sqrt(a[0][i])+r[i]*sqrt(b[i])) for i in range(nbin)]
    errorbar(xran,r,xerr=bw/2,yerr=re,fmt='ko')
    gca().set_ylim(0.01,3)
    axhline(y=1)
    figtext(0.1,0.95,r"CMS Preliminary $\sqrt s=7$ TeV, $\int\mathcal{L}$="+str(round(lumi*1./1000,1))+" fb$^{-1}$",fontsize=24)

def plotStackDataSignal(inD,samples,sig,mass,model,vrN,xleg,yleg,nbin,xrn,fName,ll=0,sc=0):
    vars=inD[0]
    is4l=False
    if isinstance(vrN,ListType):
        is4l=True
    nvar=1
    if is4l: nvar=len(vrN)
    if not ll:
        if is4l:
            ll=[1 for i in range(nvar)]
        else:
            ll=1
    if not sc:
        if is4l:
            sc=['log' for i in range(nvar)]
        else:
            sc='log'
    weight=inD[1]
    fig=figure(figsize=[8*nvar,8],facecolor='w',dpi=75)
    figtext(0.1,0.95,r"CMS Preliminary $\sqrt s=7$ TeV, $\int\mathcal{L}$="+str(round(lumi*1./1000,1))+" fb$^{-1}$",fontsize=24)
    for n in range(nvar):
        k=''
        if is4l: k=vrN[n]
        md=[]
        wd=[]
        cl=[]
        sampUsed=[]
        for s in samples:
            m=[]
            w=[]
            for f in groups[s]['sets']:
                if is4l: m+=[i[k] for i in vars[f]]
                else: m+=vars[f]
                w+=weight[f]
            if len(m):
                md.append(m)
                wd.append(w)
                cl.append(groups[s]['mpcol'])
                sampUsed.append(s)

        slist=[]
        if "ap" in sig: 
            slist.append('ap-'+str(mass))
        if "pp" in sig: 
            slist.append('pp-'+str(mass))
        smd=[]
        swd=[]
        for s in slist:
            if is4l: smd+=[i[k] for i in vars[s]]
            else: smd+=vars[s]
            swd+=weight[s]
        axstart=0.15
        axw=0.8/nvar
        xl=xleg
        yl=yleg
        xr=xrn
        nb=nbin
        scl=sc
        lloc=ll
        if is4l:
            axstart=0.15/nvar+n*(axw+0.1)
            xl=xleg[n]
            yl=yleg[n]
            xr=xrn[n]
            nb=nbin[n]
            lloc=ll[n]
            scl=sc[n]
        
        fig.add_axes([axstart,0.4,axw,0.5],xlim=xr,yscale=scl,xlabel=xl,ylabel=yl)
        b=[0.02]*nb
        a=[[0]*nb]
        for i in range(len(md)):
            a=hist(md[i],bins=nb,range=xr,histtype='bar',lw=0,rwidth=1,weights=wd[i],color=cl[i],bottom=b,label=groups[sampUsed[i]]['mplabel'])
            b=[b[i]+a[0][i] for i in range(nb)]
        vd=vars['data']
        if is4l: vd=[i[k] for i in vars['data']]
        sh=[0]
        if len(slist): sh=hist(smd,bins=nb,range=xr,histtype='step',lw=2,rwidth=1,weights=swd,color='red',label='Signal',bottom=[0.02]*nb)
        dh=[0]
        if len(vd): dh=hist(vd,bins=nb,range=xr,histtype='step',lw=0)
        bw=(xr[1]-xr[0])*1./nb
        xran=[xr[0]+(i+0.5)*bw for i in range(nb)]
        if len(vd): errorbar(xran,dh[0],xerr=bw/2,yerr=[sqrt(i)-0.021 for i in dh[0]],fmt='ko',label='Data')
        fact=3
        if scl=='linear':
            fact=1.1
        mx=fact*max(dh[0]+sh[0])
        legend(loc=lloc)
        gca().set_ylim(0.02,mx)
        gca().set_xlim(xr[0],xr[1])
        fig.add_axes([axstart,0.05,axw,0.25],xlim=xr)
        r=[0]*nb
        re=[0]*nb
        if len(vd):
            r=[dh[0][i]/b[i] for i in range(nb)]
            re=[1./b[i]*(sqrt(a[0][i])+r[i]*sqrt(b[i])) for i in range(nb)]
        rmax=max(r)+max(re)
        rmin=min(r)-max(re)
        if rmax>3: rmax=3
        if rmin<0.1: rmin=0.1
        if rmin>1: rmin=0.9
        errorbar(xran,r,xerr=bw/2,yerr=re,fmt='ko')
        gca().set_ylim(rmin,rmax)
        axhline(y=1)
    savefig(fName)
    return fig

def plotHDF(vars,samples,sig,mass,model,vrN,xleg,yleg,nbin,xrn,fName,ll=0,sc=0,showRatio=True):
    is4l=False
    if isinstance(vrN,ListType):
        is4l=True
    nvar=1
    if is4l: nvar=len(vrN)
    if not ll:
        if is4l:
            ll=[1 for i in range(nvar)]
        else:
            ll=1
    if not sc:
        if is4l:
            sc=['log' for i in range(nvar)]
        else:
            sc='log'
    weight=dict()
    for f in vars.keys():
        weight[f]=[x['weight'] for x in vars[f]]
    vsize=8
    yStart=0.4
    ySize=0.5
    if not showRatio: 
        vsize=6
        yStart=0.1
        ySize=0.8
    fig=figure(figsize=[8*nvar,vsize],facecolor='w',dpi=75)
    #figtext(0.1,0.95,r"CMS Preliminary $\sqrt s=7$ TeV, $\int\mathcal{L}$="+str(round(lumi*1./1000,1))+" fb$^{-1}$",fontsize=24)
    figtext(0.15,0.95,r"CMS $\sqrt s=7$ TeV, $\int\mathcal{L}\textrm{dt}$ = "+str(round(lumi*1./1000,1))+" fb$^{-1}$",fontsize=24)
    for n in range(nvar):
        k=''
        if is4l: k=vrN[n]
        else: k=vrN
        md=[]
        wd=[]
        cl=[]
        sampUsed=[]
        for s in samples:
            m=[]
            w=[]
            for f in groups[s]['sets']:
                m+=[i[k] for i in vars[f]]
                w+=weight[f]
            if len(m):
                md.append(m)
                wd.append(w)
                cl.append(groups[s]['mpcol'])
                sampUsed.append(s)

        slist=[]
        if "ap" in sig: 
            slist.append('ap-'+str(mass))
        if "pp" in sig: 
            slist.append('pp-'+str(mass))
        smd=[]
        swd=[]
        for s in slist:
            smd+=[i[k] for i in vars[s]]
            swd+=weight[s]
        axstart=0.15
        axw=0.8/nvar
        xl=xleg
        yl=yleg
        xr=xrn
        nb=nbin
        scl=sc
        lloc=ll
        if is4l:
            axstart=0.15/nvar+n*(axw+0.1)
            xl=xleg[n]
            yl=yleg[n]
            xr=xrn[n]
            nb=nbin[n]
            lloc=ll[n]
            scl=sc[n]
        
        fig.add_axes([axstart,yStart,axw,ySize],xlim=xr,yscale=scl,xlabel=xl,ylabel=yl)
        gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(AxNum))
        b=[0.02]*nb
        a=[[0]*nb]
        for i in range(len(md)):
            a=hist(md[i],bins=nb,range=xr,histtype='bar',lw=0,rwidth=1,weights=wd[i],color=cl[i],bottom=b,label=groups[sampUsed[i]]['mplabel'])
            b=[b[i]+a[0][i] for i in range(nb)]
        vd=[i[k] for i in vars['data']]
        sh=[0]
        if len(slist): sh=hist(smd,bins=nb,range=xr,histtype='step',lw=2,rwidth=1,weights=swd,color='red',label='Signal ('+str(mass)+' GeV)',bottom=[0.02]*nb)
        dh=[0]
        if len(vd): dh=hist(vd,bins=nb,range=xr,histtype='step',lw=0)
        bw=(xr[1]-xr[0])*1./nb
        xran=[xr[0]+(i+0.5)*bw for i in range(nb)]
        yerrL=[]
        yerrH=[]
        for i in dh[0]:
            if i: 
                yerrL.append(-0.5+sqrt(i+0.25)-0.02)
                yerrH.append(0.5+sqrt(i+0.25)-0.02)
            else:
                yerrL.append(0)
                yerrH.append(0)
        if len(vd): errorbar(xran,dh[0],xerr=bw/2,yerr=[yerrL,yerrH],fmt='ko',label='Data')
        #if len(vd): errorbar(xran,dh[0],xerr=bw/2,yerr=[[i ? -0.5+sqrt(i+0.25)-0.02 : 0 for i in dh[0]],[i ? 0.5+sqrt(i+0.25)-0.02 : 0 for i in dh[0]]],fmt='ko',label='Data')
        fact=3
        if scl=='linear':
            fact=1.1
        mx=fact*max(dh[0]+sh[0])
        handles, labels = gca().get_legend_handles_labels()
        hand,lab=[],[]
        #dhand=Line2D((0,1),(0,1),marker='o',color='k')
        #hand.append(dhand)
        #lab.append("Data")
        i=labels.index("Data")
        hand.append(handles[i])
        lab.append(labels[i])
        smp=deepcopy(samples)
        smp.reverse()
        for s in smp:
            if groups[s]['mplabel'] in labels:
                i=labels.index(groups[s]['mplabel'])
                hand.append(handles[i])
                lab.append(labels[i])
        sline=Line2D((0,1),(0,1),linewidth=2,color='r')
        hand.append(sline)
        lab.append('Signal ('+str(mass)+' GeV)')
        #l=legend(loc=lloc,numpoints=1,scatterpoints=1)
        l=legend(hand,lab,loc=lloc,numpoints=1)
        l.draw_frame(False)
        gca().set_ylim(0.02,mx)
        gca().set_xlim(xr[0],xr[1])
        if showRatio:
            fig.add_axes([axstart,0.05,axw,0.25],xlim=xr)
            r=[0]*nb
            re=[0]*nb
            if len(vd):
                r=[dh[0][i]/b[i] for i in range(nb)]
                re=[1./b[i]*(sqrt(a[0][i])+r[i]*sqrt(b[i])) for i in range(nb)]
            rmax=max(r)+max(re)
            rmin=min(r)-max(re)
            if rmax>3: rmax=3
            if rmin<0.1: rmin=0.1
            if rmin>1: rmin=0.9
            errorbar(xran,r,xerr=bw/2,yerr=re,fmt='ko')
            gca().set_ylim(rmin,rmax)
            axhline(y=1)
    savefig(fName)
    return fig

def doFilter(f,inp,cutstr,var,mod):
    r=[]
    ww=[]
    model=0
    if mod and ('ap-' in f or 'pp-' in f): model=mod
    print "Processing:",f,"mod=",model
    for ev in inp[f]:
        w=1.
        if 'puWeight' in ev.keys(): w=w*ev['puWeight']
        if 'lepWeight' in ev.keys(): w=w*ev['lepWeight']
        w=w*getScaleFactor(ev,f,model)
        passed=True
        ev['nlep']=ev['nele']+ev['nmu']+ev['ntau']
        for c in cuts:
            if c.count('=='):
                vr=c.split('==')[0]
                va=eval(c.split('==')[1])
                if ev[vr]!=va: passed=False
            if c.count('>'):
                vr=c.split('>')[0]
                va=eval(c.split('>')[1])
                if ev[vr]<va: passed=False
            if c.count('<'):
                vr=c.split('<')[0]
                va=eval(c.split('<')[1])
                if ev[vr]>va: passed=False
        if passed:
            r.append(ev[var])
            ww.append(w)

    freturn [f,r,ww]

def convertCut(cut):
    res=[]
    ct=cut.split(' && ')
    for c in ct:
        if '==' in c:
            spl=c.split('==')
            res.append({'var':spl[0],'sign':0,'val':eval(spl[1])})
        if '<' in c:
            spl=c.split('<')
            res.append({'var':spl[0],'sign':-1,'val':eval(spl[1])})
        if '>' in c:
            spl=c.split('>')
            res.append({'var':spl[0],'sign':1,'val':eval(spl[1])})
    return res

def filterByCut(flist,inp,cutstring,varList,mod=0):
    """
    This function is supposed to parse the cut string and return the filtered list
    """
    r=dict()
    weight=dict()
    #print "Got to filtering",cutstring,mod,var,len(flist)
    cut=convertCut(cutstring)
    for f in flist:
        r[f]=[]
        weight[f]=[]
        model=0
        if mod and ('ap-' in f or 'pp-' in f): model=mod
        print "Processing:",f,"mod=",model
        for ev in inp[f]:
            w=1.
            if 'puWeight' in ev.keys(): w=w*ev['puWeight']
            if 'lepWeight' in ev.keys(): w=w*ev['lepWeight']
            w=w*getScaleFactor(ev,f,model)
            passed=True
            ev['nlep']=ev['nele']+ev['nmu']+ev['ntau']
            for c in cut:
                var=c['var']
                sgn=c['sign']
                val=c['val']
                if sgn==-1:
                    if ev[var]>=val: 
                        passed=False
                        break
                if sgn==0:
                    if ev[var]!=val:
                        passed=False
                        break
                if sgn==1:
                    if ev[var]<=val:
                        passed=False
                        break
            if passed:
                if isinstance(varList,ListType):
                    resp=dict()
                    for v in varList:
                        resp[v]=ev[v]
                    r[f].append(resp)
                else:
                    r[f].append(ev[varList])
                weight[f].append(w)

    return [r,weight]

def rowToDict(ev,vlist,f,model,**kw):
    r=dict()
    for v in vlist:
        r[v]=ev[v]
    #print 'Debug in converter:',ev['puWeight'],ev['lepWeight'],getScaleFactor(ev,f,model)
    w=1.
    #w=w*ev['puWeight']
    w=w*ev['lepWeight']
    w=w*getScaleFactor(ev,f,model)
    r['weight']=w
    r['fstate']=ev['fstate']
    if 'filtFun' in kw.keys():
        rd=kw['filtFun'](ev)
        for k in rd.keys():
            r[k]=rd[k]
    return r

def filterByCutHDF(flist,cutstring,varList,mod=0,**kw):
    """
    This function is supposed to parse the cut string and return the filtered list
    """
    r=dict()
    debug=False
    if 'debug' in kw.keys():
        debug=kw['debug']
    for f in flist:
        r[f]=[]
        model=0
        if mod and ('ap-' in f or 'pp-' in f): model=mod
        if debug: print "Processing:",f
        inf=openFile(base+f+'.h5')
        tab=inf.root.events
        if 'doCandle' in kw.keys():
            tab=inf.root.candle
        if 'doLeptons' in kw.keys():
            tab=inf.root.leptons
        if 'filtFun' in kw.keys():
            r[f]=[rowToDict(x,varList,f,model,filtFun=kw['filtFun']) for x in tab.where(cutstring)]
        else:
            r[f]=[rowToDict(x,varList,f,model) for x in tab.where(cutstring)]

    return r

def filterByIsoSIP(flist,cuts,inp,isoC,sipC):
    filt=dict()
    weight=dict()
    for c in cuts:
        filt[c]=dict()
        weight[c]=dict()
    for f in flist:
        for c in cuts:
            filt[c][f]=[]
            weight[c][f]=[]

        for ev in inp[f]:
            pid=ev['pid']
            pt=ev['pt']
            eta=ev['eta']
            nlep=len(pid)
            ntau=ev['ntau']

            if ev['isoworst']<isoC and ev['sipworst']<sipC:
                w=1.
                if 'puWeight' in ev.keys(): w=w*ev['puWeight']
                if 'lepWeight' in ev.keys(): w=w*ev['lepWeight']
                for p in range(len(pid)):
                    if abs(pid[p])==11:
                        filt['ept'][f].append(pt[p])
                        filt['eeta'][f].append(eta[p])
                    if abs(pid[p])==13:
                        filt['mpt'][f].append(pt[p])
                        filt['meta'][f].append(eta[p])
                    if abs(pid[p])==15:
                        filt['tpt'][f].append(pt[p])
                        filt['teta'][f].append(eta[p])
                    if abs(pid[p]) in [11,13,15]:
                        weight[c][f].append(w)

    return [filt,weight]

def doControlRegions(filt,flist,cuts,inp):
    for f in flist:
        for c in ['ttbar', 'DY', 'QCD1', 'QCD2']:
            filt[c][f]=[]

        for ev in inp[f]:
            nlep=len(pid)
            ntau=ev['ntau']
            iso=ev['iso']
            sip=ev['sip']
            iso.sort()
            sip.sort()

            if nlep==3:
                if iso[0]+iso[1]<0.35 and iso[2]>0.3 and sip[1]<4 and sip[2]>4 and ev['met']>30 and ev['sumpt']>150 and ev['dz']>20: filt['ttbar'][f].append(ev['mp'])
                if iso[1]+iso[2]<0.35 and sip[2]<4 and ev['sumpt']>150 and ev['dz']<30: filt['DY'][f].append(ev['mp'])
                if iso[1]+iso[2]>0.35 and sip[2]>4: filt['QCD1'][f].append(ev['mp'])
                if iso[1]+iso[2]>0.35 and sip[2]>4 and ev['met']<30 and ev['sumpt']>150 and ev['dz']>80: filt['QCD2'][f].append(ev['mp'])

    return filt

def getYield(filt,samples,model):
    yld=dict()
    bgtot=0
    bgtoterr=0
    global ngen
    for s in samples:
        yld[s]=0
        yld[s+'-err']=0
        for f in groups[s]['sets']:
            fillNgen(f)
            nscaled=sum(filt[1][f])
            yld[s]+=nscaled
            nev=len(filt[0][f])
            ng=ngen[f]
            if nev<3: 
                nev=3.
                nscaled=nev*lumi*datasets[f]['xsec']/ng
            un=(ng-nev)*1./(ng*nev)*(nscaled**2)
            yld[s+'-err']+=un
        yld[s+'-err']=sqrt(yld[s+'-err'])
        if yld[s]<yld[s+'-err']: yld[s]=yld[s+'-err']
        bgtot+=yld[s]
        bgtoterr+=yld[s+'-err']**2

    if 'data' in filt[1].keys():
        yld['data']=sum(filt[1]['data'])
        yld['data-err']=[-0.5+sqrt(yld['data']+0.25),0.5+sqrt(yld['data']+0.25)]

    slist=[i for i in filt[1].keys() if i[0:2]=='pp' or i[0:2]=='ap']
    sigtot=0
    siguntot=0
    tmp=dict()
    for s in slist:
        fillNgen(s,model)
        tmp[s]=dict()
        un=0
        for i in range(21): 
            if model[i]:
                tmp[s][i]=[]
        for i in range(len(filt[0][s])):
            fs=filt[0][s][i]['fstate']-1
            if fs in tmp[s].keys(): 
                tmp[s][fs].append(filt[1][s][i])
        for i in tmp[s].keys():
            nev=len(tmp[s][i])
            nscaled=sum(tmp[s][i])
            if nev < 3: nev=3
            ng=ngen[s][i]
            un+=((ng-nev)*1./(ng*nev))*(nscaled**2)
        siguntot+=un
        un=sqrt(un)

        yld[s]=sum(filt[1][s])
        yld[s+'-err']=un
        sigtot+=yld[s]

    yld['bg']=bgtot
    yld['bg-err']=[-0.5+sqrt(bgtoterr+0.25),0.5+sqrt(bgtoterr+0.25)]
    yld['signal']=sigtot
    yld['signal-err']=sqrt(siguntot)

    return yld

def getYieldHDF(filt,samples,model,addErr=True):
    yld=dict()
    bgtot=0
    bgtoterr=0
    global ngen
    for s in samples:
        yld[s]=0
        yld[s+'-err']=0
        for f in groups[s]['sets']:
            fillNgen(f)
            nscaled=sum([x['weight'] for x in filt[f]])
            yld[s]+=nscaled
            nev=len(filt[f])
            ng=ngen[f]
            #print s,f,yld[s],nev,ng,
            un=-1
            if nev<1 and s != "ZjetsDC": 
                nev=1.
                nscaled=nev*lumi*datasets[f]['xsec']/ng
            if nev>0:
                un=(ng-nev)*1./(ng*nev)*(nscaled**2)
            if un<0: un=0
            #print un
            if f not in ignoreErrors: yld[s+'-err']+=un
        yld[s+'-err']=sqrt(yld[s+'-err'])
        if yld[s]<yld[s+'-err'] and addErr: yld[s]=yld[s+'-err']
        bgtot+=yld[s]
        bgtoterr+=yld[s+'-err']**2

    if 'data' in filt.keys():
        yld['data']=sum([x['weight'] for x in filt['data']])
        yld['data-err']=[-0.5+sqrt(yld['data']+0.25),0.5+sqrt(yld['data']+0.25)]

    slist=[i for i in filt.keys() if i[0:2]=='pp' or i[0:2]=='ap']
    sigtot=0
    siguntot=0
    tmp=dict()
    for s in slist:
        fillNgen(s,model)
        tmp[s]=dict()
        un=0
        for i in range(21): 
            if model[i]:
                tmp[s][i]=[]
        for i in range(len(filt[s])):
            fs=filt[s][i]['fstate']-1
            if fs in tmp[s].keys(): 
                tmp[s][fs].append(filt[s][i]['weight'])
        for i in tmp[s].keys():
            nev=len(tmp[s][i])
            nscaled=sum(tmp[s][i])
            if nev < 3: nev=3
            ng=ngen[s][i]
            un+=((ng-nev)*1./(ng*nev))*(nscaled**2)
        if un<0: un=0
        siguntot+=un
        un=sqrt(un)

        yld[s]=sum([x['weight'] for x in filt[s]])
        yld[s+'-err']=un
        sigtot+=yld[s]

    yld['bg']=bgtot
    bgtoterr=sqrt(bgtoterr)
    yld['bg-err']=bgtoterr #[-0.5+sqrt(bgtoterr+0.25),0.5+sqrt(bgtoterr+0.25)]
    yld['signal']=sigtot
    yld['signal-err']=sqrt(siguntot)

    return yld

def drawCutFlow(flist,samples,cuts,mass,model,outFile,plot=[0,1]):
    """
    Do a nice cut flow plot
    """
    ylds={3:[],4:[]}
    for il in plot:
        cut=""
        print "Il:",il,len(cuts[il])
        if len(cuts)>2:
            cut=cuts[2][il]
        for c in range(len(cuts[il])):
            if cut: cut+=" & "+cuts[il][c].substitute(mass=mass)
            else: cut=cuts[il][c].substitute(mass=mass)
            nl=il+3
            filt=filterByCutHDF(flist+['ap-'+str(mass),'pp-'+str(mass)],cut,["mp"],getModel(model))
            yld=getYieldHDF(filt,samples,getModel(mod),addErr=True)
            ylds[nl].append(deepcopy(yld))
    nc3=len(ylds[3])
    cutRange3=range(nc3)
    bot=[0.02]*nc3
    fig=figure(figsize=[8*len(plot),6],facecolor='w',dpi=75)
    if 0 in plot:
        xmin=0.1
        xmax=0.4
        if len(plot)==1:
            xmin=0.15
            xmax=0.8
        fig.add_axes([xmin,0.1,xmax,0.8],yscale='log')
        gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(AxNum))
        for s in samples:
            w=[]
            for c in range(nc3):
                w.append(ylds[3][c][s])
            hist(cutRange3,bins=nc3,range=[-0.5,nc3-0.5],rwidth=1,lw=0,histtype='bar',weights=w,color=groups[s]['mpcol'],bottom=bot,label=groups[s]['mplabel'])
            for b in range(nc3):
                bot[b]+=w[b]
        w=[]
        for c in range(nc3):
            w.append(ylds[3][c]['signal'])
        hist(cutRange3,bins=nc3,range=[-0.5,nc3-0.5],rwidth=1,lw=2,histtype='step',weights=w,color='red',bottom=[0.02]*nc3,label='Signal('+str(mass)+' GeV)')
        yer=[[],[]]
        for c in range(nc3):
            derr=ylds[3][c]['data-err']
            dd=ylds[3][c]['data']
            print c,dd,derr
            if not dd: derr=[0,0]
            yer[1].append(derr[1])
            if dd-derr[0]<0.02: derr[0]=dd-0.021
            yer[0].append(derr[0])
        errorbar(cutRange3,[x['data'] for x in ylds[3]],yerr=yer,fmt='ko',label='Data')
        #scatter([1,1],[1,2],c='k',marker='o')
        w=[]
        bot=[0.02]*nc3
        cLab3=[]
        for c in range(nc3):
            cLab3.append(cutLabel[0][c])
            print c,bot[c],ylds[3][c]['bg-err'],ylds[3][c]['bg']
            w.append(2*ylds[3][c]['bg-err'])
            bot[c]+=ylds[3][c]['bg']-ylds[3][c]['bg-err']
        print "Debg:",nc3,cLab3
        hist(cutRange3,bins=nc3,range=[-0.5,nc3-0.5],rwidth=1,lw=0,histtype='bar',weights=w,bottom=bot,label='MC stat. uncert.',hatch='x',color='w',fill=False)
        gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(AxNum))
        gca().set_ylim(0.02,3*ylds[3][0]['data'])
        gca().set_xlim(-0.5,nc3-0.5)
        gca().set_xticks(range(nc3))
        gca().set_xticklabels(cLab3)
        gca().set_ylabel(r'Events')
        handles, labels = gca().get_legend_handles_labels()
        hand,lab=[],[]
        idata=labels.index('Data')
        print 'idata:',idata
        hand.append(handles[idata])
        lab.append("Data")
        smp=deepcopy(samples)
        smp.reverse()
        for s in smp:
            i=labels.index(groups[s]['mplabel'])
            hand.append(handles[i])
            lab.append(labels[i])
        sline=Line2D((0,1),(0,1),linewidth=2,color='r')
        hand.append(sline)
        lab.append('Signal ('+str(mass)+' GeV)')
        istatun=labels.index('MC stat. uncert.')
        hand.append(handles[istatun])
        lab.append('MC stat. uncert.')
        l=legend(hand,lab,numpoints=1,scatterpoints=1)
        l.draw_frame(False)

    if 1 in plot:
        xmin=0.55
        xmax=0.4
        if len(plot)==1:
            xmin=0.15
            xmax=0.8
        print len(plot),xmin,xmax
        fig.add_axes([xmin,0.1,xmax,0.8],yscale='log')
        nc4=len(ylds[4])
        cutRange4=range(nc4)
        print nc4
        bot=[0.02]*nc4
        for s in samples:
            w=[]
            for c in range(nc4):
                w.append(ylds[4][c][s])
            hist(cutRange4,bins=nc4,range=[-0.5,nc4-0.5],rwidth=1,lw=0,histtype='bar',weights=w,color=groups[s]['mpcol'],bottom=bot,label=groups[s]['mplabel'])
            for b in range(nc4):
                bot[b]+=w[b]
        w=[]
        for c in range(nc4):
            w.append(ylds[4][c]['signal'])
        hist(cutRange4,bins=nc4,range=[-0.5,nc4-0.5],rwidth=1,lw=2,histtype='step',weights=w,color='red',bottom=[0.02]*nc4,label='Signal ('+str(mass)+' GeV)')
        yer=[[],[]]
        for c in range(nc4):
            derr=ylds[4][c]['data-err']
            dd=ylds[4][c]['data']
            print c,dd,derr
            if not dd: derr=[0,0]
            yer[1].append(derr[1])
            if dd-derr[0]<0.02: derr[0]=dd-0.021
            yer[0].append(derr[0])
        errorbar(cutRange4,[x['data'] for x in ylds[4]],yerr=yer,fmt='ko',label='Data')
        w=[]
        bot=[0.02]*nc4
        cLab4=[]
        for c in range(nc4):
            cLab4.append('')
            cLab4.append(cutLabel[1][c])
            w.append(2*ylds[4][c]['bg-err'])
            bot[c]+=ylds[4][c]['bg']-ylds[4][c]['bg-err']
        hist(cutRange4,bins=nc4,range=[-0.5,nc4-0.5],rwidth=1,lw=0,histtype='bar',weights=w,bottom=bot,label='MC stat. uncert.',hatch='x',color='w',fill=False)
        l=1
        if nc4==3: l=3
        print nc4,l
        handles, labels = gca().get_legend_handles_labels()
        hand,lab=[],[]
        idata=labels.index('Data')
        print 'idata:',idata
        hand.append(handles[idata])
        lab.append("Data")
        smp=deepcopy(samples)
        smp.reverse()
        for s in smp:
            i=labels.index(groups[s]['mplabel'])
            hand.append(handles[i])
            lab.append(labels[i])
        sline=Line2D((0,1),(0,1),linewidth=2,color='r')
        hand.append(sline)
        lab.append('Signal ('+str(mass)+' GeV)')
        istatun=labels.index('MC stat. uncert.')
        hand.append(handles[istatun])
        lab.append('MC stat. uncert.')
        l=legend(hand,lab,numpoints=1,loc=l,scatterpoints=1)
        l.draw_frame(False)
        gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(AxNum))
        gca().set_ylim(0.02,3*ylds[4][0]['data'])
        gca().set_xticklabels(cLab4)
        gca().set_ylabel(r'Events')
    figtext(0.15,0.95,r"CMS $\sqrt s=7$ TeV, $\int\mathcal{L}\textrm{dt}$ = "+str(round(lumi*1./1000,1))+" fb$^{-1}$",fontsize=24)
    if len(plot)==2:
        figtext(0.55,0.95,r"Model: "+texLabel[model]+" with mass "+str(mass)+" GeV",fontsize=24)
        figtext(0.2,0.01,r'3 lepton final state',fontsize=24)
        figtext(0.65,0.01,r'4 lepton final state',fontsize=24)
    savefig(outFile+'.pdf')
    savefig(outFile+'.png')
    return ylds

