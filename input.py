datasets=dict()
datasets['Wbb']=dict(name='Wbb',xsec=116.04)
datasets['Ts']=dict(name='Ts',xsec=3.19)
datasets['Tt']=dict(name='Tt',xsec=41.92)
datasets['TtW']=dict(name='TtW',xsec=7.87)
datasets['Tbars']=dict(name='Tbars',xsec=1.44)
datasets['Tbart']=dict(name='Tbart',xsec=22.65)
datasets['TbartW']=dict(name='TbartW',xsec=7.87)
datasets['WWjets']=dict(name='WWjets',xsec=4.52565056512)
datasets['WZ3l']=dict(name='WZ3l',xsec=0.59428045677) 
datasets['WZ2l2q']=dict(name='WZ2l2q',xsec=1.23754954323)
datasets['ZZ4l']=dict(name='ZZ4l',xsec=0.0806888152125)
datasets['ZZ2l2nu']=dict(name='ZZ2l2nu',xsec=0.160672154445)
datasets['ZZ2l2q']=dict(name='ZZ2l2q',xsec=0.560316280343)
datasets['TT2l2nu2b']=dict(name='TT2l2nu2b',xsec=17.32)
datasets['DYee']=dict(name='DYee',xsec=1666.0)
datasets['DYmm']=dict(name='DYmm',xsec=1666.0)
datasets['DYtata']=dict(name='DYtata',xsec=1666.0)
datasets['DYmm10']=dict(name='DYmm10',xsec=3335)
datasets['DYee10']=dict(name='DYee10',xsec=3327)
datasets['DYll10']=dict(names='DYll10',xsec=12782.63)
datasets['DYll10_bEnr']=dict(names='DYll10_bEnr',xsec=12782.63*1.25)
datasets['DYll10_bSub']=dict(names='DYll10_bSub',xsec=12782.63)
datasets['DYll50']=dict(names='DYll50',xsec=3048.)
datasets['DYll50_bEnr']=dict(names='DYll50_bEnr',xsec=3048./0.984*1.25)
datasets['DYll50_bSub']=dict(names='DYll50_bSub',xsec=3048./0.984)
datasets['QCD_BC20']=dict(name='QCD_BC20',xsec=132160)
datasets['QCD_BC30']=dict(name='QCD_BC30',xsec=136804)
datasets['QCD_BC80']=dict(name='QCD_BC80',xsec=9360)
datasets['QCD_EM20']=dict(name='QCD_EM20',xsec=2454400)
datasets['QCD_EM30']=dict(name='QCD_EM30',xsec=3866200)
datasets['QCD_EM80']=dict(name='QCD_EM80',xsec=139500)
datasets['QCD_Mu15']=dict(name='QCD_Mu15',xsec=1668096.0)
datasets['QCD_Mu20']=dict(name='QCD_Mu20',xsec=349988.0)
datasets['QCD_b15']=dict(name='QCD_b15',xsec=39014400)
datasets['QCD_b30']=dict(name='QCD_b30',xsec=3557700)
datasets['QCD_b50']=dict(name='QCD_b50',xsec=599592)
datasets['QCD_d30']=dict(name='QCD_d30',xsec=10868.0)
datasets['QCD_d40']=dict(name='QCD_d40',xsec=43571.0)

groups=dict()
groups['Wbb']=dict(sets=[ 'Wbb', ],label='W+b#bar{b}',color=400,mpcol='y',mplabel=r'$W^\pm+b\bar b$')
groups['stop']=dict(sets=[ 'Ts', 'Tt', 'TtW', 'Tbars', 'Tbart', 'TbartW', ],label='Single top',color=629,mpcol='m',mplabel=r'Single top')
groups['VVjets']=dict(sets=[ 'WWjets', 'WZ3l', 'WZ2l2q', 'ZZ4l', 'ZZ2l2nu', 'ZZ2l2q', ],label='VV+jets',color=800,mpcol='orange',mplabel=r'Diboson')
groups['TT2l2b']=dict(sets=[ 'TT2l2nu2b', ],label='t#bar{t}+jets',color=600,mpcol='b',mplabel=r'$t\bar{t}$')
groups['ZjetsHM']=dict(sets=[ 'DYee', 'DYmm', 'DYtata', ],label='Z+jets',color=416,mpcol='g',mplabel=r'$Z^0$+jets')
groups['ZjetsLM']=dict(sets=[ 'DYmm10', 'DYee10', ],label='Z+jets (LM)',color=900,mpcol='0.25',mplabel=r'$Z^0$+jets (LM)')
groups['Zjets']=dict(sets=[ 'DYee', 'DYmm', 'DYtata', 'DYmm10', 'DYee10', ],label='Z+jets',color=416,mpcol='g',mplabel=r'Drell-Yan')
groups['ZjetsDC']=dict(sets=['DYll50'], label='Z+jets',color=416,mpcol='g',mplabel=r'$Z^0$+jets')
groups['ZjetsMG']=dict(sets=['DYll10','DYll50'], label='Z+jets',color=416,mpcol='g',mplabel=r'Drell-Yan')
groups['ZjetsLight']=dict(sets=['DYll10_bSub','DYll50_bSub'], label='Z+light jets',color=416,mpcol='g',mplabel=r'$Z^0$+u/d/c jets')
groups['ZjetsBB']=dict(sets=['DYll10_bEnr','DYll50_bEnr'], label='Z+bb',color=419,mpcol='cyan',mplabel=r'$Z^0+b\bar{b}$')

groups['QCD']=dict(sets=[ 'QCD_BC20', 'QCD_BC30', 'QCD_BC80', 'QCD_EM20', 'QCD_EM30', 'QCD_EM80', 'QCD_Mu15', 'QCD_Mu20', 'QCD_b15', 'QCD_b30', 'QCD_b50', 'QCD_d30', 'QCD_d40', ],label='QCD',color=920,mpcol='0.5',mplabel=r'QCD')
samplespowheg=[ "stop","TT2l2b","VVjets","Zjets"]
samplesZbb=["stop","TT2l2b","VVjets","ZjetsBB","ZjetsLight"]
samples=["stop","TT2l2b","VVjets","ZjetsMG"]
samplesDC=["TT2l2b","VVjets","ZjetsDC"]
samplesConvert=["Wbb", "stop","TT2l2b","VVjets","Zjets","ZjetsMG"]
ignoreErrors=['DYll10','DYll10_bSub','DYll10_bEnr']
#base="inputData-v3/"
#base="inputData-v3-nomuLayers/"
base="inputData-v4/"
#lumi=4634
lumi=4928

masses = [ 130, 150, 170, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500 ]
ppN = [ "pp-"+str(masses[i]) for i in range(len(masses))]
apN = [ "ap-"+str(masses[i]) for i in range(len(masses))]
from xsinfo import *
xsPP=[xsPPNLO[massPP.index(i)]/1000 for i in masses]
xsAP=[xsAPNLO[massAP.index(i)]/1000 for i in masses]
