from input import *
import datetime
from extLibFull import getModel,getYieldHDF,filterByCutHDF
start=datetime.datetime.now()
flist=[]
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

from cutsHDF import *
models=dict()
models['prompt']=["eeee", "emem", "mmmm" ]
models['ltau']=["etet","mtmt","BP1","BP2","BP3","BP4"]
modList=texLabel.keys()
import sys
if len(sys.argv)>1:
    modList=[sys.argv[1]]
for mod in modList:
    print "\\clearpage"
    for m in masses:
        mmin=0.5*m
        mmax=1.1*m
        cMod='ltau'
        if mod in ['eeee','emem','mmmm']: 
            mmin=0.9*m
            cMod='prompt'
        print """
        \\begin{table*}[htbp]
        \\begin{center}"""
        print "\caption{Cut flow for "+texLabel[mod]+" with mass "+str(m)+"~\GeV}"
        print "\label{tab:"+mod+"-"+str(m)+"}"
        print """\scriptsize{
        \\begin{tabular}{|c|c|c|c|c|c|c|c|c|}
        \hline"""
        print "Cut & ",
        for s in samples:
            print groups[s]['mplabel']+' & ',
        print 'Total & Data & Signal \\\\'
        print '\\hline'
        for il in [0,1]:
            nl=il+3
            cut="(mp>"+str(mmin)+") & (mp<"+str(mmax)+")"
            if il:
                cut+=" & (mn>"+str(mmin)+") & (mn<"+str(mmax)+")"
            if mod in extras.keys():
                cut+=" & "+extras[mod][il]
            for c in range(len(cuts[cMod][il])):
                cut+=" & "+cuts[cMod][il][c].substitute(mass=m)
                filt=filterByCutHDF(flist+['ap-'+str(m),'pp-'+str(m)],cut,["mp"],getModel(mod))
                yld=getYieldHDF(filt,samples,getModel(mod),addErr=False)
                print str(nl)+'l-'+str(c+1),'&',
                for s in samples:
                    print "$"+str(round(yld[s],2))+"\pm"+str(round(yld[s+'-err'],2))+'$ &',
                print "$"+str(round(yld['bg'],2))+'\pm'+str(round(yld['bg-err'],2))+'$ & $'+str(yld['data'])+'\pm'+str(round(yld['data-err'],2))+'$ & $'+str(round(yld['signal'],2))+'\pm'+str(round(yld['signal-err'],2))+'$ \\\\'
            print '\\hline'
        print """
        \end{tabular}
        }
        \end{center}
        \end{table*}

        """

