from string import Template
cuts=dict()
cuts['prompt']=[
    [
        Template("(nlep==3) & (ntau==0) & (isoworst<0.35) & (sipworst<4) & ~lowMassRes & recoPass"),
        Template("(sumpt>1.1*$mass+60)"),
        Template("(dz>80)"),
        Template("(phi1<(1./600)*$mass+1.95)")
    ],[
        Template("(nlep>3) & (ntau==0) & (isoworst<0.35) & (sipworst<4) & ~lowMassRes & recoPass"),
        Template("(sumpt>0.6*$mass+130)")
    ]
]

cuts['ltau']=[
    [
        Template("(nlep==3) & (ntau<2) & (isoworst<0.35) & (sipworst<4) & ~lowMassRes & recoPass & (tWorstId>1)"),
        Template("(sumpt>0.85*$mass+125)"),
        Template("(dz>80)"),
        Template("(phi1<0.005*$mass+1.15) & (met>20)")
    ],[
        Template("(nlep>3) & (ntau<3) & (isoworst<0.35) & (sipworst<4) & ~lowMassRes & recoPass & (tWorstId>1)"),
        Template("(sumpt>$mass+100)"),
        Template("(dz>10)")
    ]
]

extras={'eeee':['(nele==3)','(nele==4)'], 'emem':['(nele>0) & (nmu>0)','(nele==2) & (nmu==2)'], 'mmmm':['(nmu==3)','(nmu==4)'], 'etet':['(nele>1)','(nele>1)'], 'mtmt':['(nmu>1)','(nmu>1)']}
