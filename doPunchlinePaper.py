from extLibFull import texLabel,mpl
mpl.rcParams.update({'legend.fontsize':16})
mpl.rcParams.update({'text.fontsize':24})
from input import lumi
from pylab import *

tevLim=(225, 210, 150, 112, 144, 128, 0, 0, 0, 0)

atlasLim=(0,0,355,0,0,0,0,0,0,0)

cms36pb=( 144,154,156,106,106,0,116,131,130,127)

cms1fb=(313, 312, 312, 253, 265, 0, 269, 296, 291, 289)

cms4fb=(445,455,457,352,369,195,380,410,406,399)
cms5fb=(444,453,459,373,375,204,383,408,403,400)

yax=[8.5-i for i in range(10)]
texLab=[
    r"$\mathbf{e^{\pm}e^{\pm}=100\%}$  ",
    r"$\mathbf{e^{\pm}\mu^{\pm}=100\%}$  ",
    r"$\mathbf{\mu^{\pm}\mu^{\pm}=100\%}$  ",
    r"$\mathbf{e^{\pm}\tau^{\pm}=100\%}$  ",
    r"$\mathbf{\mu^{\pm}\tau^{\pm}=100\%}$  ",
    r"$\mathbf{\tau^{\pm}\tau^{\pm}=100\%}$  ",
    r"BP1: normal hierarchy  ",
    r"BP2: inverse hierachy  ",
    r"BP3: degenerate masses  ",
    r"BP4: equal branchings  "
]
fig=figure(figsize=[16,8],dpi=75,facecolor='w')
ax=fig.add_axes([0.22,0.1,0.77,0.88],xlabel=r'95\% CL limit on mass of $\Phi^{\pm\pm}$ [GeV]')
#barh(yax,cms4fb,color='r',align='center',label=r'CMS $\int\mathcal{L} \textrm{dt}=4.6$ fb$^{-1}$')
barh(yax,cms5fb,color='r',align='center',label=r'CMS $\int\mathcal{L} \textrm{dt}=4.9$ fb$^{-1}$')
barh(yax,atlasLim,color='0.75',align='center',label=r'ATLAS $\int\mathcal{L} \textrm{dt}=1.6$ fb$^{-1}$')
barh(yax,tevLim,color='b',align='center',label='Tevatron')
#barh(yax,cms5fb,color='r',align='center',label='CMS $\int\mathcal{L}\text{dt}=4.91$ fb$^{-1}$')
#scatter(tevLim,yax,c='b',s=150,marker='d',label='Tevatron',zorder=3)
#scatter(atlasLim,yax,c='0.75',s=150,marker='^',label=r'ATLAS $\int\mathcal{L} \textrm{dt}=1.6$ fb$^{-1}$',zorder=3)
gca().set_yticks(yax)
gca().set_yticklabels(texLab)
gca().set_ylim(-1,9)
gca().set_xlim(90,600)
handles, labels = ax.get_legend_handles_labels()
#handles[0],handles[2] = handles[2], handles[0]
#labels[0],labels[2] = labels[2], labels[0]

#handles.insert(-3,handles[-1])
#labels.insert(-3,labels[-1])
#del handles[5],labels[5]

legend(handles,labels,loc=4,scatterpoints=1)
figtext(0.81,0.93,r"CMS $\sqrt s=7$ TeV",fontsize=25)
savefig('punchline.png')
savefig('punchline.pdf')
draw()
