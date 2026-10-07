"""Plot source-specific pair counts, retaining all nine inventory organisms."""
from collections import Counter
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collect_oma_pairwise import OUT, SPECIES, io


def main(with_bridges=False):
    dataset=OUT/'refined' if with_bridges else OUT
    io.OUT=dataset
    names=list(SPECIES)
    labels=['S. cerevisiae','C. elegans','D. melanogaster','M. musculus','S. pombe',
            'P. anserina','C. briggsae','D. rerio','M. auratus']
    pairs=io.read_csv(dataset/'longevity_pairs.csv')
    coverage={r['species']:r for r in io.read_csv(dataset/'coverage.csv')}
    counts=Counter(tuple(sorted((r['species1'],r['species2']))) for r in pairs)
    matrix=np.full((9,9),np.nan)
    data=[]
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            state='diagonal' if i==j else 'unassessed_source_gap' if not(SPECIES[a][0] and SPECIES[b][0]) else 'assessed'
            n=counts[tuple(sorted((a,b)))] if state=='assessed' else ''
            if state=='assessed':matrix[i,j]=n
            data.append(dict(species1=a,species2=b,status=state,pairs=n))
    io.write_csv('figure_pair_counts.csv',data,list(data[0]))
    assert int(np.nansum(matrix))==2*len(pairs)
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':7,'svg.fonttype':'none','pdf.fonttype':42,
        'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(183/25.4,145/25.4))
    a=fig.add_axes([.17,.38,.27,.44])
    b=fig.add_axes([.63,.38,.35,.44])
    for i,s in enumerate(names):
        r=coverage[s]
        if SPECIES[s][0]:
            pct=100*float(r['mapped_fraction'])
            a.barh(i,pct,color='#759BA9',height=.6)
            a.text(pct+2,i,f"{r['mapped_genes']}/{r['identified_genes']}",va='center',fontsize=6)
        else:a.text(2,i,'NA — source gap',va='center',fontsize=6,color='#666666')
    a.set_yticks(range(9),labels,fontsize=6)
    a.set_ylim(8.5,-.5);a.set_xlim(0,150);a.set_xticks([0,50,100])
    a.set_xlabel('Mapped identifiable genes (%)',fontsize=6)
    a.set_title('a  Exact-ID coverage',loc='left',fontweight='bold',fontsize=8,pad=10)
    cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#EEEEEE')
    b.imshow(matrix,cmap=cmap,vmin=0,vmax=np.nanmax(matrix),aspect='auto')
    for i in range(9):
        for j in range(9):
            value=matrix[i,j]
            text='—' if i==j else 'NA' if np.isnan(value) else str(int(value))
            b.text(j,i,text,ha='center',va='center',fontsize=6,
                   color='white' if value>np.nanmax(matrix)*.6 else '#222222')
    b.set_xticks(range(9),labels,rotation=55,ha='right',fontsize=6)
    b.set_yticks(range(9),labels,fontsize=6)
    b.tick_params(length=0)
    b.set_xticks(np.arange(-.5,9,1),minor=True);b.set_yticks(np.arange(-.5,9,1),minor=True)
    b.grid(which='minor',color='white',linewidth=.6);b.tick_params(which='minor',length=0)
    for spine in b.spines.values():spine.set_visible(False)
    b.set_title('b  Lifespan-evidence pairs',loc='left',fontweight='bold',fontsize=8,pad=10)
    fig.text(.05,.95,'Explicit ortholog pairs with longevity evidence',fontsize=11)
    fig.text(.05,.905,f'GenAge Build 21 • OMA All.May2026 • {len(pairs)} unique gene pairs',fontsize=8)
    fig.text(.05,.17,'All nine organisms retained; OMA represents seven. NA cells are unassessed, not biological zeros.',fontsize=7)
    fig.text(.05,.13,'Symmetric matrix: each pair appears twice. Pairs share genes and are not independent orthogroups.',fontsize=7)
    fig.text(.05,.09,'Zero = no qualifying pair recovered. Ensembl supplementary checks are not included in this OMA figure.',fontsize=7)
    folder=io.TASK/('figures/oma_refined' if with_bridges else 'figures/oma_pairwise');folder.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg','pdf','tiff']:
        kwargs={'dpi':600 if ext=='tiff' else 200}
        if ext=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(folder/('pairwise_evidence.'+ext),**kwargs)
    plt.close(fig)
    print('Figure source counts reconciled:',len(pairs))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--with-bridges',action='store_true')
    main(parser.parse_args().with_bridges)
