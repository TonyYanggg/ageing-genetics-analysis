"""Render the nine-organism source-group evidence atlas with explicit missingness."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from collect_model_orthology import OUT, TASK, SPECIES, read_csv, write_csv

SHORT = ['S. cerevisiae','C. elegans','D. melanogaster','M. musculus','S. pombe',
         'P. anserina','C. briggsae','D. rerio','M. auratus']
STATE_ORDER = ['No recovered member','Ortholog only','Catalogue only','Other','Pro','Anti','Mixed']
COLORS = ['#F2F2F2','#D7E3E8','#E4DECD','#D5D5D5','#A5CDD8','#E7BF8B','#BDB0D4']


def main():
    figures=TASK/'figures/model_orthology'
    figures.mkdir(parents=True,exist_ok=True)
    groups=read_csv(OUT/'recurrent_groups_eukaryota.csv')
    coverage=read_csv(OUT/'coverage.csv')
    members=read_csv(OUT/'recurrent_group_members.csv')
    present=defaultdict(set)
    for row in members:
        present[(row['og'].rsplit('/',1)[-1],row['organism'].rsplit('/',1)[-1])].add(row['gene'])
    cells=read_csv(OUT/'group_species_cells.csv')
    lookup={}
    rows=[]
    group_ids={r['group_id'] for r in groups}
    for cell in cells:
        if cell['group_id'] not in group_ids:
            continue
        key=(cell['group_id'],SPECIES[cell['species']])
        cell['source_member_count']=len(present[key])
        if cell['state']=='No matched GenAge evidence':
            cell['state']='Ortholog only' if present[key] else 'No recovered member'
        rows.append(cell)
        lookup[(cell['group_id'],cell['species'])]=cell
    write_csv('recurrent_atlas_cells.csv',rows,list(rows[0]))
    displayed=groups[:12]
    write_csv('figure_displayed_groups.csv',displayed,list(groups[0]))
    agreements=['Pro','Anti','Discordant','Unresolved']
    distribution=Counter((int(r['lifespan_species']),r['direction_agreement']) for r in groups)
    dist_rows=[dict(species=n,agreement=a,groups=distribution[(n,a)]) for n in range(2,10) for a in agreements]
    write_csv('figure_recurrence_distribution.csv',dist_rows,list(dist_rows[0]))
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,
        'axes.spines.right':False,'axes.linewidth':.6,'legend.frameon':False})
    fig=plt.figure(figsize=(183/25.4,220/25.4))
    grid=fig.add_gridspec(2,2,height_ratios=[1,1.5],left=.18,right=.985,bottom=.19,top=.90,
                        wspace=.48,hspace=.68)
    a,c,b=fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])
    cov={r['species']:r for r in coverage}
    vals=[100*float(cov[s]['mapped_fraction']) for s in SPECIES]
    a.barh(range(9),vals,color='#658F9F',height=.62)
    a.set_yticks(range(9),SHORT,fontsize=6.5)
    a.invert_yaxis()
    a.set_xlim(0,140)
    a.set_xticks([0,50,100])
    a.set_xlabel('Genes mapped to common-level groups (%)',fontsize=6.5)
    for i,s in enumerate(SPECIES):
        r=cov[s]
        a.text(vals[i]+2,i,f"{r['mapped_genes']}/{r['unique_identified_genes']}",va='center',fontsize=6)
    a.set_title('a   Mapping coverage',loc='left',fontweight='bold',pad=9)
    bottoms=np.zeros(8)
    for agreement,color in zip(agreements,[COLORS[4],COLORS[5],COLORS[6],COLORS[3]]):
        counts=np.array([distribution[(n,agreement)] for n in range(2,10)])
        c.bar(range(2,10),counts,bottom=bottoms,color=color,width=.7,label=agreement)
        bottoms+=counts
    for n,total in zip(range(2,10),bottoms):
        c.text(n,total,str(int(total)),ha='center',va='bottom',fontsize=6)
    c.set_ylim(0,max(bottoms)*1.22)
    c.set_xticks(range(2,10))
    c.set_xlabel('Species with lifespan evidence',fontsize=6.5)
    c.set_ylabel('Recurrent ancestral groups')
    c.set_title('c   Recurrence and direction',loc='left',fontweight='bold',pad=9)
    c.legend(fontsize=6,ncol=2,loc='upper right')
    matrix=[[STATE_ORDER.index(lookup[(g['group_id'],s)]['state']) for s in SPECIES] for g in displayed]
    b.imshow(matrix,cmap=ListedColormap(COLORS),vmin=0,vmax=6,aspect='auto',interpolation='nearest')
    b.set_yticks(range(len(displayed)),[g['group_id'].replace('at2759','') for g in displayed],fontsize=6.5)
    b.set_xticks(range(9),SHORT,rotation=45,ha='left',rotation_mode='anchor',fontsize=6.2)
    b.xaxis.tick_top()
    b.tick_params(length=0,pad=4)
    code={'No recovered member':'?','Ortholog only':'O','Catalogue only':'C','Other':'L',
          'Pro':'P','Anti':'A','Mixed':'M'}
    truncated=0
    for i,g in enumerate(displayed):
        for j,s in enumerate(SPECIES):
            cell=lookup[(g['group_id'],s)]
            symbols=cell['model_symbols'].split(';') if cell['model_symbols'] else []
            label=code[cell['state']]
            if symbols:
                symbol=symbols[0]
                if len(symbol)>11:
                    symbol=symbol[:9]+'…'
                    truncated+=1
                label+='\n'+symbol
                if len(symbols)>1:
                    label+=f' +{len(symbols)-1}'
            b.text(j,i,label,ha='center',va='center',fontsize=5.4)
    b.set_xticks(np.arange(-.5,9,1),minor=True)
    b.set_yticks(np.arange(-.5,len(displayed),1),minor=True)
    b.grid(which='minor',color='white',linewidth=1)
    b.tick_params(which='minor',length=0)
    for spine in b.spines.values():
        spine.set_visible(False)
    b.set_ylabel('OrthoDB group ID (suffix: at2759)',fontsize=6.5)
    b.set_title('b   Selected recurrent groups: model genes and influence',loc='left',fontweight='bold',pad=61)
    fig.suptitle('Lifespan evidence across nine GenAge model organisms',x=.05,ha='left',fontsize=10,y=.977)
    fig.text(.05,.948,'SOURCE AUDIT ONLY • OrthoDB 12.2 • orthology validation pending',fontsize=7)
    fig.text(.06,.135,'P  Pro-longevity     A  Anti-longevity     M  Both directions     L  Other lifespan evidence',fontsize=7)
    fig.text(.06,.112,'O  Source member, no linked lifespan record     C  Catalogue record without lifespan effect',fontsize=6.6)
    fig.text(.06,.089,'?  No member recovered; not evidence of biological absence.   +n  Additional model genes in the cell.',fontsize=6.6)
    fig.text(.06,.055,'Group 3459473 has a cross-source conflict (Coq2/HSP104). Counts are NOT validated longevity orthogroups.',fontsize=6.5)
    fig.text(.06,.033,f'Shown: {len(displayed)} of {len(groups)} recurrent groups, ranked by species count, then fewer matched genes.',fontsize=6.5)
    for ext in ['png','svg','pdf','tiff']:
        options={'dpi':600 if ext=='tiff' else 200}
        if ext=='tiff':options['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(figures/('model_orthology_atlas.'+ext),**options)
    plt.close(fig)
    print(json.dumps(dict(displayed=len(displayed),total=len(groups),truncated_symbol_labels=truncated)))


if __name__=='__main__':
    main()
