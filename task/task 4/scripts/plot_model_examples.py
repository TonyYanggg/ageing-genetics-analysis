"""Display selected explicit gene pairs; never fill missing edges transitively.

Figure contract: qualitative pair table illustrating annotation agreement,
discordance and sparse-organism contributions. Python-only, 183 x 153 mm,
editable SVG/PDF and 600-dpi TIFF. No effect-size or mechanistic inference.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collect_oma_pairwise import OUT, io

GROUPS=[('S6 kinase: annotation agreement',{'578','583','1737'},False),
        ('SOD: annotation discordance',{'984','621','622'},False),
        ('Sparse-organism contributions',{'2025','1792'},True)]
SHORT={'Caenorhabditis elegans':'C. elegans','Drosophila melanogaster':'D. melanogaster',
       'Mus musculus':'M. musculus','Saccharomyces cerevisiae':'S. cerevisiae',
       'Caenorhabditis briggsae':'C. briggsae','Danio rerio':'D. rerio'}
COLORS={'Anti-Longevity':'#DBE8EE','Pro-Longevity':'#F1E5CD',
        'discordant_or_mixed':'#E9DDEB','unresolved':'#EEEEEE'}


def selected_pairs(pairs):
    selected=[]
    for group,ids,any_endpoint in GROUPS:
        for row in pairs:
            left=set(row['genage_ids1'].split(';'))
            right=set(row['genage_ids2'].split(';'))
            keep=bool((left|right)&ids) if any_endpoint else bool(left&ids) and bool(right&ids)
            if keep:selected.append(dict(row,example_group=group))
    return selected


def main():
    io.OUT=OUT/'refined'
    rows=selected_pairs(io.read_csv(io.OUT/'longevity_pairs.csv'))
    io.write_csv('figure_example_pairs.csv',rows,list(rows[0]))
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
                         'font.size':7,'svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(183/25.4,153/25.4))
    ax=fig.add_axes([.04,.23,.92,.62]);ax.set_axis_off()
    labels={'Anti-Longevity':'Anti / Anti','Pro-Longevity':'Pro / Pro',
            'discordant_or_mixed':'Discordant','unresolved':'Unresolved'}
    table=[];group_rows=[];data_rows=[]
    for group,_,_ in GROUPS:
        group_rows.append(len(table)+1);table.append([group,'','',''])
        for row in rows:
            if row['example_group']!=group:continue
            table.append([SHORT[row['species1']]+'\n'+row['symbols1'],
                          SHORT[row['species2']]+'\n'+row['symbols2'],
                          row['pair_rel_types'],labels[row['direction_agreement']]])
            data_rows.append((len(table),row['direction_agreement']))
    t=ax.table(cellText=table,colLabels=['Model gene 1','Model gene 2','OMA relation','GenAge labels'],
               colWidths=[.32,.32,.14,.22],cellLoc='left',colLoc='left',bbox=[0,0,1,1])
    t.auto_set_font_size(False);t.set_fontsize(6.5)
    for (r,c),cell in t.get_celld().items():
        cell.set_edgecolor('#DDDDDD');cell.set_linewidth(.35)
        if r==0:cell.set_text_props(weight='bold');cell.set_facecolor('#E8EFF2')
        if r in group_rows:
            cell.set_facecolor('#F4F4F4');cell.set_text_props(weight='bold',fontsize=6.3)
    for r,category in data_rows:t[r,3].set_facecolor(COLORS[category])
    fig.text(.04,.95,'Model genes behind selected ortholog pairs',fontsize=11)
    fig.text(.04,.907,'Illustrative selection • explicit source edges only • not independent gene families',fontsize=7.5)
    fig.text(.04,.175,'S6-kinase experiments differ: worm inhibition, fly dominant-negative expression, mouse deletion.',fontsize=6.6)
    fig.text(.04,.14,'Mouse benefit is female-specific; linked erratum remains unchecked. SOD labels are catalogue-level only.',fontsize=6.6)
    fig.text(.04,.105,'No worm–mouse S6-kinase edge was recovered; the two displayed edges do not imply a third.',fontsize=6.6)
    fig.text(.04,.07,'All nine organisms and missing-source states are shown in Figure 1. Relation n:1 retains source multiplicity.',fontsize=6.6)
    folder=io.TASK/'figures/oma_refined';folder.mkdir(exist_ok=True,parents=True)
    for ext in ['png','svg','pdf','tiff']:
        kwargs={'dpi':600 if ext=='tiff' else 200}
        if ext=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(folder/('model_gene_examples.'+ext),**kwargs)
    plt.close(fig)
    print('Example pairs:',len(rows))


if __name__=='__main__':main()
