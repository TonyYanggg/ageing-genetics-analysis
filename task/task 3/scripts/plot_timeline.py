"""Plot an explicitly provisional timeline after the bibliography coverage audit."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/verified_timeline'


def plot(review=False):
    audit = json.loads((DATA / 'audit.json').read_text())
    if not audit.get('preview_ready') or audit['missing_pilot_pairs'] or audit['missing_reference_years'] or audit['pilot_identity_conflicts']:
        raise ValueError('Resolve missing pilot links, conflicting gene identities, or publication years before plotting')
    annual = pd.read_csv(DATA / 'annual_entry_timeline.csv')
    papers = pd.read_csv(DATA / 'annual_reference_timeline.csv')
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'svg.fonttype': 'none',
        'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.7, 'legend.frameon': False})
    out = ROOT / 'figures'
    if review:
        out = out / 'review'
    out.mkdir(parents=True, exist_ok=True)
    suffix = 'review' if review else 'exploratory'
    extensions = ['png', 'svg', 'pdf'] if review else ['png', 'svg']
    total = annual.groupby('year').cumulative_entries.sum()
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 5.3), sharex=True,
                            gridspec_kw={'height_ratios': [1.5, 1]}, layout='constrained')
    axes[0].step(total.index, total.values, where='post', color='#396A93', linewidth=1.8)
    axes[0].set_ylabel('Cumulative GenAge entries')
    axes[0].set_title('a   Earliest dated associated publication per entry', loc='left', weight='bold')
    axes[1].bar(papers.year, papers.unique_libage_references, width=0.85, color='#658B9A')
    axes[1].set_ylabel('Unique LibAge references')
    axes[1].set_xlabel('Publication year')
    axes[1].set_title('b   Publications linked to mapped GenAge entries', loc='left', weight='bold')
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    heading = 'GenAge reference timeline' if review else 'Exploratory timeline'
    fig.suptitle(f'{heading} | {audit["dated_entries"]:,}/{audit["total_entries"]:,} entries dated', fontsize=10)
    for ext in extensions:
        fig.savefig(out / f'timeline_{suffix}.{ext}', dpi=300 if review else 220)
    plt.close(fig)
    organisms = annual.organism.unique()
    coverage = pd.read_csv(DATA / 'entry_coverage.csv') if review else None
    fig, axes = plt.subplots(3, 3, figsize=(7.1, 6), sharex=True, layout='constrained')
    for ax, organism in zip(axes.flat, organisms):
        frame = annual[annual.organism.eq(organism)]
        ax.step(frame.year, frame.cumulative_entries, where='post', color='#396A93', linewidth=1.5)
        ax.set_title(organism, fontsize=8, fontstyle='italic')
        if review:
            count = len(coverage[coverage.organism.eq(organism)])
            ax.text(0.03, 0.9, f'{int(frame.cumulative_entries.iloc[-1])}/{count} dated',
                    transform=ax.transAxes, fontsize=7, va='top')
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=3))
    fig.supylabel('Cumulative GenAge entries (separate y-axis scales)')
    fig.supxlabel('Earliest dated associated publication year')
    fig.suptitle('GenAge reference timelines by organism' if review else 'Exploratory organism timelines', fontsize=10)
    for ext in extensions:
        fig.savefig(out / f'timeline_by_organism_{suffix}.{ext}', dpi=300 if review else 220)
    plt.close(fig)
    inputs = [DATA / name for name in ['annual_entry_timeline.csv', 'annual_reference_timeline.csv', 'audit.json']]
    if review:
        inputs.append(DATA / 'entry_coverage.csv')
    outputs = [out / (stem + '_' + suffix + '.' + ext) for stem in ['timeline', 'timeline_by_organism'] for ext in extensions]
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        script=str(Path(__file__).relative_to(ROOT)), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        files={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs + outputs},
        status='Exploratory drafts; see verified_timeline_report.md for exclusions and interpretation')
    (out / 'timeline_provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'{suffix} figures written; see coverage audit before interpretation.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--review', action='store_true', help='Write separate PNG/SVG/PDF review exports')
    plot(review=parser.parse_args().review)
