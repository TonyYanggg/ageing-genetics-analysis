"""Independent review-package consistency, provenance and export checks.

Does not import analysis or plotting code. Complements raw-cache reconstruction.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET
from PIL import Image
import fitz

TASK=Path(__file__).resolve().parents[1]
BASE=TASK/'results/oma_pairwise'
OUT=BASE/'refined'


def read(path,delimiter=','):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream,delimiter=delimiter))


def key(row):
    return tuple(sorted(((row['species1'],row['entrez1']),(row['species2'],row['entrez2']))))


def main():
    raw=read(TASK.parents[1]/'GenAge Project/models_genes/genage_models.csv')
    raw_by_id={r['GenAge ID']:r for r in raw}
    ledger=read(OUT/'resolution_ledger.csv')
    assert len(ledger)==len(raw)==len({r['genage_id'] for r in ledger})==2205
    assert set(raw_by_id)=={r['genage_id'] for r in ledger}
    identities={r['genage_id']:r for r in read(BASE/'gene_identity_audit.csv')}
    for row in ledger:
        original=raw_by_id[row['genage_id']]
        assert (row['species'],row['entrez_id'],row['symbol'])==(original['organism'],original['entrez gene id'],original['symbol'])
        if row['status']=='mapped':
            assert row['resolution_category']=='mapped_model_protein'
            assert identities[row['genage_id']]['identity_status']=='matched_taxon'
        else:assert row['resolution_category']!='mapped_model_protein'
    gaps=read(OUT/'remaining_gaps.csv')
    assert gaps==[r for r in ledger if r['status']!='mapped']
    resolution=json.loads((OUT/'resolution_summary.json').read_text())
    assert resolution['resolution_counts']==dict(Counter(r['resolution_category'] for r in ledger))
    assert resolution['mapped_records']+resolution['remaining_records']==len(raw)

    pairs=read(OUT/'longevity_pairs.csv')
    summary=json.loads((OUT/'analysis_summary.json').read_text())
    for pair in pairs:
        evidence=[[raw_by_id[x] for x in pair['genage_ids'+str(side)].split(';')
                   if raw_by_id[x]['lifespan effect'] in {'Increase','Decrease','Increase and Decrease'}]
                  for side in (1,2)]
        assert all(evidence)
        signs=[{r['longevity influence'] for r in group} for group in evidence]
        known={'Anti-Longevity','Pro-Longevity'}
        union=set.union(*signs)
        expected=('discordant_or_mixed' if len(union&known)>1 else
                  next(iter(union)) if len(union)==1 and union<=known else 'unresolved')
        assert expected==pair['direction_agreement']
        assert pair['both_directional']==str(all(bool(s&known) for s in signs))
    assert summary['agreement']==dict(Counter(r['direction_agreement'] for r in pairs))
    baseline={key(r) for r in read(BASE/'longevity_pairs.csv')}
    assert baseline<={key(r) for r in pairs}

    examples=read(OUT/'figure_example_pairs.csv')
    lookup={key(r):r for r in pairs}
    expected_examples=set()
    for row in pairs:
        a=set(row['genage_ids1'].split(';'));b=set(row['genage_ids2'].split(';'))
        if (a|b)&{'2025','1792'} or any(a&ids and b&ids for ids in ({'578','583','1737'},{'984','621','622'})):
            expected_examples.add(key(row))
    assert len(examples)==len({key(r) for r in examples})
    assert {key(r) for r in examples}==expected_examples
    for row in examples:
        for field,value in lookup[key(row)].items():assert row[field]==value

    review=read(TASK/'notes/experiment_review.tsv','\t')
    references=read(BASE/'experiment_review_sources.csv')
    primary={(r['genage_id'],r['pmid']) for r in references}
    assert len(review)==len({r['genage_id'] for r in review})
    for row in review:
        assert row['genage_id'] in raw_by_id
        assert all((row['genage_id'],pmid) in primary for pmid in row['pmid'].split(';'))
        assert all(row[x] for x in ['review_level','intervention','context','lifespan_finding','identity_decision','claim_limit'])
    reference_file=TASK.parents[1]/'task/task 3/results/verified_timeline/entry_references.csv'
    assert {r['reference_table_sha256'] for r in references}=={hashlib.sha256(reference_file.read_bytes()).hexdigest()}

    exports=[]
    for name,height in [('pairwise_evidence',145),('model_gene_examples',153)]:
        folder=TASK/'figures/oma_refined'
        svg=ET.parse(folder/(name+'.svg')).getroot()
        text=[''.join(e.itertext()) for e in svg.iter() if e.tag.endswith('}text')]
        assert text
        if name=='model_gene_examples':
            for r in examples:
                assert any(r['symbols1'] in s for s in text)
                assert any(r['symbols2'] in s for s in text)
        else:assert any(str(len(pairs))+' unique gene pairs' in s for s in text)
        with Image.open(folder/(name+'.tiff')) as im:
            assert all(abs(float(d)-600)<1 for d in im.info['dpi'])
            assert abs(im.width/600*25.4-183)<.1
            assert abs(im.height/600*25.4-height)<.1
        with Image.open(folder/(name+'.png')) as im:im.verify()
        pdf=(folder/(name+'.pdf')).read_bytes()
        assert pdf.startswith(b'%PDF-') and b'/FontFile2' in pdf
        with fitz.open(folder/(name+'.pdf')) as document:
            assert len(document)==1
            page=document[0]
            assert abs(page.rect.width/72*25.4-183)<.1
            assert abs(page.rect.height/72*25.4-height)<.1
            assert page.get_text().strip()
            for word in page.get_text('words'):
                assert page.rect.contains(fitz.Rect(word[:4])), ('Clipped text',word)
            qa=folder/'qa';qa.mkdir(exist_ok=True)
            page.get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(qa/(name+'_pdf.png'))
        exports.append(dict(name=name,editable_svg=True,pdf_embedded_truetype=True,tiff_dpi=600))
    report=dict(passed=True,records=len(ledger),remaining_records=len(gaps),
                classified_pair_directions=len(pairs),baseline_pairs_retained=len(baseline),
                added_pairs=len(pairs)-len(baseline),verified_example_pairs=len(examples),
                selected_review_records=len(review),exports=exports,
                limitation='Consistency and file-format checks; not exhaustive experimental validation or proof of conservation.')
    (OUT/'review_package_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':main()
