"""Classify remaining gaps without rewriting raw annotations or inventing IDs."""
from collections import Counter
import csv
import json
from collect_oma_pairwise import OUT, io


def main():
    records=io.read_csv(OUT/'refined/record_mapping_audit.csv')
    identities={r['genage_id']:r for r in io.read_csv(OUT/'gene_identity_audit.csv')}
    with (io.TASK/'notes/experiment_review.tsv').open(encoding='utf-8',newline='') as f:
        review={r['genage_id']:r for r in csv.DictReader(f,delimiter='\t')}
    rows=[]
    for r in records:
        identity=identities[r['genage_id']]
        decision=review.get(r['genage_id'],{}).get('identity_decision','')
        if r['status']=='mapped':resolution='mapped_model_protein'
        elif decision in {'catalogue_gene_paper_conflict','catalogue_taxon_conflict'}:resolution=decision
        elif decision in {'confirmed_heterologous_construct','construct_species_unresolved'}:resolution=decision
        elif identity['identity_status']=='taxon_mismatch':resolution='taxon_mismatch_unresolved'
        elif identity['identity_status']=='discontinued_or_replaced':resolution='discontinued_identifier'
        elif identity['gene_type'] in {'ncRNA','miscRNA'}:resolution='noncoding_not_assessed_by_protein_orthology'
        elif identity['gene_type']=='PSEUDO':resolution='current_pseudogene_annotation'
        elif r['status']=='ambiguous_identifier':resolution='ambiguous_gene_to_protein_identity'
        elif r['status']=='missing_id':resolution='missing_stable_identifier'
        elif r['status']=='source_species_unavailable':resolution='OMA_source_species_gap'
        else:resolution='protein_identifier_gap_after_stable_bridge_search'
        rows.append(dict(r,ncbi_identity_status=identity['identity_status'],current_ncbi_type=identity['gene_type'],
            current_ncbi_symbol=identity['ncbi_symbol'],resolution_category=resolution,
            experiment_review_decision=decision,review_pmid=review.get(r['genage_id'],{}).get('pmid','')))
    io.OUT=OUT/'refined'
    io.write_csv('resolution_ledger.csv',rows,list(rows[0]))
    gaps=[r for r in rows if r['status']!='mapped']
    io.write_csv('remaining_gaps.csv',gaps,list(rows[0]))
    summary=dict(records=len(rows),mapped_records=sum(r['status']=='mapped' for r in rows),
        remaining_records=len(gaps),resolution_counts=dict(Counter(r['resolution_category'] for r in rows)),
        limitation='Current gene-type metadata is queried for baseline unmapped records only. Noncoding loci remain in the inventory but are not assessed by OMA protein orthology. No source record has been corrected in place.')
    (OUT/'refined/resolution_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))


if __name__=='__main__':main()
