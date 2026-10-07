"""Propose exact stable-ID bridges for native protein-coding mapping gaps.

This produces an audited proposal table; it does not silently replace baseline mappings.
"""
from collections import defaultdict
import json
from collect_oma_pairwise import OUT, SPECIES, io


def main():
    identity={r['genage_id']:r for r in io.read_csv(OUT/'gene_identity_audit.csv')}
    baseline=io.read_csv(OUT/'record_mapping_audit.csv')
    metadata={}
    for path in io.CACHE.glob('ncbi_unmapped_dataset_*.json'):
        if path.name.endswith('.meta.json'):continue
        for report in json.loads(path.read_text())['reports']:
            gene=report.get('gene',{})
            if gene.get('gene_id'):metadata[gene['gene_id']]=gene
    candidates=[]
    for r in baseline:
        audit=identity[r['genage_id']]
        if r['status']!='no_exact_xref' or audit['identity_status']!='matched_taxon' or audit['gene_type']!='PROTEIN_CODING':continue
        gene=metadata[r['entrez_id']]
        refs=[('UniProt',x) for x in gene.get('swiss_prot_accessions',[])]+[('Ensembl',x) for x in gene.get('ensembl_gene_ids',[])]
        authority=gene.get('nomenclature_authority',{})
        if r['species']=='Saccharomyces cerevisiae' and gene.get('locus_tag'):
            refs.append(('NCBI_locus_tag',gene['locus_tag']))
        if authority.get('authority') in {'FLYBASE','WormBase','SGD','MGI','PomBase'} and authority.get('identifier'):
            refs.append((authority['authority'],authority['identifier']))
        for kind,ref in refs:
            candidates.append(dict(genage_id=r['genage_id'],species=r['species'],entrez_id=r['entrez_id'],
                                   symbol=r['symbol'],identifier_type=kind,identifier=ref))
    io.write_csv('identifier_bridge_candidates.csv',candidates,list(candidates[0]))
    identifiers=sorted({r['identifier'] for r in candidates})
    matches=defaultdict(set)
    xref_audit=[]
    for number,identifier in enumerate(identifiers):
        data=io.cached('oma_xref_bridge_'+identifier.replace(':','_'),'https://omabrowser.org/api/xref/',
                       {'search':identifier})
        for row in data:
            # Search can return substrings; accept literal stable-ID equality only.
            if row['xref']!=identifier:continue
            matches[(identifier,str(row['genome']['taxon_id']))].add(row['omaid'])
            xref_audit.append(dict(identifier=identifier,oma_protein=row['omaid'],taxon=row['genome']['taxon_id'],
                                  source=row['source'],sequence_match=row['seq_match']))
        print('Bridge identifiers',number+1,'/',len(identifiers),flush=True)
    io.write_csv('identifier_bridge_xref_audit.csv',xref_audit,['identifier','oma_protein','taxon','source','sequence_match'])
    output=[]
    for c in candidates:
        ps=matches[(c['identifier'],str(SPECIES[c['species']][1]))]
        output.append(dict(c,oma_proteins=';'.join(sorted(ps)),status='exact_bridge_found' if ps else 'no_bridge_hit'))
    io.write_csv('identifier_bridge_results.csv',output,list(output[0]))
    print('Genes with candidate recovery',len({(r['species'],r['entrez_id']) for r in output if r['oma_proteins']}))


if __name__=='__main__':main()
