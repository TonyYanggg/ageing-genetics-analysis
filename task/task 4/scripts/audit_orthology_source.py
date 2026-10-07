"""Freeze a specific cross-resource disagreement; not a full orthology validation."""
import json
from collect_model_orthology import cached, OUT


def main():
    odb = cached('validation_coq2_search', 'https://data.orthodb.org/v12/search',
                 {'gid': '40988', 'level': '2759'})
    cached('validation_coq2_details', 'https://data.orthodb.org/v12/ogdetails',
           {'id': '7227_0:000d2a'})
    cached('validation_chaperone_group', 'https://data.orthodb.org/v12/group',
           {'id': '3459473at2759'})
    fly = cached('validation_oma_coq2', 'https://omabrowser.org/api/protein/DROME19540/')
    yeast = cached('validation_oma_hsp104', 'https://omabrowser.org/api/protein/P31539/')
    assert fly['roothog_id'] != yeast['roothog_id']
    result = dict(status='source_audit_only', accepted_biological_group_count=None,
                  concern='OrthoDB Eukaryota group combines Coq2 and HSP104; OMA assigns different root HOGs.',
                  orthodb_group='3459473at2759', orthodb_search=odb,
                  oma_coq2_root=fly['roothog_id'], oma_hsp104_root=yeast['roothog_id'],
                  limitation='One independent cross-check; neither resource is presumed ground truth. Full group validation remains pending.')
    (OUT/'source_validation.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
