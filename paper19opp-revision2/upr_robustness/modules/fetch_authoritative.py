#!/usr/bin/env python3
"""Preserve authoritative human GO annotations and UniProt symbol mappings."""
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

OUT = Path(__file__).resolve().parent
RAW = OUT / 'authoritative_raw'
RAW.mkdir(parents=True, exist_ok=True)
records = []

def get(url, filename, accept='application/json'):
    path = RAW / filename
    if path.exists():
        body = path.read_bytes()
        headers = {'cache_reuse': 'existing exact response'}
    else:
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'Accept': accept})
                with urllib.request.urlopen(request, timeout=60) as response:
                    body = response.read()
                    headers = dict(response.headers)
                path.write_bytes(body)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 + attempt)
    records.append(dict(file=str(path.relative_to(OUT)), url=url,
                        accessed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                        response_headers=headers))
    return body

terms = ['GO:0034975', 'GO:0036503', 'GO:0045047', 'GO:0140467', 'GO:0030968']
term_url = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/' + ','.join(terms)
get(term_url, 'QuickGO_terms.json')
for endpoint, filename in [('ontology/go/about', 'QuickGO_ontology_about.json'),
                           ('annotation/about', 'QuickGO_annotation_about.json')]:
    try:
        get('https://www.ebi.ac.uk/QuickGO/services/' + endpoint, filename)
    except Exception as error:
        records.append(dict(url=endpoint, error=str(error)))

annotations = []
for term in terms:
    page = 1
    while True:
        query = urllib.parse.urlencode(dict(taxonId=9606, goId=term,
                    goUsage='descendants', goUsageRelationships='is_a,part_of',
                    limit=200, page=page))
        body = get('https://www.ebi.ac.uk/QuickGO/services/annotation/search?' + query,
                   'QuickGO_' + term.replace(':', '') + '_page' + str(page) + '.json')
        result = json.loads(body)
        for row in result['results']:
            row = dict(row)
            row['requested_parent_go_id'] = term
            annotations.append(row)
        print(term, 'page', page, 'hits', result['numberOfHits'], flush=True)
        if page >= result['pageInfo']['total']:
            break
        page += 1

get_ids = set()
for row in annotations:
    if row['geneProductId'].startswith('UniProtKB:'):
        get_ids.add(row['geneProductId'].split(':', 1)[1].split('-')[0])
reactome = json.loads((RAW / 'Reactome_R-HSA-381119_participants.json').read_text())
for participant in reactome:
    for entity in participant.get('refEntities', []):
        if entity.get('schemaClass') == 'ReferenceGeneProduct' and entity.get('stId', '').startswith('uniprot:'):
            get_ids.add(entity['identifier'].split('-')[0])

mapping = []
ids = sorted(get_ids)
for start in range(0, 0, 60):
    batch = ids[start:start+60]
    query = ' OR '.join('accession:' + value for value in batch)
    params = urllib.parse.urlencode(dict(query='(' + query + ')',
                    fields='accession,gene_primary,xref_geneid', format='tsv'))
    body = get('https://rest.uniprot.org/uniprotkb/stream?' + params,
               'UniProt_primary_mapping_' + str(start//60+1) + '.tsv', 'text/tab-separated-values')
    mapping.extend(csv.DictReader(io.StringIO(body.decode()), delimiter='\t'))
    print('UniProt mapped batch', start//60+1, flush=True)

(OUT / 'authoritative_annotations.json').write_text(json.dumps(annotations, indent=2))
(OUT / 'authoritative_uniprot_mapping.json').write_text(json.dumps(mapping, indent=2))
(OUT / 'download_manifest.json').write_text(json.dumps(records, indent=2))
print('Annotations', len(annotations), 'UniProt identifiers', len(ids), 'mapping rows', len(mapping), flush=True)
print('Membership uses official QuickGO protein-product symbols; accession-only symbols remain unmapped. UniProt bulk mapping was unavailable and is not imputed.', flush=True)
