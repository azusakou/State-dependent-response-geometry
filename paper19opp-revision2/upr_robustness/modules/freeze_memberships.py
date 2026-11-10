#!/usr/bin/env python3
"""Outcome-independent GO/Reactome memberships, HGNC mapping and fixed variants."""
import csv
import collections
import datetime
import hashlib
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
RAW = OUT / 'authoritative_raw'

def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))

def write_tsv(path, rows, fields):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

hgnc = [r for r in read_tsv(RAW/'HGNC_complete_set.tsv') if r['status']=='Approved']
approved = {r['symbol']:r for r in hgnc}
aliases, proteins = collections.defaultdict(set), collections.defaultdict(set)
for row in hgnc:
    for symbol in (row['alias_symbol']+'|'+row['prev_symbol']).split('|'):
        if symbol:
            aliases[symbol].add(row['symbol'])
    for protein in row['uniprot_ids'].split('|'):
        if protein:
            proteins[protein].add(row['symbol'])

def resolve(symbol, protein=''):
    if protein and len(proteins[protein.split('-')[0]])==1:
        return next(iter(proteins[protein.split('-')[0]])), 'HGNC_UniProt'
    if symbol in approved:
        return symbol, 'HGNC_approved_symbol'
    if len(aliases[symbol])==1:
        return next(iter(aliases[symbol])), 'HGNC_unique_previous_or_alias_symbol'
    return '', 'unmapped_or_ambiguous'

universe = [r['gene'] for r in read_tsv(OUT/'primary_gene_universe.tsv')]
assay = []
canonical_to_assay = collections.defaultdict(list)
for symbol in universe:
    canonical, method = resolve(symbol)
    assay.append(dict(assay_gene=symbol, canonical_gene=canonical, mapping_method=method))
    if canonical:
        canonical_to_assay[canonical].append(symbol)
assert all(len(v)==1 for v in canonical_to_assay.values()), 'Alias collision requires count-level aggregation, not duplicate score weighting'
write_tsv(OUT/'assay_HGNC_mapping.tsv', assay, list(assay[0]))

terms = {'GO:0034975':'ER_folding', 'GO:0036503':'ERAD',
         'GO:0045047':'ER_targeting', 'GO:0140467':'ISR',
         'GO:0030968':'GO_ER_UPR'}
term_info = {r['id']:r for r in json.load(open(RAW/'QuickGO_terms.json'))['results']}
sets = collections.defaultdict(set)
annotation_map=[]
for row in json.load(open(OUT/'authoritative_annotations.json')):
    include = row['geneProductId'].startswith('UniProtKB:') and 'NOT' not in row['qualifier'].split('|')
    protein = row['geneProductId'].split(':',1)[-1] if include else ''
    canonical, method = resolve(row.get('symbol',''), protein) if include else ('','excluded_nonprotein_or_NOT')
    if include and canonical:
        sets[terms[row['requested_parent_go_id']]].add(canonical)
    annotation_map.append(dict(module=terms[row['requested_parent_go_id']],
            parent_id=row['requested_parent_go_id'], annotation_go_id=row['goId'],
            gene_product_id=row['geneProductId'], annotated_symbol=row.get('symbol',''),
            canonical_gene=canonical, mapping_method=method, included=bool(include and canonical),
            qualifier=row['qualifier'], evidence=row['goEvidence'], annotation_date=row['date'],
            reference=row['reference']))
write_tsv(OUT/'GO_annotation_mapping.tsv',annotation_map,list(annotation_map[0]))

reactome_map=[]
for participant in json.load(open(RAW/'Reactome_R-HSA-381119_participants.json')):
    for entity in participant.get('refEntities',[]):
        if entity.get('schemaClass')!='ReferenceGeneProduct' or not entity.get('stId','').startswith('uniprot:'):
            continue
        symbol=entity['displayName'].split(' ',1)[1]
        canonical,method=resolve(symbol,entity['identifier'])
        if canonical:
            sets['Reactome_UPR'].add(canonical)
        reactome_map.append(dict(physical_entity_id=participant['peDbId'],reference_entity_id=entity['dbId'],
                protein_id=entity['identifier'],display_name=entity['displayName'],canonical_gene=canonical,
                mapping_method=method,included=bool(canonical)))
write_tsv(OUT/'Reactome_protein_reference_mapping.tsv',reactome_map,list(reactome_map[0]))

# Original Hallmark measurement definition is preserved, with its authoritative snapshot copied to a compact list.
import gzip
with gzip.open('/Users/wang/Documents/paper6yhn/resource/networks/hallmark_human.tsv.gz','rt') as handle:
    hallmark_rows=list(csv.DictReader(handle,delimiter='\t'))
hallmark_raw=sorted({r['target'] for r in hallmark_rows if r['source']=='UNFOLDED_PROTEIN_RESPONSE'})
assert len(hallmark_raw)==113
hallmark_canonical={resolve(g)[0] for g in hallmark_raw if resolve(g)[0]}
write_tsv(OUT/'original_Hallmark_UPR_full.tsv',[dict(gene=g,canonical_gene=resolve(g)[0])for g in hallmark_raw],['gene','canonical_gene'])

functional=['ER_folding','ERAD','ER_targeting','ISR']
definitions=[]
for module in functional:
    definitions.extend([(module+'_full',module,'full',sets[module]),
        (module+'_no_Hallmark',module,'Hallmark_disjoint',sets[module]-hallmark_canonical),
        (module+'_unique',module,'between_functional_module_unique',sets[module]-set().union(*(sets[k]for k in functional if k!=module)))])
for module in ['Reactome_UPR','GO_ER_UPR']:
    definitions.extend([(module+'_full',module,'full',sets[module]),
        (module+'_no_Hallmark',module,'Hallmark_disjoint',sets[module]-hallmark_canonical)])

members=[]
for axis,module,variant,genes in definitions:
    for g in sorted(genes):
        measured=canonical_to_assay.get(g,[])
        members.append(dict(axis=axis,parent_module=module,variant=variant,canonical_gene=g,
                            assay_gene=measured[0]if measured else'',measured=bool(measured)))
write_tsv(OUT/'locked_memberships.tsv',members,['axis','parent_module','variant','canonical_gene','assay_gene','measured'])

definitions_rows=[]
for axis,module,variant,genes in definitions:
    source_id=next((k for k,v in terms.items()if v==module),'R-HSA-381119.4')
    label=term_info[source_id]['name']if source_id in term_info else'Unfolded Protein Response (UPR)'
    definitions_rows.append(dict(axis=axis,parent_module=module,variant=variant,source_id=source_id,
       official_label=label,full_canonical_members=len(genes),measured_members=sum(g in canonical_to_assay for g in genes),
       source_url=('https://amigo.geneontology.org/amigo/term/'+source_id)if source_id.startswith('GO:')else'https://reactome.org/content/detail/R-HSA-381119'))
write_tsv(OUT/'locked_axis_definitions.tsv',definitions_rows,list(definitions_rows[0]))
base_rows=[dict(module=k,canonical_gene=g,assay_gene=canonical_to_assay.get(g,[''])[0],measured=g in canonical_to_assay)for k in functional+['Reactome_UPR','GO_ER_UPR']for g in sorted(sets[k])]
write_tsv(OUT/'authoritative_full_memberships.tsv',base_rows,list(base_rows[0]))
overlap=[]
for left in functional+['Reactome_UPR','GO_ER_UPR','Hallmark_UPR']:
    a=hallmark_canonical if left=='Hallmark_UPR'else sets[left]
    for right in functional+['Reactome_UPR','GO_ER_UPR','Hallmark_UPR']:
        b=hallmark_canonical if right=='Hallmark_UPR'else sets[right]
        overlap.append(dict(left=left,right=right,left_members=len(a),right_members=len(b),
            overlap=len(a&b),jaccard=len(a&b)/len(a|b),measured_overlap=sum(g in canonical_to_assay for g in a&b),
            left_subset_of_right=a<=b,right_subset_of_left=b<=a))
write_tsv(OUT/'hierarchy_and_overlap.tsv',overlap,list(overlap[0]))
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
protocol=f'''Membership lock: {stamp}
Definitions use authoritative GO human annotations retrieved from QuickGO, ontology release 2026-09-27, annotation snapshot 2026-07-28 10:08. Propagation includes is_a and part_of only, excluding regulates. Only UniProtKB protein gene products are included; NOT annotations are excluded; all remaining evidence classes, including IEA, are retained. RNAcentral RNA and ComplexPortal complex records do not define individual protein-coding gene readouts.
Gene products are mapped to current approved HGNC genes by unique HGNC UniProt mapping, exact approved symbol, then unique previous/alias symbol. Ambiguous or unmapped products are retained in mapping tables but do not generate invented genes. Current assay symbols are mapped by the same exact/unique-alias procedure. No canonical gene has multiple columns within the fixed assay universe. The complete HGNC snapshot and all raw API results are retained with hashes.
The independent Reactome parent is R-HSA-381119.4, using its recursively expanded participants' UniProt ReferenceGeneProduct entities. GO:0030968 is an ontology sensitivity. These sets are independent definitions of the pre-existing Hallmark readout, not evidence of a newly established biochemical mechanism.
Functional modules are ER protein folding (GO:0034975), ERAD quality control pathway (GO:0036503), protein targeting to ER (GO:0045047), and integrated stress response signaling (GO:0140467). ER targeting is specifically signal-directed co/post-translational targeting and ER import machinery; it does not include the general ER-to-Golgi vesicle-transport process. ISR includes PERK, HRI, GCN2 and PKR ontology descendants. Expression of their members is not a measurement of protein-folding, ERAD, translocation or kinase activity.
Full authoritative memberships are primary. Hallmark-disjoint variants remove all original113 Hallmark-UPR canonical genes before measurement intersection. Module-unique variants remove every member of the other three full functional modules before intersection; no genes are added to compensate. Every empty/single/few-member set is retained and labeled; no coverage threshold is tuned on outcomes.
Normalization uses the unchanged6780 measured-gene CP10k denominator. Metadata-only intersection of the four genetic studies has6783 genes; the extra3 do not alter this primary denominator. Primary real readouts exclude each target's measured raw count from the denominator and exclude that target from the numerator if a member, in treated and target-specific matched controls. Fixed-denominator readouts with the same numerator exclusion provide a separate sensitivity and the fair calibration against matched random modules. The original Hallmark90 readout is rescored in both modes. No approximate nonlinear denominator rescaling is used.
RNA gene-set scores are equal-weight mean log1pCP10k of remaining measured members. Empty post-exclusion membership or nonpositive/nonfinite target-specific controlSD is unestimable; fewer than3 remaining members is flagged as a few-transcript proxy. All support failures remain in outputs. Fixed panels are438 K562/RPE1 genes,30 Adamson genes,19 Frangieh Control/IFNg genes and17 co-culture genes. These memberships may also annotate target-function strata descriptively; expression readouts and perturbation-target annotations are separate objects.
All functional, parent and sensitivity definitions above are fixed without inspecting their new validation associations. Source-only control-expression/detection matching for random modules follows this definition lock and is frozen before their response associations. The inference families and unavailable-test handling are specified independently in ../specificity/inference_protocol.txt.
'''
(OUT/'membership_protocol.txt').write_text(protocol)
files=['locked_memberships.tsv','locked_axis_definitions.tsv','authoritative_full_memberships.tsv','assay_HGNC_mapping.tsv','hierarchy_and_overlap.tsv','original_Hallmark_UPR_full.tsv','primary_gene_universe.tsv','membership_protocol.txt','freeze_memberships.py','fetch_authoritative.py','GO_annotation_mapping.tsv','Reactome_protein_reference_mapping.tsv']
files+= [str(p.relative_to(OUT))for p in sorted(RAW.iterdir())if p.is_file()]
manifest=[dict(file=f,sha256=hashlib.sha256((OUT/f).read_bytes()).hexdigest(),bytes=(OUT/f).stat().st_size)for f in files]
write_tsv(OUT/'membership_lock_manifest.tsv',manifest,['file','sha256','bytes'])
print(stamp)
for row in definitions_rows:
    print(row['axis'],row['full_canonical_members'],row['measured_members'])
print('manifest_sha256',hashlib.sha256((OUT/'membership_lock_manifest.tsv').read_bytes()).hexdigest())
