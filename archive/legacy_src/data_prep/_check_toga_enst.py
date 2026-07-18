import sys, re, urllib.request, ssl, json, gzip
sys.path.insert(0,'src')
from pathlib import Path
from collections import defaultdict, Counter

DATA_DIR = Path('data')
REF_DIR = DATA_DIR / 'reference'
CURATED_DIR = DATA_DIR / 'annotations' / 'curated'
BASE_URL = 'https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments'

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

req = urllib.request.Request(f'{BASE_URL}/')
with urllib.request.urlopen(req, context=ctx) as r:
    html = r.read().decode()

pattern = re.compile(r'href="(ENST\d+)\.([^.]+)\.fasta\.gz"')
toga_by_gene = defaultdict(list)
for enst, gene in pattern.findall(html):
    toga_by_gene[gene].append(enst)

mane_map = {}
with gzip.open(REF_DIR / 'MANE_GRCh38_v1.5.txt.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.split('\t')
        if len(cols) < 8: continue
        nm = cols[5].strip().split('.')[0]
        enst = cols[7].strip().split('.')[0]
        if nm.startswith('NM_') and enst.startswith('ENST'):
            mane_map[nm] = enst

nm_counts = defaultdict(Counter)
with open(CURATED_DIR / 'nucDNA_annotations.json') as f:
    variants = json.load(f)
for v in variants:
    gene = v.get('locus','')
    tx = v.get('transcript_id','').split('.')[0]
    if gene and tx.startswith('NM_'):
        nm_counts[gene][tx] += 1
clinvar_nm = {g: counts.most_common(1)[0][0] for g, counts in nm_counts.items()}

problem_genes = ['NDUFS1','NDUFV2','NDUFA10','NDUFS7','NDUFA6','NDUFS2','NDUFS6','ATP5MC2','NDUFB11','COX7A2']
print(f"{'Gene':<15} {'MANE ENST':<20} {'On TOGA?':<10} {'All TOGA ENSTs'}")
print('-'*80)
for gene in problem_genes:
    nm = clinvar_nm.get(gene,'?')
    mane_enst = mane_map.get(nm,'?')
    available = toga_by_gene.get(gene, [])
    on_server = mane_enst in available
    print(f"{gene:<15} {mane_enst:<20} {str(on_server):<10} {available}")
