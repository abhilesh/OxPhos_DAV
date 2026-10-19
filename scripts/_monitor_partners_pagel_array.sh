#!/bin/bash
set -uo pipefail
JOBID=1840345
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"
JOBS_DIR="$HPC_WORK/data/phylo/partners_pagel_jobs"
RESULTS_DIR="$HPC_WORK/results/structural/partners_pagel_results"
cd /Users/ad2347/Documents/OxPhos_DAV

echo "=== Monitoring SLURM array job $JOBID (60 tasks) on lambda ==="
for i in $(seq 1 270); do
  N_LEFT=$(ssh -o BatchMode=yes -o ConnectTimeout=10 lambda "squeue -j $JOBID -h -o %T" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$N_LEFT" = "0" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): all array tasks left the queue"
    break
  fi
  echo "$(date '+%Y-%m-%d %H:%M:%S'): $N_LEFT / 60 tasks still in queue"
  sleep 90
done

echo ""
echo "=== Per-task final states ==="
ACCT=$(ssh -o BatchMode=yes lambda "sacct -j $JOBID --format=JobID,State,ExitCode,Elapsed --noheader -P 2>&1")
echo "$ACCT" | grep -v '\.batch\|\.extern'

N_COMPLETED=$(echo "$ACCT" | grep -v '\.batch\|\.extern' | awk -F'|' '$2=="COMPLETED"' | wc -l | tr -d ' ')
N_FAILED=$(echo "$ACCT" | grep -v '\.batch\|\.extern' | awk -F'|' '$2=="FAILED"' | wc -l | tr -d ' ')
echo ""
echo "=== Summary: $N_COMPLETED completed, $N_FAILED failed (of 60) ==="

echo ""
echo "=== Result files actually present on WORK ==="
ssh -o BatchMode=yes lambda "ls $RESULTS_DIR/ 2>&1 | wc -l"
ssh -o BatchMode=yes lambda "wc -l $RESULTS_DIR/results_*.tsv 2>&1 | tail -5"

if [ "$N_COMPLETED" = "0" ]; then
  echo "=== NO tasks completed successfully — not copying anything back or merging ==="
  exit 1
fi

echo ""
echo "=== Copying results back ==="
mkdir -p results/structural/partners_pagel_results
rsync -av lambda:$RESULTS_DIR/ results/structural/partners_pagel_results/ 2>&1 | tail -10

echo ""
echo "=== Merging (Docker) ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/structural/01_find_compensating_partners.py --merge-pagel-jobs \
    --jobs-dir data/phylo/partners_pagel_jobs --results-dir results/structural/partners_pagel_results 2>&1

echo ""
echo "=== Comparison: old (2026-05-19 backup) vs newly merged ==="
python3 -c "
import csv

def load(path):
    with open(path) as f:
        return list(csv.DictReader(f))

old = load('/tmp/compensatory_partners_backup_20260519/compensatory_partners.csv')
new = load('results/structural/compensatory_partners.csv')
old_all = load('/tmp/compensatory_partners_backup_20260519/all_tested_pairs.csv')
new_all = load('results/structural/all_tested_pairs.csv')

print(f'compensatory_partners.csv: old={len(old)} rows, new={len(new)} rows')
print(f'all_tested_pairs.csv:      old={len(old_all)} rows, new={len(new_all)} rows')

def sig_counts(rows, col):
    return sum(1 for r in rows if r.get(col) not in (None, '', 'None') and float(r[col]) <= 0.10)

def nonnull_counts(rows, col):
    return sum(1 for r in rows if r.get(col) not in (None, '', 'None'))

for col in ('fisher_fdr', 'pagel_fdr', 'branch_cooccur_fdr'):
    o = sig_counts(old_all, col)
    n = sig_counts(new_all, col)
    on = nonnull_counts(old_all, col)
    nn = nonnull_counts(new_all, col)
    print(f'{col}: non-null old={on} new={nn} | <=0.10 old={o} new={n}')

targets = [('MT-ATP6','192','189'), ('SDHC','146','144'), ('SDHC','46','47'), ('COX6C','24','25')]
print()
print('Named candidate check (dar_gene, dar_pos, contact_pos):')
for gene,pos,cpos in targets:
    m = [r for r in new if r.get('dar_gene')==gene and r.get('dar_aa_coord')==pos and r.get('contact_refseq_pos')==cpos]
    print(f'  {gene}:{pos}<->{cpos}  {\"FOUND\" if m else \"still NOT present\"} in new compensatory_partners.csv')
"

echo ""
echo "=== DONE ==="
