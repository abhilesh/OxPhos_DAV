#!/bin/bash
set -uo pipefail
JOBID=1840338
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"
cd /Users/ad2347/Documents/OxPhos_DAV

echo "=== Monitoring SLURM job $JOBID on lambda ==="
for i in $(seq 1 270); do
  STATE=$(ssh -o BatchMode=yes -o ConnectTimeout=10 lambda "squeue -j $JOBID -h -o %T" 2>/dev/null)
  if [ -z "$STATE" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): job $JOBID left the queue (finished or failed)"
    break
  fi
  echo "$(date '+%Y-%m-%d %H:%M:%S'): job $JOBID state=$STATE"
  sleep 60
done

echo ""
echo "=== Final job accounting ==="
ACCT=$(ssh -o BatchMode=yes lambda "sacct -j $JOBID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS --noheader -P 2>&1")
echo "$ACCT"

# Only look at the main job line (not the .batch/.extern sub-steps) for the gate decision
MAIN_STATE=$(echo "$ACCT" | awk -F'|' -v jid="$JOBID" '$1==jid {print $3}')
echo ""
echo "=== Gate: main job state = '$MAIN_STATE' ==="

echo ""
echo "=== Tail of stdout/stderr logs ==="
ssh -o BatchMode=yes lambda "tail -n 80 $HPC_WORK/logs/partners_rerun_${JOBID}.out 2>&1"
echo "--- stderr ---"
ssh -o BatchMode=yes lambda "tail -n 80 $HPC_WORK/logs/partners_rerun_${JOBID}.err 2>&1"

if [ "$MAIN_STATE" != "COMPLETED" ]; then
  echo ""
  echo "=== JOB DID NOT COMPLETE SUCCESSFULLY (state=$MAIN_STATE) — NOT touching local files, NOT running any comparison ==="
  exit 1
fi

echo ""
echo "=== Job completed. Verifying remote output files actually exist and are fresh before copying ==="
ssh -o BatchMode=yes lambda "ls -la $HPC_WORK/results/structural/compensatory_partners.csv $HPC_WORK/results/structural/all_tested_pairs.csv $HPC_WORK/results/structural/concordance_summary.csv 2>&1"

echo ""
echo "=== Copying results back ==="
rsync -av lambda:$HPC_WORK/results/structural/compensatory_partners.csv \
          lambda:$HPC_WORK/results/structural/all_tested_pairs.csv \
          lambda:$HPC_WORK/results/structural/concordance_summary.csv \
          results/structural/ 2>&1

echo ""
echo "=== Comparison: old (backed up 2026-05-19) vs new ==="
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
echo "=== DONE — pre-fix backups still at /tmp/compensatory_partners_backup_20260519/ pending explicit cleanup ==="
