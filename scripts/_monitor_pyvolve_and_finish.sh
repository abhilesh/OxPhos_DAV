#!/bin/bash
set -uo pipefail
JOBID=1840740
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"
cd /Users/ad2347/Documents/OxPhos_DAV

echo "=== Monitoring SLURM array job $JOBID (50 tasks, Pyvolve permutation) on lambda ==="
for i in $(seq 1 480); do
  N_LEFT=$(ssh -o BatchMode=yes -o ConnectTimeout=10 lambda "squeue -j $JOBID -h -o %T" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$N_LEFT" = "0" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): all array tasks left the queue"
    break
  fi
  echo "$(date '+%Y-%m-%d %H:%M:%S'): $N_LEFT / 50 tasks still in queue"
  sleep 90
done

echo ""
echo "=== Per-task final states ==="
ACCT=$(ssh -o BatchMode=yes lambda "sacct -j $JOBID --format=JobID,State,ExitCode,Elapsed --noheader -P 2>&1")
echo "$ACCT" | grep -v '\.batch\|\.extern'

N_COMPLETED=$(echo "$ACCT" | grep -v '\.batch\|\.extern' | awk -F'|' '$2=="COMPLETED"' | wc -l | tr -d ' ')
N_FAILED=$(echo "$ACCT" | grep -v '\.batch\|\.extern' | awk -F'|' '$2=="FAILED"' | wc -l | tr -d ' ')
echo ""
echo "=== Summary: $N_COMPLETED completed, $N_FAILED failed (of 50) ==="

if [ "$N_COMPLETED" = "0" ]; then
  echo "=== NO tasks completed successfully — stopping, not touching local files ==="
  exit 1
fi

echo ""
echo "=== Copying chunk results back ==="
rsync -av "lambda:$HPC_WORK/results/phylo/conditional_permissiveness_chunk*.csv" results/phylo/ 2>&1 | tail -20
N_CHUNKS_BACK=$(ls results/phylo/conditional_permissiveness_chunk*.csv 2>/dev/null | wc -l | tr -d ' ')
echo "Chunks copied back: $N_CHUNKS_BACK"

if [ "$N_CHUNKS_BACK" = "0" ]; then
  echo "=== No chunk files found — stopping ==="
  exit 1
fi

echo ""
echo "=== Merging Pyvolve chunks (05_merge_perm_chunks.py) ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/phylo/05_merge_perm_chunks.py 2>&1

if [ ! -f results/phylo/conditional_permissiveness.csv ]; then
  echo "=== conditional_permissiveness.csv not produced — stopping ==="
  exit 1
fi
echo ""
echo "=== conditional_permissiveness.csv row count and mtime ==="
wc -l results/phylo/conditional_permissiveness.csv
ls -la results/phylo/conditional_permissiveness.csv

echo ""
echo "=== Running 06b_temporal_ordering_mitonuc.py (fixed baseline, fresh conditional_permissiveness.csv) ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/phylo/06b_temporal_ordering_mitonuc.py 2>&1

echo ""
echo "=== Running 07_compartment_analysis.py ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/phylo/07_compartment_analysis.py 2>&1

echo ""
echo "=== Running 09_pyvolve_rate_controlled.py ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/phylo/09_pyvolve_rate_controlled.py 2>&1

echo ""
echo "=== Running 11_tier1_extraction.py ==="
docker run --rm -v "$(pwd)":/app oxphos_dav_analysis conda run -n oxphos_dav \
    python src/phylo/11_tier1_extraction.py 2>&1

echo ""
echo "=== Final results/phylo/ file listing ==="
ls -la results/phylo/*.csv

echo ""
echo "=== DONE ==="
