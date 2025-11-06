#!/bin/bash
set -e
# 需要EMBOSS，python和Bio包
# 使用方法：1.cds_cleaner.sh your_cds.fasta
module load EMBOSS
module load python

input="$1"
if [ -z "$input" ]; then
    echo "Usage: $0 cds.fasta"
    exit 1
fi

base="${input%.*}"

echo ">> Translating CDS → Protein"
transeq -sequence "$input" -outseq "${base}.pep.tmp" -table 1 -frame 1

echo ">> Filtering"
python3 - "$input" "${base}.pep.tmp" <<'PY'
import sys
from Bio import SeqIO
from Bio.Seq import Seq

cds_file = sys.argv[1]
pep_file = sys.argv[2]

cds = {r.id: str(r.seq).upper() for r in SeqIO.parse(cds_file,"fasta")}
pep = {r.id: str(r.seq) for r in SeqIO.parse(pep_file,"fasta")}

valid_ids = []
bad_len = []
bad_stop = []
bad_x = []

for cid, seq in cds.items():
    # must multiple of 3
    if len(seq) % 3 != 0:
        bad_len.append(cid)
        continue

    prot = str(Seq(seq).translate())

    # internal stop
    if "*" in prot[:-1]:
        bad_stop.append(cid)
        continue

    # X
    if "X" in prot:
        bad_x.append(cid)
        continue

    valid_ids.append(cid)

print("GOOD:", len(valid_ids))
print("len%3!=0:", len(bad_len))
print("internal_stop:", len(bad_stop))
print("has X:", len(bad_x))

with open("valid_ids.lst","w") as f:
    f.write("\n".join(valid_ids))
PY

echo ">> Extracting valid sequences"
seqkit grep -f valid_ids.lst "$input" > "${base}.clean.cds.fasta"
seqkit grep -f valid_ids.lst "${base}.pep.tmp" > "${base}.clean.pep.fasta"

echo "DONE:"
echo "  ${base}.clean.cds.fasta"
echo "  ${base}.clean.pep.fasta"
