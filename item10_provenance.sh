#!/usr/bin/env bash
# Item 10: record the missing provenance details. Run on YOUR analysis machine (needs internet).
# Usage: bash item10_provenance.sh [saved_Q15116.fasta] [saved_Q9NZQ7.fasta]
echo "== Interpreter used for the phase scripts (confirm this is /usr/local/bin/python3, Python 3.14.7) =="
which python3; python3 --version
python3 -c "import Bio, numpy; print('Biopython', Bio.__version__, '| NumPy', numpy.__version__)"
echo; echo "== Shebang line of each script (shows which interpreter it asks for) =="
for f in *.py; do [ -f "$f" ] && echo "$f: $(head -1 "$f")"; done
echo; echo "== UniProt release currently served (record this, and the date you ran it) =="
date -u +"run at %Y-%m-%d %H:%M UTC"
for acc in Q15116 Q9NZQ7; do
  echo "-- $acc"
  curl -sI "https://rest.uniprot.org/uniprotkb/$acc.fasta" | grep -iE "x-uniprot-release"
  if [ -n "$1" ] && [ "$acc" = "Q15116" ]; then saved="$1"; elif [ -n "$2" ] && [ "$acc" = "Q9NZQ7" ]; then saved="$2"; else saved=""; fi
  if [ -n "$saved" ]; then
    curl -s "https://rest.uniprot.org/uniprotkb/$acc.fasta" > "/tmp/$acc.now.fasta"
    grep -v '^>' "$saved" > "/tmp/$acc.saved.seq"; grep -v '^>' "/tmp/$acc.now.fasta" > "/tmp/$acc.now.seq"
    if cmp -s "/tmp/$acc.saved.seq" "/tmp/$acc.now.seq"; then echo "sequence identical to your saved file"; else echo "SEQUENCE DIFFERS from your saved file"; fi
  fi
done
echo; echo "== Repository: after pushing, record this commit hash =="
git rev-parse HEAD 2>/dev/null || echo "(not inside a git repository yet)"
