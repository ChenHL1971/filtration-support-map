"""Added for the code release. After run_reproduce_results.sh, compares every regenerated output with the SHA-256 of the
file that produced the reported results (docs/expected_outputs.sha256). Figure files are compared too; they match only
with the matplotlib version in requirements.txt.
A documented known difference is accepted only if the regenerated file matches the exact SHA-256 recorded for that
difference (docs/known_differences.sha256); any other change to that file counts as unexpected.
Exit status: 0 if every file is identical or matches its recorded known difference; 1 otherwise.
Usage: python3 scripts/verify_reproduction.py"""
import hashlib, os, sys
NOTE = {"out/aim1/aim2_mgi.json": "upper bootstrap bound of the unreported MGI L1 difference is -0.0002, saved -0.0003 (README, Known differences)"}
known_hash = dict(reversed(l.split()) for l in open("docs/known_differences.sha256") if l.strip())
identical, known, bad = 0, 0, 0
for line in open("docs/expected_outputs.sha256"):
    h, p = line.split()
    if not os.path.exists(p):
        print(f"MISSING           {p}"); bad += 1; continue
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if got == h:
        print(f"IDENTICAL         {p}"); identical += 1
    elif known_hash.get(p) == got:
        print(f"KNOWN DIFFERENCE  {p}  [{NOTE.get(p, 'documented')}]"); known += 1
    else:
        print(f"UNEXPECTED DIFF   {p}"); bad += 1
print(f"\n{identical} byte-identical; {known} known difference(s); {bad} unexpected difference(s) or missing file(s)")
sys.exit(1 if bad > 0 else 0)
