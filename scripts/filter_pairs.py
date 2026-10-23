"""Keep the read pairs listed in a file. Usage: filter_pairs.py names.txt in_1 in_2 out_1 out_2"""
import gzip, sys
names, in1, in2, out1, out2 = sys.argv[1:]
keep = set(line.strip() for line in open(names))
kept = total = 0
with gzip.open(in1, "rt") as i1, gzip.open(in2, "rt") as i2, gzip.open(out1, "wt", compresslevel=9) as o1, gzip.open(out2, "wt", compresslevel=9) as o2:
    while True:
        r1 = [i1.readline() for _ in range(4)]
        r2 = [i2.readline() for _ in range(4)]
        if not r1[0]:
            break
        total += 1
        if r1[0][1:].split()[0] in keep:
            o1.writelines(r1)
            o2.writelines(r2)
            kept += 1
print(out1, kept, "of", total, file=sys.stderr)
