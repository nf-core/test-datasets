#!/usr/bin/env python3
"""Build a tiny Bakta DB for CI by subsetting a light DB to the entries hit by test genomes.

Usage: build_bakta_db_subset.py <light_db_dir> <out_dir> <run_dir>...
Each run_dir is a `bakta --debug --tmp-dir <run>.tmp --output <run>` run against the light DB;
both <run> (json/log) and <run>.tmp (raw hit tables) are read.
"""
import re
import shutil
import sqlite3
import subprocess as sp
import sys
from pathlib import Path

import pyhmmer

light, out = Path(sys.argv[1]), Path(sys.argv[2])
runs = [Path(r) for r in sys.argv[3:]]
tmps = [t for r in runs for t in Path(f'{r}.tmp').iterdir()]
assert runs and tmps, 'no bakta runs given'


def cols(name, col, skip_comments=False):
    ids = set()
    for t in tmps:
        p = t / name
        if p.exists():
            for line in p.read_text().splitlines():
                if line and not (skip_comments and line.startswith('#')):
                    ids.add(re.split(r'\t' if not skip_comments else r'\s+', line)[col])
    return ids


def grep(pattern, glob):
    return {m for r in runs for f in r.glob(glob) for m in re.findall(pattern, f.read_text())}


def run(*cmd, **kw):
    sp.run(cmd, check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL, **kw)


def subset_fasta(src, dst, keep):
    """keep(header_line) -> bool"""
    with open(src) as fi, open(dst, 'w') as fo:
        write = False
        for line in fi:
            if line.startswith('>'):
                write = keep(line)
            if write:
                fo.write(line)


def subset_dmnd(name, ids):
    faa = out / f'{name}.faa'
    run('diamond', 'getseq', '-d', str(light / f'{name}.dmnd'), '--out', str(faa))
    subset_fasta(faa, faa.with_suffix('.sub.faa'), lambda h: h[1:].split()[0] in ids)
    run('diamond', 'makedb', '--in', str(faa.with_suffix('.sub.faa')), '-d', str(out / name))
    faa.unlink(), faa.with_suffix('.sub.faa').unlink()


def subset_hmm(name, keep):
    txt = out / name
    with pyhmmer.plan7.HMMFile(str(light / name)) as fi, open(txt, 'wb') as fo:
        n = 0
        for hmm in fi:
            if keep(hmm):
                hmm.write(fo, binary=False)
                n += 1
    assert n, f'no {name} models kept'
    run('hmmpress', '-f', str(txt))
    txt.unlink()


def subset_cm(name, accs):
    txt = out / name
    with open(txt, 'w') as fo:
        sp.run(['cmconvert', str(light / f'{name}.i1m')], check=True, stdout=fo)
    run('cmfetch', '--index', str(txt))
    keys = out / f'{name}.keys'
    keys.write_text('\n'.join(sorted(accs)) + '\n')
    with open(f'{txt}.sub', 'w') as fo:
        sp.run(['cmfetch', '-f', str(txt), str(keys)], check=True, stdout=fo)
    Path(f'{txt}.sub').replace(txt)
    Path(f'{txt}.ssi').unlink(), keys.unlink()
    run('cmpress', '-F', str(txt))
    txt.unlink()


out.mkdir(parents=True, exist_ok=True)
for f in ['version.json', 'rfam-go.tsv', 'orit.fna', 'sorf.dmnd'] + [f'rRNA.i1{x}' for x in 'fimp']:
    shutil.copy(light / f, out / f)  # already tiny

pscc = cols('diamond.pscc.tsv', 1)
subset_dmnd('pscc', pscc)
subset_dmnd('expert-protein-sequences', cols('diamond.cds.expert.tsv', 1))

db = out / 'bakta.db'
db.unlink(missing_ok=True)
con = sqlite3.connect(db)
con.execute(f"ATTACH '{light / 'bakta.db'}' AS src")
for (sql,) in con.execute("SELECT sql FROM src.sqlite_master WHERE type='table'").fetchall():
    con.execute(sql)
con.execute('CREATE TEMP TABLE keep (id TEXT PRIMARY KEY)')
con.executemany('INSERT INTO keep VALUES (?)', [(i,) for i in pscc | cols('diamond.sorf.tsv', 1)])  # sORF hits are pscc keys too
con.execute('INSERT INTO pscc SELECT * FROM src.pscc WHERE uniref50_id IN (SELECT id FROM keep)')
con.commit()
con.execute('DETACH src')
con.execute('VACUUM')
con.close()

pfams = grep(r'"(PF\d{5}\.\d+)"', '*.json')
subset_hmm('pfam', lambda h: h.accession in pfams)
antifams = grep(r'homology=(Spurious_ORF_\d+)', '*.log')
subset_hmm('antifam', lambda h: h.name in antifams or h.accession == 'ANF00001')  # keep >=1 model

subset_cm('ncRNA-genes', cols('ncrna-genes.tsv', 1, skip_comments=True))
subset_cm('ncRNA-regions', cols('ncrna-regions.tsv', 1, skip_comments=True))

# blastn query file: keep hits plus a few so the search still runs
oric_hits = cols('ori.blastn.tsv', 0)
n = 0
def keep_oric(h):
    global n
    n += 1
    return n <= 10 or h[1:].split()[0] in oric_hits
subset_fasta(light / 'oric.fna', out / 'oric.fna', keep_oric)

# AMRFinderPlus: keep the reference proteins/CDS and HMMs of the genes found, then reindex
amr_src = (light / 'amrfinderplus-db' / 'latest').resolve()
amr_dst = out / 'amrfinderplus-db' / amr_src.name
shutil.rmtree(out / 'amrfinderplus-db', ignore_errors=True)
shutil.copytree(amr_src, amr_dst, ignore=shutil.ignore_patterns('*.n??', '*.p??', '*.h3?', 'AMR.LIB', 'AMRProt.fa', 'AMR_CDS.fa'))
(out / 'amrfinderplus-db' / 'latest').symlink_to(amr_src.name)
amr_rows = [l.split('\t') for t in tmps if (t / 'amrfinder.tsv').exists() for l in (t / 'amrfinder.tsv').read_text().splitlines()[1:]]
amr_genes = {r[1] for r in amr_rows}
amr_refs = {r[14] for r in amr_rows}
fam = [l.split('\t') for l in (amr_src / 'fam.tsv').read_text().splitlines() if not l.startswith('#')]
amr_hmms = {r[3] for r in fam if r[2] in amr_genes and r[3] != '-'}
keep_amr = lambda h: h[1:].split('|')[0] in amr_refs or not amr_genes.isdisjoint(h[1:].split('|'))
subset_fasta(amr_src / 'AMRProt.fa', amr_dst / 'AMRProt.fa', keep_amr)
subset_fasta(amr_src / 'AMR_CDS.fa', amr_dst / 'AMR_CDS.fa', keep_amr)
with pyhmmer.plan7.HMMFile(str(amr_src / 'AMR.LIB')) as fi, open(amr_dst / 'AMR.LIB', 'wb') as fo:
    for hmm in fi:
        if hmm.accession in amr_hmms:
            hmm.write(fo, binary=False)
run('amrfinder_index', str(amr_dst))

print(f'pscc={len(pscc)} pfam={len(pfams)} antifam={len(antifams)} amr_genes={sorted(amr_genes)} amr_hmms={len(amr_hmms)}')
