# Data for funcprofiler

## humann

In order to understand the database compatibility network, the following can be run:
```
for tag in 3.6.1 3.7 v3.8 v3.9 4.0.0.alpha.1 4.0.0.alpha.1-final ; do git checkout $tag ; cat humann/config.py | grep -e "metaphlan_v4_db_version\|metaphlan_v3_db_version\|metaphlan_v3_db_matching_uniref\|matching_uniref" | sed "s|^|$tag\t|g" ; done  > version_table
```


- Clone the humann repo

```
mkdir -p ./data/database/humann/v3
mkdir -p ./data/database/humann/v4

#checkout 3.6.1 in human repo
cp -r <path to clone>/humann/humann/data/ ./data/database/humann/v3
rm -r data/database/humann/v3/pathways
# rename so it doesn't warn about us using demo databases
mv ./data/database/humann/v3/chocophlan_DEMO/ ./data/database/humann/v3/chocophlan_nfDEMO/
mv ./data/database/humann/v3/uniref_DEMO/ ./data/database/humann/v3/uniref_nfDEMO/

# remove excess stuff from the utility database, any file that is not in the full db at https://g-227ca.190ebd.75bc.data.globus.org/humann_data/full_mapping_v201901b.tar.gz
for f in map_EC_to_triplet_AC_U50_U90_Swissprot_and_Trembl.txt.gz KeggOrgId2OrgNameTable.txt map_infogo1000_uniref50.txt.gz map_infogo1000_uniref90.txt.gz  map_kegg-pwy_name.txt.gz   map_metacyc-pwy_name.txt.gz map_metacyc-rxn_name.txt.gz map_transporter_uniref50.txt.gz mpa_vJan21_CHOCOPhlAnSGB_202103.tsv uniref50-tol-lca.dat.gz uniref90-tol-lca.dat.gz ; do
rm ./data/database/humann/v3/misc/$f
done

# shrink the remaining
find ./data/database/humann/v3/misc/ -name "*.gz" | while read f ; do tmpfile="tmp_$(basename $f)"; mv $f $tmpfile ; zcat < $tmpfile | head -n 100 | gzip --best > $f ; rm $tmpfile ; done
find ./data/database/humann/v3/misc/ -name "*.bz2" | while read f ; do tmpfile="tmp_$(basename $f)"; mv $f $tmpfile ; bzcat < $tmpfile | head -n 100 | bzip2 --best > $f ; rm $tmpfile ; done



tar czf data/database/humann/v3/chocophlan_nfDEMO.tar.gz data/database/humann/v3/chocophlan_nfDEMO
tar czf data/database/humann/v3/uniref_nfDEMO.tar.gz data/database/humann/v3/uniref_nfDEMO
tar czf data/database/humann/v3/utility_nfDEMO.tar.gz data/database/humann/v3/misc

git checkout 4.0.0.alpha.1-final

cp -r ../humann-orig/humann/data/ ./data/database/humann/v4/
mv  ./data/database/humann/v4/chocophlan_DEMO/ ./data/database/humann/v4/chocophlan_nfDEMO/
mv ./data/database/humann/v4/uniref_DEMO/ ./data/database/humann/v4/uniref_nfDEMO/
mv data/database/humann/v4/utility_DEMO/ data/database/humann/v4/utility_nfDEMO/

# strategic shrinking of utility dbs
mv data/database/humann/v4/utility_nfDEMO/map_eggnog_uniref90.txt.gz data/database/humann/v4/utility_nfDEMO/map_eggnog_uniref90.txt.gz_full
zcat < data/database/humann/v4/utility_nfDEMO/map_eggnog_uniref90.txt.gz_full | head -n 500 | gzip --best > data/database/humann/v4/utility_nfDEMO/map_eggnog_uniref90.txt.gz
rm data/database/humann/v4/utility_nfDEMO/map_eggnog_uniref90.txt.gz_full

mv data/database/humann/v4/utility_nfDEMO/map_ko_uniref90.txt.gz data/database/humann/v4/utility_nfDEMO/map_ko_uniref90.txt.gz_full
zcat < data/database/humann/v4/utility_nfDEMO/map_ko_uniref90.txt.gz_full | head -n 500 | gzip --best > data/database/humann/v4/utility_nfDEMO/map_ko_uniref90.txt.gz
rm data/database/humann/v4/utility_nfDEMO/map_ko_uniref90.txt.gz_full


tar czf data/database/humann/v4/utility_nfDEMO.tar.gz data/database/humann/v4/utility_nfDEMO/
tar czf data/database/humann/v4/chocophlan_nfDEMO.tar.gz data/database/humann/v4/chocophlan_nfDEMO/
tar czf data/database/humann/v4/uniref_nfDEMO.tar.gz data/database/humann/v4/uniref_nfDEMO/


```


## fmhfunprofiler

```
mkdir data/database/fmhfunprofiler/

cd data/database/fmhfunprofiler/
wget https://zenodo.org/records/10045253/files/KOs_sketched_scaled_1000.sig.zip

unzip KOs_sketched_scaled_1000.sig.zip

# find the 30s ribosomal protein s8 for a nice essential gene to test with
mkdir -p tmp/signatures

cat SOURMASH-MANIFEST.csv | grep "K02994\|moltype" > tmp/SOURMASH-MANIFEST.csv
tail -n+2 tmp/SOURMASH-MANIFEST.csv |cut -f 1 -d, |while read x ; do cp $x tmp/signatures/$(basename $x) ; done

cd tmp

zip ../KOs_sketched_scaled_1000_demo.sig.zip SOURMASH-MANIFEST.csv signatures/*
cd ../../../../
```


## RGI

```
wget https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2

bunzip2 -c <  broadstreet-v4.0.1.tar.bz2 | gzip -c > broadstreet-v4.0.1.tar.gz

mkdir ./data/database/rgi
mv broadstreet-v4.0.1.tar.gz data/database/rgi/



```



<!-- # this raw file is 66mb so we subset it to test data size -->
<!-- mv data/database/humann/v3/pathways/metacyc_reactions_level4ec_only.uniref.bz2 data/database/humann/v3/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full -->
<!-- bzcat data/database/humann/v3/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full | head -n 250  | bzip2 > data/database/humann/v3/pathways/metacyc_reactions_level4ec_only.uniref.bz2 -->
<!-- rm data/database/humann/v3/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full -->
<!-- tar czf data/database/humann.tar.gz data/database/humann/ && mv data/database/humann.tar.gz data/database/humann/ -->
