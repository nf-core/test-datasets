# Data for funcprofiler

## humann

In order to understand the database compatibility network, the following can be run:
```
for tag in 3.6.1 3.7 v3.8 v3.9 4.0.0.alpha.1 4.0.0.alpha.1-final ; do git checkout $tag ; cat humann/config.py | grep -e "metaphlan_v4_db_version\|metaphlan_v3_db_version\|metaphlan_v3_db_matching_uniref\|matching_uniref" | sed "s|^|$tag\t|g" ; done  > version_table
```


- Clone the humann repo

```
mkdir -p ./data/database/humann/v3/
mkdir -p ./data/database/humann/v3

checkout 3.6.1 in human repo
cp -r <path to clone>/humann/humann/data/ ./data/database/humann/v3
rm -r data/database/humann/v3/misc
mv ./data/database/humann/v3/chocophlan_DEMO/ ./data/database/humann/v3/chocophlan_nfDEMO/
mv ./data/database/humann/v3/uniref_DEMO/ ./data/database/humann/v3/uniref_nfDEMO/


tar czf data/database/humann/v3/chocophlan_nfDEMO.tar.gz data/database/humann/v3/chocophlan_nfDEMO
tar czf data/database/humann/v3/uniref_nfDEMO.tar.gz data/database/humann/v3/uniref_nfDEMO


# this raw file is 66mb so we subset it to test data size
#mv data/database/humann/pathways/metacyc_reactions_level4ec_only.uniref.bz data/database/humann/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full
#bzcat data/database/humann/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full | head -n 250  | bzip2 > data/database/humann/pathways/metacyc_reactions_level4ec_only.uniref.bz2
#rm data/database/humann/pathways/metacyc_reactions_level4ec_only.uniref.bz2_full
#tar czf data/database/humann.tar.gz data/database/humann/ && mv data/database/humann.tar.gz data/database/humann/


git checkout 4.0.0.alpha.1-final
cp -r ../humann-orig/humann/data/ ./data/database/humann/v4/
mv  ./data/database/humann/v4/chocophlan_DEMO/ ./data/database/humann/v4/chocophlan_nfDEMO/
mv ./data/database/humann/v4/uniref_DEMO/ ./data/database/humann/v4/uniref_nfDEMO/
mv data/database/humann/v4/utility_DEMO/ data/database/humann/v4/utility_nfDEMO/
tar czf data/database/humann/v4/utility.tar.gz data/database/humann/v4/utility_nfDEMO/
tar czf data/database/humann/v4/nucleotide.tar.gz data/database/humann/v4/chocophlan_nfDEMO/
tar czf data/database/humann/v4/protein.tar.gz data/database/humann/v4/uniref_nfDEMO/


```


## fmhfunprofiler

```
mkdir data/database/fmhfunprofiler/

wget https://zenodo.org/records/10045253/files/KOs_sketched_scaled_1000.sig.zip

unzip KOs_sketched_scaled_1000.sig.zip

# interupt after a few seconds, we only need the first 500
zip data/database/fmhfunprofiler/KOs_sketched_scaled_1000_demo.sig.zip signatures/*

```


## RGI

```
wget https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2

bunzip2 -c <  broadstreet-v4.0.1.tar.bz2 | gzip -c > broadstreet-v4.0.1.tar.gz

mkdir ./data/database/rgi
mv broadstreet-v4.0.1.tar.gz data/database/rgi/



```
