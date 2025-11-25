# Population Genetics Test Files

Files in this folder as primarily intended to test population genetics applications.
It also contains a phased panel of reference that can be used for imputation and phasing.

## Simulated data

The data have been simulated using PLINK with the following command

```{bash}
plink --simulate nfcore_plink_params.sim --simulate-ncases 100 --simulate-ncontrols 100 --out plink_simulated
```

The parameter file was created with the following characteristics

```{bash}
200 null 0.05 0.8 1 1
20 dis 0.01 0.08 2 3
```

Details about the format of this file can be found [here](https://www.cog-genomics.org/plink/1.9/input#simulate).

## Format conversion

The files created with the previous command are generated in PLINK binary format.
In order to convert them to VCF the following command has been used

```{bash}
plink --bfile plink_simulated --recode vcf bgz --out plink_simulated
```

And to convert them to BCF the following command has been used

```{bash}
bcftools view -O b -o plink_simulated.bcf.gz plink_simulated.vcf.gz
```

And to convert them to PLINK 2 binary the following command has been used

```{bash}
plink2 --bfile plink_simulated --make-pgen --out plink_simulated
```

## Phased dataset

This folder contains the data from the [1000 Genome Project](https://www.internationalgenome.org/).
The data are phased and can be used as a panel of normal.

The data has been generated as follow :

```bash
# Download phased panel for chr21 - chr22
for chr in chr21 chr22; do

    for ext in "" ".tbi"; do
        wget -c \
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_$chr.filtered.shapeit2-duohmm-phased.vcf.gz${ext}" \
        -O "data/genomics/homo_sapiens/popgen/1000GP.${chr}.full.vcf.gz${ext}"
    done

    PANEL_FILE=data/genomics/homo_sapiens/popgen/1000GP.$chr
    REGION=$chr:16570000-16610000
    # Normalize, select the region, keep only biallelic SNPs, rename variants ID and compute allele frequency
    bcftools norm -m +any ${PANEL_FILE}.full.vcf.gz \
            --regions ${REGION} --threads 4 -Ov | \
        bcftools view \
            -m 2 -M 2 -v snps --threads 4 -Ov | \
        bcftools annotate --threads 4 -Ov \
            --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' | \
        vcffixup - | bgzip -c > ${PANEL_FILE}.vcf.gz

    # Index the panel file
    bcftools index -f ${PANEL_FILE}.vcf.gz --threads 4

    # Convert to hap legend samples file set
    bcftools convert --haplegendsample 1000GP.$chr ${PANEL_FILE}.vcf.gz -Oz -o ${PANEL_FILE}

    # Convert to posfile
    zcat ${PANEL_FILE}.legend.gz | awk 'NR>1 {
        split($1,a,/[:_]/);
        printf "%s\t%s\t%s\t%s\n", a[1], a[2], a[3], a[4]
    }' > ${PANEL_FILE}.posfile

    # Keep only variants informations
    bcftools view -G -m 2 -M 2 -v snps ${PANEL_FILE}.vcf.gz -Oz -o ${PANEL_FILE}.sites.vcf.gz
    bcftools index -f ${PANEL_FILE}.sites.vcf.gz

    # Chunk the panel for easy parallelization
    GLIMPSE_chunk \
        --input ${PANEL_FILE}.vcf.gz --region ${REGION} \
        --window-size 10000 --window-count 400 --buffer-size 5000 --buffer-count 30 \
        --output ${PANEL_FILE}.chunks.txt

    rm ${PANEL_FILE}.full.vcf.gz*
done
```
