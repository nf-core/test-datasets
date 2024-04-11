# Population Genetics Test Files

Files in this folder as primarily intended to test population genetics applications.

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
