# References created for sarek

Using the nf-core/references pipeline

```bash
nextflow run nf-core/references -profile docker \
    -r 8112ae8 \
    --input https://raw.githubusercontent.com/nf-core/references-assets/41545a3631addaf491d22751b17607149b8512ac/assets/test/sarek/GRCh38_chr22.yml \
    --tools bwamem1,createsequencedictionary,faidx,intervals,tabix \
    --outdir references
```
