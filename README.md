# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

This branch contains new test data and references for the 
[nf-core/rnadnavar](https://github.com/nf-core/rnadnavar) 
pipeline.

## Content

>*Any content in this repo is to be used for testing 
purposes only.*

- `data/tcrb` contains minified fastq files from the 
  [Texas Cancer Research Biobank (TCRB)](http://stegg.hgsc.bcm.edu/open.html). We thank the TCRB and the 
  Texas cancer patient who donated their samples to make 
  open research truly accesible to everyone[^1]. By using 
  this data you agree to not attempt to re-identify 
  participants. The data is only to be used for testing purposes.
- `reference/chr7_hg38` contains reference files for 
  the analysis of a region in chromosome 7.
- `resources/vep` contains files for vep cache. This 
  cache has been modified to be as small as possible, 
  thus shouldn't be used for actual annotation, only for 
  testing purposes.


## Downloading test data

Due the large number of large files in the test-datasets 
repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^2]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Support

For further information or help, don't hesitate to get 
in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging). There is 
a dedicated channel to ask about [#test-data](https://join.slack.com/share/enQtNjE3MDg2MDc3Mzc4Mi01MjBhNTdmYzYyNzEzZjExMGMxYTA1YWY4ZWVkMjc3YzI4ZmRiOTE2YWI5ZDYxYjU5OGU0NTMxMDZiNTg0MDZh) and for 
[#rnadnavar](https://join.slack.com/share/enQtNjE3MDg2MDc3Mzc4Mi01MjBhNTdmYzYyNzEzZjExMGMxYTA1YWY4ZWVkMjc3YzI4ZmRiOTE2YWI5ZDYxYjU5OGU0NTMxMDZiNTg0MDZh)

[^1]: See publication [here](https://www.nature.com/articles/sdata201610)
[^2]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
