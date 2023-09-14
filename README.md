# test-datasets: rangeland

This branch contains test data to be used for automated testing with the nf-core/rangeland pipeline.

## Data description

The dataset comprises a collection of Landsat data
derived based on the [Landsat Collection 2](https://www.usgs.gov/landsat-missions/landsat-collection-2) scenes native to the 181036 Landsat tile and acquired
between 01-01-1987 and 31-12-1989. This data is located at `Landsat_collection2/`.

Additionally, the dataset contains the followind data:

- `dem/`: digital elevation model derived from [copernicus](https://www.copernicus.eu/en)
- `endmember/`: endmember definition, obtained by [Hostert et al.](<https://doi.org/10.1016/S0034-4257(03)00145-7>)
- `wvp/`: water vapor database obtained [here](https://zenodo.org/record/4468701)
- `vector/`: vector data representing the targeted area
- `datacube/`: datacube definition
- `reference/` reference data for pipeline results

The data included in this repository conforms with the naming conventions of their respective products.
However, the spatial extend of each dataset is only a spatial subset of the original scene.
As such, the datasets in question should not be treated as 'real data' and _must not_
be used for research activities going beyond workflow testing.

The authors of this dataset accepts no responsibility for errors or omissions in this work
and shall not be liable for any damage caused by these.

## Funding Information

This dataset was developed as a part of research activity supported by  
the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project-ID 414984028 – SFB 1404 FONDA.

# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1.  [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2.  [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
