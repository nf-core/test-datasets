# test-datasets: rangeland

This branch contains test data to be used for automated testing with the nf-core/rangeland pipeline.

## Data description

The dataset comprises a collection of Landsat data
derived based on the [Landsat Collection 2](https://www.usgs.gov/landsat-missions/landsat-collection-2) scenes native to the 181036 Landsat tile and acquired
between 01-01-1987 and 31-12-1989.

Many remote sensing tools rely on information encoded in the file structure. This workflow is no exception and thus, we provide the `landsat_dem_wvp.tar.gz` archive. This archive contains:

- `Landsat_collection2/`: Landsat collection 2 imagery.
- `dem/`: digital elevation model derived from [copernicus](https://www.copernicus.eu/en).
- `wvp/`: water vapor database obtained [here](https://zenodo.org/record/4468701).

Additionally, the dataset contains the following data:

- `endmember/`: endmember definition, obtained by [Hostert et al.](<https://doi.org/10.1016/S0034-4257(03)00145-7>).
- `vector/`: vector data representing the area of interest.
- `datacube/`: datacube definition.
- `reference/` reference data for pipeline results.

The data included in this repository conforms with the naming conventions of their respective products.
However, the spatial extend of each dataset is only a spatial subset of the original scene.
As such, the datasets in question should not be treated as 'real data' and _must not_
be used for research activities going beyond workflow testing.

The authors of this dataset accepts no responsibility for errors or omissions in this work
and shall not be liable for any damage caused by these.

## Funding Information

This dataset was developed as a part of research activity supported by  
the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project-ID 414984028 – SFB 1404 FONDA.
