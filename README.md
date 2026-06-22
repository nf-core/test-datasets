# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

# test-datasets: stableexpression
This branch contains test data to be used for automated testing with the nf-core/stableexpression pipeline.

## Content of this repository

All the data contained here were subsampled from datasets collected from Expression. In some cases, data were also generated randomly.

<details markdown="1">
<summary>Folders</summary>

- `test_data/`: contains files for different pipeline use cases
  - `drosophila_simulans/`: two real life datasets (produced in our lab) for drosophila simulans, used for testing that the pipeline output matches the expected results from the publication to come
  - `mus_musculus/`: a link to a dataset used in the `nf-core/differentialabundance` pipeline, together with its design
  - `sampled/`: multiple mini-datasets used for testing (real gene ids but random values)
    - `solanum_tuberosum`
    - `beta vulgaris`
    - `prunus_persica`

</details>
