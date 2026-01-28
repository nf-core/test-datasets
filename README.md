# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Test data and samplesheets for Sopa

Some samplesheets are available under the `samplesheet` directory:

- `toy.csv`: a toy samplesheet (under the hood, Sopa will create a synthetic dataset via [`sopa.io.toy_dataset`](https://gustaveroussy.github.io/sopa/api/readers/#sopa.io.toy_dataset)). This dataset is composed of two spatial elements (an image and a transcript dataframe) representing cells generated uniformly in a square.
