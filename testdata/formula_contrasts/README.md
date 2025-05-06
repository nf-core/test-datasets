## Experimental Contrasts

- The file [`SRP254919.contrasts_string.yaml`](./SRP254919.contrasts_string.yaml) defines contrast configurations for analyzing the effect of treatment between *mCherry* and *hND6*, using either formula or comparison, also testing blocking variables. This file is based on `SRP254919.contrasts.yaml`.

- The contrasts defined in [`rnaseq_complex_contrast.yaml`](./rnaseq_complex_contrast.yaml) are used to test and compare different model formulas and contrast strings for RNA-seq differential expression analysis. Each contrast represents a distinct way of modeling the data, ranging from simple treatment effects to complex interaction terms involving genotype and time.
