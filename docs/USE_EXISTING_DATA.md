# How to use an existing test dataset

Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things requested when adding a new test dataset.

 - [ ] Check [here](https://github.com/nf-core/test-datasets/branches/all) to find the branch corresponding to the test dataset you want to use (i.e., `modules`, or a given pipeline)
   - [ ] Alternatively use nf-core tools to search and explore, e.g. `nf-core test-datasets list` (press tab for options), and then `nf-core test-datasets search -b modules fasta` (press tab for options)
   - [ ] You can also use the flag `--generate-nf-path` to generate config ready strings
 - [ ] Specify in the `conf/test.config` the path to the files from the test dataset 
 - [ ] Set up your CI tests following the nf-core best practices
   - [ ] For pipelines: [tests/default.nf.test](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/tests/default.nf.test))
   - [ ] For modules: [tests/main.nf.test](https://github.com/nf-core/tools/blob/main/nf_core/module-template/tests/main.nf.test.j2)
