# How to add new test data

<<<<<<< HEAD
- collect/generate your test dataset
- make sure your test data is working with the module(s) it is designed for
- make a fork of the `nf-core/test-dataset` repository
- add the test data to your fork (checkout and use the `modules` branch)
- make a pull request to the `modules` branch from your fork of `nf-core/test-datasets`
- once the pull request is accepted, you can change the paths to the test data for your module(s) on `nf-core/modules`
=======
Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things requested when adding a new test dataset.

 - [ ] Check [here](https://github.com/nf-core/test-datasets/branches/all) that there isn't already a branch containing data that could be used
   - If this is the case, follow the [documentation on how to use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)
 - [ ] Fork the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets) to your GitHub account
 - [ ] Create a new branch on your fork
 - [ ] Add your test dataset
   - [ ] If you clone it locally use `git clone <url> --branch <branch> --single-branch`
 - [ ] Make a PR on a new branch with a relevant name
 - [ ] Wait for the PR to be merged
 - [ ] Use this newly created branch for your tests
>>>>>>> 13c222f2b2276f98271901b0b51dc84a34e5dc09
