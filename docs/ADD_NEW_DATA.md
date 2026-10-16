# How to add and use new test dataset

Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things requested when adding a new test dataset.

- [ ] Familiarise yourself with the test-data [specifications](https://nf-co.re/docs/specifications/test-data/overview).
- [ ] Check [here](https://github.com/nf-core/test-datasets/branches/all) that there isn't already a branch containing data that could be used.
  - If this is the case, follow the [documentation on how to use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)
- [ ] Fork the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets) to your GitHub account.
- [ ] Create a new branch on your fork.
- [ ] Add your test dataset.
  - [ ] Warning: this repository is extremely large. If you clone locally, ensure you use the command `git clone <url> --branch <branch> --single-branch` to clone only the branch of interest.
- [ ] Make a PR on a new branch with a relevant name.
	- [ ] In case you accidentally commited large datasets in a PR, drop the specific commits using git rebase.
- [ ] Wait for the PR to be merged. Once merged, use your new test data as described on the main [README](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md).
