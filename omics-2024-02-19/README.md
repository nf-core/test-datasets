# Amazon Omics

This zipfile contains API models for the AWS CLI and Python (boto3).

## Installation

```
unzip omics-2024-02-19.zip
aws configure add-model --service-name omics --service-model file://omics-2022-11-28.json
cp omics-2022-11-28.paginators.json ~/.aws/models/omics/2022-11-28/paginators-2.json
cp omics-2022-11-28.waiters2.json ~/.aws/models/omics/2022-11-28/waiters-2.json
aws omics help
```
