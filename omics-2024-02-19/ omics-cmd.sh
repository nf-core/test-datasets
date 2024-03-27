# https://s3.console.aws.amazon.com/s3/buckets/omics-beta-assets

aws configure add-model --service-name omics --service-model file://omics-2022-11-28.json
cp omics-2022-11-28.paginators.json ~/.aws/models/omics/2022-11-28/paginators-2.json
cp omics-2022-11-28.waiters2.json ~/.aws/models/omics/2022-11-28/waiters-2.json


arn:aws:iam::850787717197:role/omics-quilt-omicsquiltomicsserviceroleB2864563-badZYM0haNEV # us-east-1


AWS_PROFILE=sales aws omics list-reference-stores
            "arn": "arn:aws:omics:us-east-1:850787717197:referenceStore/8277373997",
            "id": "8277373997",
            "name": "Reference store",
            "description": "Default reference store ",
            "creationTime": "2024-03-12T23:12:26.219000+00:00"


AWS_PROFILE=sales aws omics start-reference-import-job \
  --role-arn arn:aws:iam::850787717197:role/omics-quilt-omicsquiltomicsserviceroleB2864563-badZYM0haNEV \
  --reference-store-id 8277373997 \
  --sources sourceFile=s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta,name=OmicsBetaReferenceStore

# Reference store
## hg38

2881303957

omics://850787717197.storage.us-east-1.amazonaws.com/8277373997/reference/2881303957

arn:aws:omics:us-east-1:850787717197:referenceStore/8277373997/reference/2881303957

# Sequence Store omics-beta-assets
## 1867753048

arn:aws:omics:us-east-1:850787717197:sequenceStore/1867753048

AWS_PROFILE=sales aws omics get-sequence-store --id 1867753048
AWS_PROFILE=sales aws omics list-read-sets --sequence-store-id 1867753048
AWS_PROFILE=sales aws omics get-read-set-metadata --id 5091960702 --sequence-store-id 1867753048 | jq -r .files.source1.s3Access.s3Uri
AWS_PROFILE=sales aws s3 ls s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5091960702/

AWS_PROFILE=sales aws s3 ls s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/ --recursive

AWS_PROFILE=sales aws s3 ls s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5091960702/HG00405.final.cram

## S3API

AWS_PROFILE=sales aws s3api get-object --bucket 850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias --key 850787717197/sequenceStore/1867753048/readSet/5091960702/HG00405.final.cram.crai HG00405.final.cram.crai

```json
{
    "AcceptRanges": "bytes",
    "LastModified": "2024-03-12T23:46:42+00:00",
    "ContentLength": 1334307,
    "ETag": "\"5c981ab5c452b538fafabb6574e66c82-1\"",
    "VersionId": "kDds52DrVLhFe24U38o6BIxIo8oL9kgS",
    "ContentType": "binary/octet-stream",
    "ServerSideEncryption": "aws:kms",
    "Metadata": {},
    "SSEKMSKeyId": "arn:aws:kms:us-east-1:559620149354:key/c46b5818-b660-4e2b-ad4a-721b24d60d3b",
    "BucketKeyEnabled": true,
    "TagCount": 4
}
```
# list-objects-v2
# list-object-versions

AWS_PROFILE=sales aws s3api list-objects-v2 --bucket 850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias 

AWS_PROFILE=sales aws s3api list-object-versions --bucket 850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias --key 850787717197/sequenceStore/1867753048/readSet/5091960702/HG00405.final.cram.crai