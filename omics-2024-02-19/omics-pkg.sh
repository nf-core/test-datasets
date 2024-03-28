# Omics Packaging

# S3 URIs

s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5447294294/U0a_CGATGT_L001_R1_004.fastq.gz?5gMOaqcXAL9H9pYJIobTOudSH0x_hC4m


s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5447294294/U0a_CGATGT_L001_R2_004.fastq.gz?R62nn0rPU7BbVhuUYToncMt4VZK1Rl9T

## Copying to create Hash

## Commands

RUST_LOG=debug AWS_PROFILE=sales quilt-rs package s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5447294294 --target "quilt+s3://quilt-demos#package=test/omics"


AWS_PROFILE=sales aws s3 ls s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/850787717197/sequenceStore/1867753048/readSet/5447294294/

2024-03-12 16:46:16  225803169 U0a_CGATGT_L001_R1_004.fastq.gz
2024-03-12 16:46:30  245026681 U0a_CGATGT_L001_R2_004.fastq.gz


AWS_PROFILE=sales aws s3api get-object --bucket 850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias --key 850787717197/sequenceStore/1867753048/readSet/5447294294/U0a_CGATGT_L001_R1_004.fastq.gz U0a_CGATGT_L001_R1_004.fastq.gz 
{
    "AcceptRanges": "bytes",
    "LastModified": "2024-03-12T23:46:16+00:00",
    "ContentLength": 225803169,
    "ETag": "\"8ad6a45f2bbb583ab7f9cdc9eb8ae170-3\"",
    "VersionId": "5gMOaqcXAL9H9pYJIobTOudSH0x_hC4m",
    "ContentType": "binary/octet-stream",
    "ServerSideEncryption": "aws:kms",
    "Metadata": {},
    "SSEKMSKeyId": "arn:aws:kms:us-east-1:559620149354:key/c46b5818-b660-4e2b-ad4a-721b24d60d3b",
    "BucketKeyEnabled": true,
    "TagCount": 4
}

AWS_PROFILE=sales aws s3api get-object --bucket 850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias --key 850787717197/sequenceStore/1867753048/readSet/5447294294/U0a_CGATGT_L001_R2_004.fastq.gz U0a_CGATGT_L001_R2_004.fastq.gz

{
    "AcceptRanges": "bytes",
    "LastModified": "2024-03-12T23:46:30+00:00",
    "ContentLength": 245026681,
    "ETag": "\"60692c6fa187e42f2c65b75390254847-3\"",
    "VersionId": "R62nn0rPU7BbVhuUYToncMt4VZK1Rl9T",
    "ContentType": "binary/octet-stream",
    "ServerSideEncryption": "aws:kms",
    "Metadata": {},
    "SSEKMSKeyId": "arn:aws:kms:us-east-1:559620149354:key/c46b5818-b660-4e2b-ad4a-721b24d60d3b",
    "BucketKeyEnabled": true,
    "TagCount": 4
}


