#/usr/bin/env bash

SRC=s3://850787717197-1867753-fepwgrx9iujr5b9pkjudkhpgxwbuhuse1b-s3alias/
DEST=s3://quilt-sales-raw//
KEYS=omics-keys.txt

# For each KEY in omics-keys.txt, copy SRC/KEY to DEST

while read -r KEY; do
  echo "Copying $KEY"
  AWS_PROFILE=sales aws s3 cp $SRC$KEY $DEST
done < $KEYS

