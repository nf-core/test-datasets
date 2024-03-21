#!/usr/bin/env sh -x

FILE=HASHME.md
SCALE_BUCKET=data-yaml-spec-tests
SCALE_URI=s3://$SCALE_BUCKET
ROOT=scale
HASH=SHA256

# set N to argv 1 or 0
# set PREFIX to argv 2 or 'test/e0'

N=${1:-0}
PREFIX=${2:-test/e0}

# IF N is 0, put file in bucket with that PREFIX
# else call call with N-1 and PREFIX/e$N-1

if [ $N -eq 0 ]; then
    KEY=$ROOT/$PREFIX.txt
    aws s3api put-object --key $KEY --bucket $SCALE_BUCKET --body $FILE --checksum-algorithm $HASH
else
  ./make32.sh $(($N-1)) $PREFIX/e$(($N-1))
  KEY=$PREFIX/e$N/$FILE
fi

