#!/bin/bash

GITHUB_REPO_URL=$1
BRANCH_NAME=$2
COMMIT_HASH=$3

export REF=$3

git clone --single-branch --branch $BRANCH_NAME $GITHUB_REPO_URL /media/ci/genx-relatives-snakemake-$COMMIT_HASH > log_$COMMIT_HASH.txt
cd /media/ci/genx-relatives-snakemake-$COMMIT_HASH/
docker build -t genx-relatives-snakemake-$COMMIT_HASH:latest -f containers/snakemake/Dockerfile . >> log_$COMMIT_HASH.txt
envsubst < /media/ag3r/funnel/bin/simulation-template.json > /media/ag3r/funnel/bin/simulation-$COMMIT_HASH.json
/media/ag3r/funnel/bin/funnel task create /media/ag3r/funnel/bin/simulation-$COMMIT_HASH.json