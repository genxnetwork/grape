#!/bin/bash

export BITBUCKET=$BITBUCKET_COMMIT

git clone git@bitbucket.org:genxglobal/genx-relatives-snakemake.git /media/ci/genx-relatives-snakemake-${BITBUCKET_COMMIT}
cd /media/ci/genx-relatives-snakemake-${BITBUCKET_COMMIT}/
docker build -t genx-relatives-snakemake-${BITBUCKET_COMMIT}:latest -f containers/snakemake/Dockerfile .
envsubst < /media/ag3r/funnel/bin/simulation-dry-run-template.json > /media/ag3r/funnel/bin/simulation-${BITBUCKET_COMMIT}.json
/media/ag3r/funnel/bin/funnel task create /media/ag3r/funnel/bin/simulation-${BITBUCKET_COMMIT}.json