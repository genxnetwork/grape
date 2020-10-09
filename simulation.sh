#!/bin/bash

HASH=$1

git clone git@bitbucket.org:genxglobal/genx-relatives-snakemake.git /media/ci/genx-relatives-snakemake-${HASH}
cd /media/ci/genx-relatives-snakemake-${HASH}/
docker build -t genx-relatives-snakemake-${HASH}:latest -f containers/snakemake/Dockerfile .
sed -e "s/\${HASH}/$HASH/" /media/ag3r/funnel/bin/simulation-dry-run-template.json > /media/ag3r/funnel/bin/simulation-${HASH}.json
/media/ag3r/funnel/bin/funnel task create /media/ag3r/funnel/bin/simulation-${HASH}.json

#docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx-relatives-snakemake-${HASH}:latest launcher.py simulate --directory /media/ci/sim-${HASH} --real-run