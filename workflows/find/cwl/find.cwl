#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE find"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE find workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: genxnetwork/grape:v1.4
  NetworkAccess:
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference)
      - $(inputs.data_bed)
      - $(inputs.data_bim)
      - $(inputs.data_fam)
inputs:
    []
outputs:
  dry_run:
    type: stdout
  result:
    type: File
    outputBinding:
        glob: "results/relatives.tsv"
baseCommand:
    [
        "/src/repo/launcher.py",
        "find",
        "--assembly", "hg37",
        "--conda-prefix", "/tmp",
        "--configfile", "/src/repo/config.yaml",
        "--ref-directory", "media",
        "--real-run",
        "--flow", "ibis"
    ]
