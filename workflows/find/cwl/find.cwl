#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE find"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE find workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: grape:cwl
  NetworkAccess:
    networkAccess: true
  EnvVarRequirement:
    envDef:
      SINGULARITY_TMPDIR: /tmp
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference)
      - $(inputs.workdir)
inputs:
    []
outputs:
  dry_run:
    type: stdout
  results:
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
        "--directory", "workdir",
        "--singularity-args", "-W tmp",
        "--real-run",
        "--flow", "ibis"
    ]
