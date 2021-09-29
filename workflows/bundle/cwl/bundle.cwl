#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE bundle downloading"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE reference downloading. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: genxnetwork/grape:v1.4
  NetworkAccess:
    networkAccess: true
  EnvVarRequirement:
    envDef:
      SINGULARITY_TMPDIR: /tmp
inputs:
    []
outputs:
  reference:
    type: Directory
      outputBinding:
        glob: "ref"
baseCommand:
    [
        "/src/repo/launcher.py",
        "bundle",
        "--assembly", "hg37",
        "--conda-prefix", "/tmp",
        "--configfile", "/src/repo/config.yaml",
        "--ref-directory", "ref",
        "--phase",
        "--real-run"
    ]
