#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE simulation"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE simulation workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: genxnetwork/grape:v1.4
  NetworkAccess:
    networkAccess: true
  EnvVarRequirement:
    envDef:
      SINGULARITY_TMPDIR: /tmp
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference)
inputs:
    []
outputs:
  simulation_result:
    type: Directory
      outputBinding:
        glob: "results"
baseCommand:
    [
        "/src/repo/launcher.py",
        "simulate",
        "--flow", "ibis",
        "--conda-prefix", "/tmp",
        "--configfile", "/src/repo/config.yaml",
        "--assembly", "hg37",
        "--ref-directory", "ref",
        "--real-run"
    ]