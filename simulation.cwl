#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE simulation"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE simulation workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
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
inputs:
    []
    # command:
    #     type: string
    #     inputBinding:
    #         position: 1
outputs:
  dry_run:
    type: stdout
baseCommand:
    ["/src/repo/launcher.py", "simulate", "--flow", "ibis", "--conda-prefix", "/tmp", "--configfile", "workflows/pedsim/config.yaml", "--assembly", "hg37", "--ref-directory", "media", "--real-run", "--use-singularity"]