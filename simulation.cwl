#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE simulation"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE simulation workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  - class: DockerRequirement
    dockerPull: genx_relatives:cwl
  - class: InitialWorkDirRequirement
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
    ["/src/repo/launcher.py", "simulate", "--flow", "ibis", "--conda-prefix", "/tmp", "--configfile", "workflows/pedsim/config.yaml", "--assembly", "hg37", "--ref-directory", "media"]