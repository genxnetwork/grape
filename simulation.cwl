#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE simulation"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE simulation workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: genxnetwork/grape:cwl
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
  log:
    type: File
    outputBinding:
        glob: "outputs/output.txt"
  accuracy:
    type: File
    outputBinding:
        glob: "results/accuracy.png"
  overlaps:
    type: File
    outputBinding:
        glob: "results/overlaps.tsv"
  relatives:
    type: File
    outputBinding:
        glob: "results/relatives.tsv"
  precision_recall:
    type: File
    outputBinding:
        glob: "results/precision_recall.png"
  confusion_matrix:
    type: File
    outputBinding:
        glob: "results/confusion_matrix.png"
  updated_relatives:
    type: File
    outputBinding:
        glob: "results/updated_relatives.tsv"
  pedigree_plot:
    type: File
    outputBinding:
        glob: "results/pedigree_plot.png"
baseCommand:
    ["/src/repo/launcher.py", "simulate", "--flow", "ibis", "--conda-prefix", "/tmp", "--configfile", "workflows/pedsim/config.yaml", "--assembly", "hg37", "--ref-directory", "/home/sashaaero/ref", "--real-run", "--use-singularity"]