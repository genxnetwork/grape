#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE preprocessing"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE preprocess workflow. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
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
      - $(inputs.vcf)
inputs:
    []
outputs:
  dry_run:
    type: stdout
  data_bed:
    type: File
    outputBinding:
        glob: "preprocessed/data.bed"
  data_bim:
    type: File
    outputBinding:
        glob: "preprocessed/data.bim"
  data_fam:
    type: File
    outputBinding:
        glob: "preprocessed/data.fam"
baseCommand:
    [
        "/src/repo/launcher.py",
        "preprocess",
        "--assembly", "hg37",
        "--conda-prefix", "/tmp",
        "--configfile", "/src/repo/config.yaml",
        "--ref-directory", "media",
        "--vcf-file", "vcf_file",
        "--singularity-args", "-W tmp",
        "--real-run"
    ]
