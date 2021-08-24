#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE bundle downloading"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE reference downloading. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
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
      -
inputs:
    []
outputs:
  dry_run:
    type: stdout
  affymetrix_chip:
    type: Directory
    outputBinding:
        glob: "1000genome/affymetrix_chip"
  allele_info:
    type: Directory
    outputBinding:
        glob: "1000genome/allele_info"
  bcf:
    type: Directory
    outputBinding:
        glob: "1000genome/bcf"
  genetic_map_b37:
    type: Directory
    outputBinding:
        glob: "genetic_map_b37"
  genetic_map_GRCh37:
    type: Directory
    outputBinding:
        glob: "genetic_map_GRCh37"
  hg38ToHg19:
    type: File
    outputBinding:
        glob: "hg38ToHg19.over.chain.gz"
  fasta:
    type: File
    outputBinding:
        glob: "human*"
  refined_genetic_map:
    type: Directory
    outputBinding:
        glob: "Refined_genetic_map_b37"
  simmap:
    type: Directory
    outputBinding:
        glob: "refined_mf.simmap"
  tables:
    type: Directory
    outputBinding:
        glob: "tables"
baseCommand:
    [
        "/src/repo/launcher.py",
        "bundle",
        "--assembly", "hg37",
        "--conda-prefix", "/tmp",
        "--configfile", "/src/repo/config.yaml",
        "--ref-directory", "media",
        "--singularity-args", "-W tmp",
        "--real-run"
    ]
