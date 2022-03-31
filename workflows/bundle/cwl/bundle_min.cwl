#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "GRAPE"
label: "GRAPE bundle minimal downloading"
cwlVersion: v1.1
doc: |
    A Docker container for the GRAPE reference downloading. See the [GRAPE](https://github.com/genxnetwork/grape) GitHub repo for more information.
requirements:
  DockerRequirement:
    dockerPull: genxnetwork/grape:v1.4
  NetworkAccess:
    networkAccess: true
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
        "--real-run"
    ]
