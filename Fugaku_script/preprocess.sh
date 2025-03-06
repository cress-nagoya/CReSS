#!/bin/bash

#PJM -L node=1
#PJM -L rscgrp=small
#PJM -L elapse=24:00:00
#PJM -g <group name>
#PJM -S
#PJM --rsc-list "retention_state=0"
#PJM -j

# Set paths
crsdir="/path/to/output/dir"
bin_dir="/path/to/dir"
conf_path="/path/to/file"

# Make directory to output
mkdir -p "${crsdir}"

# Run pre-processes
"${bin_dir}/terrain.exe" -Wl,-T < "${conf_path}"
"${bin_dir}/surface.exe" -Wl,-T < "${conf_path}"
"${bin_dir}/gridata.exe" -Wl,-T < "${conf_path}"