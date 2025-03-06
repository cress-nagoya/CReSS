#!/bin/bash

#PJM -L node=16
#PJM -L rscgrp=small
#PJM -L elapse=24:00:00
#PJM -g <group name>
#PJM -S
#PJM --rsc-list "retention_state=0"
#PJM --mpi "max-proc-per-node=1"
#PJM -j

# Set paths
crsdir="/path/to/output/dir"
conf_path="/path/to/file"
unite_path="/path/to/file"

# Link check file
pushd "${crsdir}" > /dev/null
for file in 0/*.check.txt; do
    if [ -e "$file" ]; then
        base_name=$(basename "$file")
        ln -sf "$file" .
    fi
done
popd > /dev/null

# Run unite
mpiexec -np ${PJM_NODE} -stdin "${conf_path}" "${unite_path}" -Wl,-T