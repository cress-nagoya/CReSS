#!/bin/bash

#PJM -L "node=8x8x8:torus:strict"
#PJM -L rscgrp=large
#PJM -L elapse=0:10:00
#PJM -g <group name>
#PJM -S
#PJM --mpi "max-proc-per-node=4"
#PJM --rsc-list "retention_state=0"
#PJM -j
#PJM --llio async-close=on

# Recommended configuration: 4 processes, 12 threads per process
export OMP_NUM_THREADS=12

# Set paths
vcoord_file_path="/path/to/file"
conf_path="/path/to/file"
solver_path="/path/to/file"

# Transfer shared files
llio_transfer "${vcoord_file_path}"
llio_transfer "${conf_path}"
llio_transfer "${solver_path}"

# Solve with 4x512 = 2048 processes
mpiexec -np 2048 --vcoordfile "${vcoord_file_path}" -mca mpi_print_stats 1 -stdin "${conf_path}" "${solver_path}" -Wl,-T