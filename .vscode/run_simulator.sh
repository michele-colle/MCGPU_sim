#!/bin/bash
cd "${WORKSPACE_FOLDER}/example_simulations" 2>/dev/null || cd "$(dirname "$0")/../example_simulations"
./MC-GPU_v1.5b.x MC-GPU_v1.5b_scattered_phantom_mammo_fast.in | tee output_fast.out
