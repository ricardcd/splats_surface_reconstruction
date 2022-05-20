#!/bin/bash

# User-tunable variables: set them to your environment
BIN_DIR="" # Path to the binaries of the (the install path if you run make install). Leave it blank if the binaries are already in the path. Otherwise, make sure it ends with a trailing slash '/'.
DEMO_DATA_DIR="./data" # Path to the demo data directory
RESULTS_DIR="./results" # Path where to store the results of the demo

# Auxiliary variables (just to print colors)
RED=$'\033[31m'
GREEN=$'\033[32m'
YELLOW=$'\033[33m'
BLUE=$'\033[34m'
NC=$'\033[0m' # no color


# Auxiliary functions
print_step() {
    printf "\n${GREEN}*** $1 ***\n${NC}"
}

print_command() {
    printf ${BLUE}"Command:\n"
    printf "\t$1\n${NC}"
}

# mesh_point_set
print_step "Meshing the point set directly using mesh_point_set"
mkdir -p "${RESULTS_DIR}/mesh_point_set/"
cmd="${BIN_DIR}mesh_point_set -i ${DEMO_DATA_DIR}/sphere.xyz -o ${RESULTS_DIR}/mesh_point_set/sphere.off"
print_command "$cmd"
eval $cmd

# compute_splats
print_step "Computing the splats structure using compute_splats"
mkdir -p "${RESULTS_DIR}/splats/"
cmd="${BIN_DIR}compute_splats -i ${DEMO_DATA_DIR}/sphere.xyz -o ${RESULTS_DIR}/splats/sphere.splat"
print_command "$cmd"
eval $cmd

# mesh_splats
print_step "Meshing the splats using mesh_splats"
mkdir -p "${RESULTS_DIR}/mesh_splats/"
cmd="${BIN_DIR}mesh_splats -i ${RESULTS_DIR}/splats/sphere.splat -o ${RESULTS_DIR}/mesh_splats/sphere.off"
print_command "$cmd"
eval $cmd

# orient_splats
print_step "Globally orienting the splats using orient_splats"
mkdir -p "${RESULTS_DIR}/splats/"
cmd="${BIN_DIR}orient_splats -i ${RESULTS_DIR}/splats/sphere.splat -o ${RESULTS_DIR}/splats/sphere_oriented.splat"
print_command "$cmd"
eval $cmd

# mesh_oriented_splats
print_step "Meshing the oriented splats using mesh_oriented_splats"
mkdir -p "${RESULTS_DIR}/mesh_oriented_splats/"
cmd="${BIN_DIR}mesh_oriented_splats -i ${RESULTS_DIR}/splats/sphere_oriented.splat -o ${RESULTS_DIR}/mesh_oriented_splats/sphere.off"
print_command "$cmd"
eval $cmd

# mesh_splats_gcut
print_step "Meshing the splats using mesh_splats_gcut"
mkdir -p "${RESULTS_DIR}/mesh_splats_gcut/"
cmd="${BIN_DIR}mesh_splats_gcut -i ${RESULTS_DIR}/splats/sphere.splat -o ${RESULTS_DIR}/mesh_splats_gcut/sphere.off"
print_command "$cmd"
eval $cmd

# mesh_splats_ncut
print_step "Meshing the splats using mesh_splats_ncut"
mkdir -p "${RESULTS_DIR}/mesh_splats_ncut/"
cmd="${BIN_DIR}mesh_splats_ncut -i ${RESULTS_DIR}/splats/sphere.splat -o ${RESULTS_DIR}/mesh_splats_ncut/sphere.off"
print_command "$cmd"
eval $cmd

