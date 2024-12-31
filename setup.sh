#!/bin/bash

# Check if the Mamba environment name is provided as an argument, otherwise use default
ENV_NAME=${1:-B0nav-LL}

# Step 1: Create the new Mamba environment from the YAML file
echo "Creating Mamba environment: $ENV_NAME"
mamba env create --name $ENV_NAME --file environment_nocuda.yml

# Step 2: Get the Python executable path from the activated environment
ENV_PATH=$(mamba env list | grep "^$ENV_NAME" | awk '{print $NF}')
PYTHON_PATH=$ENV_PATH/bin/python

echo "Python path in the environment: $PYTHON_PATH"

# Step 3: Open Julia and configure the environment
echo "Configuring Julia environment..."

# Run Julia commands to activate, instantiate, and set the Python environment
julia -e "
using Pkg
Pkg.activate(\".\")
Pkg.instantiate()

# Set the Python environment for PyCall
ENV[\"PYTHON\"] = \"$PYTHON_PATH\"
Pkg.build(\"PyCall\")

println(\"Julia environment configured and PyCall built.\")
"