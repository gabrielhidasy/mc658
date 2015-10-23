#!/bin/bash

export PATH=$PATH:$HOME/opt/openmpi/bin 
export GUROBI_HOME="/opt/gurobi602/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib"

