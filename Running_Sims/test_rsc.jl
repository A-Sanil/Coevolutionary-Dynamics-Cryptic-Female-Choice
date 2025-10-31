# Test script for RSC implementation
using Distributed
if !haskey(ENV, "NO_AUTO_ADDPROCS")
  addprocs(23)
end

using Dates

# Load packages across all cores
@everywhere using Random, Distributions, StatsBase, GLM, DataFrames, CSV, SharedArrays

# Set environment variable
ENV["RUN_SINGLE_TEST"] = "1"

# Include the main model
include("RunModel.jl")
