
module TwoTemperatureProblems
    using Reexport
    using Parameters, Kwonly
    @reexport using DifferentialEquations
    using Statistics, StatsBase, LinearAlgebra
    using Makie

    include("Simulation.jl")
    include("AnalysisTools.jl")
    include("Visualization.jl")

end
