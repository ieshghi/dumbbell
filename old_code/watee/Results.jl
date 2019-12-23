push!(LOAD_PATH, "/Users/wsrinin/Documents/Research/2-Projects/TwoTemperature")
using Revise
import TwoTemperatureProblems: ParametricProblemSolution, velocity_and_statistic, visualize_drift
using DifferentialEquations, Makie


#BASE PARAMETER FOR PURE HARMONICS
parametersH  = TWOTEMPParameters(k = 10^(-18), γ₁ = 2*10^(-14), γ₂ = 2*10^(-14), T₁ = 300, T₂ = 300)
#BASE PARAMETER FOR HARMONICS + RAMP
parametersHR = TWOTEMPParameters(k = 0, γ₁ = 2*10^(-14), γ₂ = 2*10^(-14), T₁ = 300, T₂ = 300, α = 1.0, λ₁ = 0.5, λ = 5, n = 1)
#BASE PARAMETER FOR HARMONICS +
parametersFR = TWOTEMPParameters(k = 10^(-18), γ₁ = 2*10^(-14), γ₂ = 2*10^(-14), T₁ = 300, T₂ = 300, α = 1.0, λ₁ = 0.5, λ = 2, n = 2, r₀ = 0.5, rmax = 2)

################################################################################
# The output of is of the ParametricProblemSolution is a namedtuple of the form#
# (T₁, T₂, ΔT, x₁, x₂, t, parameter, solutions)                                #
#  ΔT, x₁, x₂, t I extract them from parameter and solutions
# The follwoing 4 function is all you need I think the code should be self evident.                #                                                                         #
################################################################################

solutionHarmonicsRamp       = ParametricProblemSolution(10, 0.1; trange = (0.0, 3600.0), α = 10^(-17), λ₁ = 1, λ = 5, n = 3)
visualize_RampE(-10:0.1:10, solutionHarmonicsRamp[1].parameter)
animatation_realspace(solutionHarmonicsRamp[1],  "HarmonicsRamp.gif",  solutionHarmonicsRamp[1].parameter ; skip = 100)
vs = velocity_and_statistic.(solutionHarmonicsRamp)


################################################################################
