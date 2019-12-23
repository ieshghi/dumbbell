###########################Problems Components################################
using  DifferentialEquations
import DifferentialEquations:solve

const kB = 1.381*10^(-17) #N⋅μm/K

function HarmonicsE(x1, x2, p) #This function also takes a pair of vector of x₁ and x₂, the unit of energy is Kelvin
  @unpack k = p
  (1/2) .* (k/kB) .* (x1 .- x2).^2
end

function Harmonics(dx, x, p, t)
  @unpack k, γ₁, γ₂ = p
  dx[1] = k*(-x[1] + x[2])/γ₁
  dx[2] = k*(-x[2] + x[1])/γ₂
end

function FENEE(x1, x2, p)     #This function also takes a pair of vector of x₁ and x₂, the unit of energy is Kelvin
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = p
    r = x1 .- x2
    -0.5 .* rmax^2 .* (k/kB) .* log.(1 .- ((r .- r₀) ./ rmax).^2) #The unit of energy is Kelvin
end

function FENEE(r, p)
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = p
    -0.5 .* rmax^2 .* (k/kB) .* log.(1 .- ((r .- r₀) ./ rmax).^2)
end

function FENE(dX, X, p, t)
  @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = p #extract problem parameters from parameter array

  r = abs(X[1]-X[2])
  f = (k/kB)*(r - r₀)/(1 - ((r - r0)^2 / rmax^2))

  dX[1] = (f/r)*(X[1]-X[2])/γ₁
  dX[2] = (f/r)*(X[2]-X[1])/γ₂
end

function RampE(x, p) #apply to single coordinate of the particle or vector of that coordinates
  @unpack α, λ₁, λ, n = p
  α = α/kB
  u1 = mapreduce(n ->  (α * λ^3)/(2*n^2*π^2*(λ-λ₁)*λ₁)*(sin((n*π*λ₁)/λ))^2*cos((2*n*π*x)/λ), +, 1:n)
  u2 = mapreduce(n -> -(α * λ^3)/(4*n^2*π^2*(λ-λ₁)*λ₁)*(sin((2*n*π*λ₁)/λ))*sin((2*n*π*x)/λ), +, 1:n)
  U  = α/2 + (2/λ)*(u1 + u2)
end

function RampF(dX, X, p, t)
  @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = p

  dX[1] = (α/γ₁)*mapreduce(n -> (2*λ)/(n*π*(λ-λ₁)*λ₁)*sin((2*n*π*X[1])/λ)*(sin((n*π*λ₁)/λ))^2, +, 1:n)
        + (α/γ₁)*mapreduce(n -> (2*λ)/(n*π*(λ-λ₁)*λ₁)*cos((2*n*π*X[1])/λ)*(sin((2*n*π*λ₁)/λ)), +, 1:n)

  dX[2] = (α/γ₂)*mapreduce(n -> (2*λ)/(n*π*(λ-λ₁)*λ₁)*sin((2*n*π*X[2])/λ)*(sin((n*π*λ₁)/λ))^2, +, 1:n)
        + (α/γ₂)*mapreduce(n -> (2*λ)/(n*π*(λ-λ₁)*λ₁)*cos((2*n*π*X[2])/λ)*(sin((2*n*π*λ₁)/λ)), +, 1:n)
end

function σ(dX, X, p, t)
  @unpack γ₁, γ₂, T₁, T₂ = p
  dX[1] = sqrt(2*kB*T₁/γ₁)
  dX[2] = sqrt(2*kB*T₂/γ₂)
end

HarmonicsRamp(dX, X, p, t) = RampF(dX, X, p, t) + Harmonics(dX, X, p, t)
FENERamp(dX, X, p, t)      = RampF(dX, X, p, t) + FENE(dX, X, p, t)

###########################Problems################################

#parameter = [k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax]
#The distance is measured in the unit of μm the rest is standard SI units
# γ   ≈ 1.88*10^(-8)  Ns/(m) = 1.88*10^(-14) Ns/(μm)
# kB  ≈ 1.381*10^(-23) J/ K   = 1.381*10^(-17) N⋅μm/K
# k   ≈ 10^(-18)      N/μm
# k/kB ≈ 0.072 K/(μm)²
# k/γ ≈ 5.31*10^(-5)  1/s
# (kB)*T/γ ≈ ( 0.000735 μm²/(s K)) * T
# D = 0.221 μm²/s
# α = we choose value of α by saying that when the spring is stretched over
# one period of ramp it gain enough energy to overcome the barier
# α ≈ kλ² = which in the unit that we expressed energy in term of tempature
# and length in term of micron this has the unit of temperature it is roughly of order
# 1, harder to overcome if α is large.
# α/γ ≈

mutable struct TWOTEMPParameters
  k
  γ₁ ; γ₂
  T₁ ; T₂
  α  ; λ₁ ; λ ; n
  r₀ ; rmax
  @add_kwonly TWOTEMPParameters(k, γ₁, γ₂, T₁, T₂, α = 0.0, λ₁ =0.0, λ = 0.0, n = 0 , r₀ = 0.0, rmax = 0.0) = new(k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax)
end

struct TWOTEMPProblem
    parameters  ::TWOTEMPParameters
    SDE
    Flag
end

function TWOTEMPProblem(x₀, tspan, parameters::TWOTEMPParameters)
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = parameters
    if rmax == 0 && α == 0
        SDE = SDEProblem(Harmonics, σ, x₀, tspan, parameters)
        TWOTEMPProblem(parameters, SDE, "Harmonics Problems")
    elseif rmax == 0 && α != 0
        SDE = SDEProblem(HarmonicsRamp, σ, x₀, tspan, parameters)
        TWOTEMPProblem(parameters, SDE, "Harmonics + Ramp Problems")
    elseif rmax != 0 && α == 0
        SDE = SDEProblem(FENE, σ, x₀, tspan, parameters)
        TWOTEMPProblem(parameters, SDE, "FENE")
    elseif rmax != 0 && α != 0
        SDE = SDEProblem(FENERamp, σ, x₀, tspan, parameters)
        TWOTEMPProblem(parameters, SDE, "FENE + Ramp Problems")
    end
end

function DifferentialEquations.solve(problem::TWOTEMPProblem, N; arg...) #N is the number of realization, example of arg is solve(problem, EM(), dt = 0.01, saveat = 0.2)
    solutions = Vector(undef, N)
    for i = 1:N
        solutions[i] = DifferentialEquations.solve(problem.SDE; arg...)
    end
    x₁   = [solutions[j][1, i] for i = 1:size(solutions[1].t,1), j = 1:N]
    x₂   = [solutions[j][2, i] for i = 1:size(solutions[1].t,1), j = 1:N]
    time =  solutions[N].t #Here we assume that the time for each trajectory is equal
    return (x₁ = x₁, x₂ = x₂, t = vec(time), solutions = solutions) #it returns named tuples
end

function ParametricProblemSolution(ΔT, dt; trange = (0.0, 100.0), T = 300, α = 0, λ₁ = 0, λ = 0, n = 0, r₀ = 0, rmax = 0, k = 10^(-18), γ₁ = 2*10^(-14), γ₂ = 2*10^(-14) )
    δT = collect(-ΔT:dt:ΔT)
    solutions = Vector(undef, size(δT, 1))
    for i = 1:size(δT, 1)
        parameter    = TWOTEMPParameters(k = k, γ₁ = γ₁, γ₂ = γ₂, T₁ = T, T₂ = T + δT[i], α = α, λ₁ = λ₁, λ = λ, n = n, r₀ = r₀, rmax = rmax)
        problem      = TWOTEMPProblem([0.0, 0.0], trange , parameter) #This is the standard SDE problems
        solution     = solve(problem, 10;  alg = EM(), dt = 0.01, saveat = 0.2) #play with this line if you want to change shit!
        solutions[i] = (T₁ = T, T₂ = T + δT[i], ΔT = δT[i], solution..., parameter = parameter)
    end
    return solutions
end
################################################################################
