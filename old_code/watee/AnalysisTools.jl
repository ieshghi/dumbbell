##############################Trajectory Analysis###############################
# Output Format: Always NamedTuple
# Cx₁x₁(τ) = ⟨(x₁(t+τ)- ⟨x₁⟩)(x₁(t)- ⟨x₁⟩)⟩ time average
# Cx₁x₁(t, t') = ⟨(x₁(t)- ⟨x₁⟩)(x₁(t')- ⟨x₁⟩)⟩  ensemblem average
# Cx₂x₂(τ) = ⟨(x₂(t+τ)- ⟨x₂⟩)(x₂(t)- ⟨x₂⟩)⟩ time average
# Cx₂x₂(t, t') = ⟨(x₂(t)- ⟨x₂⟩)(x₂(t')- ⟨x₂⟩)⟩  ensemblem average
# Cx₁x₂(τ) = ⟨(x₁(t+τ)- ⟨x₁⟩)(x₂(t)- ⟨x₂⟩)⟩ time average
# Cx₁x₂(t, t') = ⟨(x₁(t)- ⟨x₁⟩)(x₂(t')- ⟨x₂⟩)⟩  ensemblem average
# ⟨(x₁(t)- x₁(0))²⟩ = ⟨(x₁(t)- x₁(0))(x₁(t)- x₁(0))⟩ := Cx₁x₁(τ = 0)

function changevariables(x₁, x₂, parameters::TWOTEMPParameters; Options = :ShuraJF)
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = parameters

    if Options == :ShuraJF
        T = [γ₁T₂/(γ₁T₂ + γ₂T₁) γ₂T₁/(γ₁T₂ + γ₂T₁) ; 1 -1]
    elseif Options == :Michael
        T = [γ₁/(γ₁ + γ₂) γ₂/(γ₁ + γ₂) ; 1 -1]
    else
        T = Options
    end
    R = T[1, 1].*x₁ + T[1, 2].*x₂
    r = T[2, 1].*x₁ + T[2, 2].*x₂
    return (R = R, r = r)
end

function _velocity(t, x₁, x₂) #calculate velicity using finite difference of all trajectorues it does not have to be x1 and x2 coordinates
    dx₁ = diff(x₁, dims = 1)
    dx₂ = diff(x₂, dims = 1)
    dt  = diff(t)
    v₁  = broadcast(/, dx₁, dt)
    v₂  = broadcast(/, dx₂, dt)
    return (t = t[1:size(t,1)-1], v₁ = v₁, v₂ = v₂)
end

function velocity_and_statistic(solution) #The argument of this should be the solution that comes out of ParametricProblemSolution
    vtemp  = _velocity(solution.t, solution.x₁, solution.x₂)
    vstat1 = mean_and_std(vtemp.v₁)
    vstat2 = mean_and_std(vtemp.v₂)
    return (ΔT = solution.ΔT, vstat₁ = vstat1, vstat₂ = vstat2, v₁ = vtemp.v₁, v₂ = vtemp.v₂ )
end


function covariance_time(t, x₁, x₂; demean = true) #timeaverage convariance for each path in the ensemblem
  τ = collect(0:size(t,1)-1)
  Cx₁x₁ = autocov(x₁, τ; demean = demean)
  Cx₂x₂ = autocov(x₂, τ; demean = demean)
  Cx₁x₂ = map(n -> crosscov(x₁[:, n], x₂[:, n], τ; demean = demean), collect(1:size(x₁, 2)))
  Cx₁x₂ = [Cx₁x₂[j][i] for i = 1:size(τ, 1), j = 1:size(x₁, 2)]
  return (mode = :τ, τ = t[τ[2:end]] , Cx₁x₁ = Cx₁x₁, Cx₂x₂ = Cx₂x₂, Cx₁x₂ = Cx₁x₂) #return a named tuples
end

function covariance_ensm(t, x₁, x₂; demean = true) #ensemble average convariance out put is 1 path

  if demean == true
    for i = 1:size(t,1)-1
     x₁[i, :] .= x₁[i, :] .- mean(x₁[i, :])
     x₂[i, :] .= x₂[i, :] .- mean(x₂[i, :])
    end
  end

  Cx₁x₁ = [ mean(x₁[i, :].*x₁[j, :]) for i =  1:length(t), j =  1:length(t) ]
  Cx₂x₂ = [ mean(x₂[i, :].*x₂[j, :]) for i =  1:length(t), j =  1:length(t) ]
  Cx₁x₂ = [ mean(x₁[i, :].*x₂[j, :]) for i =  1:length(t), j =  1:length(t) ]

  return (mode = :ensemble, t = t, Cx₁x₁ = Cx₁x₁, Cx₂x₂ = Cx₂x₂, Cx₁x₂ = Cx₁x₂)
end

#For the above Cx₁x₁[i, j] = ⟨(x₁(t[i])- ⟨x₁⟩)(x₁(t[j])- ⟨x₁⟩)⟩
# Example: Cx₁x₁[:, 1] = ⟨(x₁(t[:])- ⟨x₁⟩)(x₁(0)- ⟨x₁⟩)⟩ := Cx₁x₁(t, 0)


################################################################################

#Alternative methods for calculting velcity only works on individual traj for now fix later.

#=function velocity(t, x₁, x₂, parameters::TWOTEMPParameters) #calculate velocity using exact formula
    @unpack k, γ₁, γ₂, T₁, T₂, α, λ₁, λ, n, r₀, rmax = parameters
    v₁Ramp = Vector(undef, size(t, 1) - 1)
    v₂Ramp = Vector(undef, size(t, 1) - 1)

    for i = 2:size(t,1)
        v₁Ramp[i] = (1/γ₁)*mapreduce(n -> (α * λ^2)/(n*π*(λ-λ₁)*λ₁)*sin((2*n*π*x₁[i])/λ)*(sin((n*π*λ₁)/λ))^2, +, 1:n)
                  + (1/γ₁)*mapreduce(n -> (α * λ^2)/(n*π*(λ-λ₁)*λ₁)*cos((2*n*π*x₁[i])/λ)*(sin((2*n*π*λ₁)/λ)), +, 1:n)
        v₂Ramp[i] = (1/γ₂)*mapreduce(n -> (α * λ^2)/(n*π*(λ-λ₁)*λ₁)*sin((2*n*π*x₂[i])/λ)*(sin((n*π*λ₁)/λ))^2, +, 1:n)
                  + (1/γ₂)*mapreduce(n -> (α * λ^2)/(n*π*(λ-λ₁)*λ₁)*cos((2*n*π*x₂[i])/λ)*(sin((2*n*π*λ₁)/λ)), +, 1:n)
    end

    if rmax == 0
        v₁Harmonics = [k*(-x₁[i] + x₂[i])/γ₁ for i = 2:size(t, 1)]
        v₂Harmonics = [k*(-x₂[i] + x₁[i])/γ₂ for i = 2:size(t, 1)]
        v₁ = v₁Ramp .+ v₁Harmonics
        v₂ = v₂Ramp .+ v₂Harmonics
    else
        r = [abs(x₁[i] - x₂[i]) for i = 2:size(t, 1)]
        f = [k*(r[i] - r₀)/(1 - ((r[i] - r0)^2 / rmax^2)) for i = 2:size(t,1)]
        v₁FENE = [(f[i]/r[i])*(x₁[i] - x₂[i])/γ₁]
        v₂FENE = [(f[i]/r[i])*(x₂[i] - x₁[i])/γ₂]
        v₁ = v₁Ramp .+ v₁FENE
        v₂ = v₂Ramp .+ v₂FENE
    end
    return (t = t[1:size(t,1)-1], v₁ = v₁, v₂ = v₂)
end=#

################################################################################
