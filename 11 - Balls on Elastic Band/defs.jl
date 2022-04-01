module Defs

export Ball, Elastic, Environment, System, @ep

struct Ball
    # Parameters
    m::Float64  # Ball mass
    μ::Float64  # Ball friction coefficient
    R::Float64  # Ball radius
    I::Float64  # Ball moment of inertia
    # Cd::Float64 # Ball drag coefficient

    # State Vector
    X::Array{Any} # (x, y, ̇x, ̇y, ϕ, ̇ϕ)
end

struct Elastic
    # Parameters
    m::Float64            # Elastic mass
    L₀::Float64           # Elastic natural length
    dynamic::Function     # Elastic dynamic funtion
    # dissipation::Function # Elastic dissipation function

    # State Vector
    X::Array{Any} # (d, ḋ, φ, ̇φ)
end

struct Environment
    # Parameters
    g::Float64 # Gravity
    # ρ::Float64 # Air density
end

struct System
    # Parameters
    μ::Float64 # System reduced mass

    # State Vector
    X::Array{Any} # (x_CM, y_CM, ̇x_CM, ̇y_CM, θ, ̇θ)
end

# Empty placeholder macro
macro ep(n)
    return :(Array{Any}(undef, $n, 100000))
end

end
