module Defs

export Ball, Elastic, System, @ep

struct Ball
    # State
    # x, y, vx, vy, ϕ
    X::Array{Any}

    # Parameters
    μ::Float64  # Friction Coefficient
    R::Float64  # Ball Radius
    m::Float64  # Mass
end

struct Elastic
    # State
    # L, ϕ
    X::Array{Any}

    # Parameters
    L₀::Float64 # Initial Length
    m::Float64  # Mass

    # Elastic Forces and Torques
    dynamic::Function
end

struct System
    # Members
    ball_1::Ball
    ball_2::Ball
    elastic::Elastic

    # State
    # CM, θ
    X::Array{Any}
end

macro ep()
    return :(Vector{Float64}(undef, 100000))
end

end
