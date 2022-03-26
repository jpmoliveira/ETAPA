module Defs

export Ball, Elastic, System, @ep

struct Ball
    # Postition
    x::Vector{Float64}
    y::Vector{Float64}

    # Velocity
    vx::Vector{Float64}
    vy::Vector{Float64}

    # Rotation
    ϕ::Vector{Float64}

    # Parameters
    μ::Float64  # Friction Coefficient
    R::Float64  # Ball Radius
    m::Float64  # Mass
end

struct Elastic
    # State
    L::Vector{Float64}
    ϕ::Vector{Float64}

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

    # Position
    CM::Vector{Float64}
    θ::Vector{Float64}
end

macro ep()
    return :(Vector{Float64}(undef, 100000))
end

end
