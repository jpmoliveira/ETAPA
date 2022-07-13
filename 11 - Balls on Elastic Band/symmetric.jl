include("commons.jl")

function eq_motion(u, p, t)
    θ, l, ϕ, θ_dot, l_dot, ϕ_dot = u
    m, l₀, R, I, μ, k, K, g = p

    dθ = θ_dot; dl = l_dot; dϕ = ϕ_dot

    fel = Fel(l - 2 * R, l₀, 2 * ϕ, k)
    tel = -Tel(l - 2 * R, 2 * ϕ, K)

    N = m * g
    ∂γ = θ_dot / 2 - R * ϕ_dot / l
    vd = sqrt((l_dot / 2) ^ 2 + (∂γ * l) ^ 2)
    if vd != 0
        Fat_mod = μ * N / vd
        Fat_x = -Fat_mod * (l_dot / 2 * cos(θ) - ∂γ * l * sin(θ))
        Fat_y = -Fat_mod * (l_dot / 2 * sin(θ) + ∂γ * l * cos(θ))
    else
        Fat_x = Fat_y = 0
    end

    dθ_dot = (2 * (-Fat_x * sin(θ) + Fat_y * cos(θ)) / m - 2 * l_dot * θ_dot) / l
    dl_dot = 2 * (Fat_x * cos(θ) + Fat_y * sin(θ) - fel) / m + l * θ_dot ^ 2
    dϕ_dot = (tel + R * ∂γ * l) / I

    SA[dθ, dl, dϕ, dθ_dot, dl_dot, dϕ_dot]
end

function build_p_normal(m::Real, l₀::Real, R::Real, μ::Real, k::Real, K::Real)
    return [m, l₀, R, (2 / 5) * m * R ^ 2, μ, k, K, 9.8]
end

function build_p_special(m::Real, l₀::Real, R::Real, I::Real, μ::Real, k::Real, K::Real, g::Real)
    return [m, l₀, R, I, μ, k, K, g]
end

function build_p(p::Vector)
    if length(p) == 6
        return build_p_normal(p[1], p[2], p[3], p[4], p[5], p[6])
    elseif length(p) == 8
        return build_p_special(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])
    end
end

function terminate_affect!(integrator)
    terminate!(integrator)
end

function terminate(u, t, integrator)
    (u[4] > 0.1 && t > 0.1)
end

terminate_cb = DiscreteCallback(terminate, terminate_affect!)

function clean_data(data::DataFrame, p::Vector)
    R = p[3]
    rename!(data, [:timestamp, :value1, :value2, :value3, :value4, :value5, :value6]
                .=> [:t, :θ, :l, :ϕ, :θ_dot, :l_dot, :ϕ_dot])
    data.∂γ = data.θ_dot ./ 2 - R * data.ϕ_dot ./ data.l
    data.rx1 = -data.l .* cos.(data.θ)
    data.ry1 = -data.l .* sin.(data.θ)
    data.rx2 = -data.rx1; data.ry2 = -data.ry1
    return data
end

function solution(p::Vector, u0::Vector)
    tspan = (0., 100.)
    p_built = build_p(p)
    ode_problem = ODEProblem(eq_motion, u0, tspan, p_built, callback = terminate_cb)
    sol = solve(ode_problem)
    data = DataFrame(sol)
    return clean_data(data, p)
end

function full_solution(p::Vector, u0::Vector)
    tspan = (0., 100.)
    p_built = build_p(p)
    ode_problem = ODEProblem(eq_motion_same, u0, tspan, p_built)
    sol = solve(ode_problem, saveat=0.:0.01:100.)
    data = DataFrame(sol)
    return clean_data(data, p)
end
