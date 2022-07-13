include("commons.jl")

function eq_motion(u, p, t)
    rx, ry, θ, l, ϕ₁, ϕ₂, rx_dot, ry_dot, θ_dot, l_dot, ϕ₁_dot, ϕ₂_dot = u
    m₁, m₂, l₀, μ, μ_m1, μ_m2, M, R₁, R₂, ΔR, ΣR, I₁, I₂, μ₁, μ₂, k, K, g = p

    drx = rx_dot; dry = ry_dot
    dθ = θ_dot; dl = l_dot
    dϕ₁ = ϕ₁_dot; dϕ₂ = ϕ₂_dot

    d = sqrt(l ^ 2 + ΔR ^ 2)
    fel = Fel(d - ΣR, l₀, ϕ₁ + ϕ₂, k)
    adjusted_fel = fel * l / d
    tel = -Tel(d - ΣR, ϕ₁ + ϕ₂, K)

    N₁ = m₁ * g - fel * ΔR / d
    ∂γ₁ = μ_m1 * θ_dot - R₁ * ϕ₁_dot / d
    vd₁ = sqrt((μ_m1 * l_dot) ^ 2 + (∂γ₁ * l) ^ 2)
    if vd₁ != 0
        Fat1_mod = μ₁ * N₁ / vd₁
        Fat1_x = -Fat1_mod * (μ_m1 * l_dot * cos(θ) - ∂γ₁ * l * sin(θ))
        Fat1_y = -Fat1_mod * (μ_m1 * l_dot * sin(θ) + ∂γ₁ * l * cos(θ))
    else
        Fat1_x = Fat1_y = 0
    end

    N₂ = m₂ * g + fel * ΔR / d
    ∂γ₂ = μ_m2 * θ_dot - R₂ * ϕ₂_dot / d
    vd₂ = sqrt((μ_m2 * l_dot) ^ 2 + (∂γ₂ * l) ^ 2)
    if vd₂ != 0
        Fat2_mod = μ₂ * N₂ / vd₂
        Fat2_x = -Fat2_mod * (-μ_m2 * l_dot * cos(θ) + ∂γ₂ * l * sin(θ))
        Fat2_y = -Fat2_mod * (-μ_m2 * l_dot * sin(θ) - ∂γ₂ * l * cos(θ))
    else
        Fat2_x = Fat2_y = 0
    end

    drx_dot = (Fat1_x + Fat2_x) / M; dry_dot = (Fat1_y + Fat2_y) / M
    dθ_dot = (((m₁ * drx_dot - Fat1_x) * sin(θ) + (Fat1_y - m₁ * dry_dot) * cos(θ)) / μ- 2 * l_dot * θ_dot) / l
    dl_dot = ((Fat1_x - m₁ * drx_dot) * cos(θ) + (Fat1_y - m₁ * dry_dot) * sin(θ) - adjusted_fel) / μ + l * θ_dot ^ 2
    dϕ₁_dot = ((m₁ * (l_dot * cos(θ) - rx_dot) - μ_m2 * m₁ * l * θ_dot * sin(θ)) * ry_dot + (m₁ * (ry_dot - l_dot * sin(θ))
              - μ_m2 * m₁ * l * θ_dot * cos(θ)) * rx_dot - (m₁ * rx + μ * l * cos(θ)) * dry_dot + (m₁ * ry + μ * l * sin(θ))
              * drx_dot + ((dl_dot - l * θ_dot ^ 2) * (ry * cos(θ) - rx * sin(θ)) - (2 * l_dot * θ_dot + l * dθ_dot) *
              (ry * sin(θ) + rx * cos(θ))) * m₁ - 2 * μ_m1 * l * l_dot * θ_dot - μ * μ_m1 * l ^ 2 * dθ_dot) * d / (I₁ * ΔR)
    dϕ₂_dot = I₁ / I₂ * dϕ₁_dot - μ * ΔR / I₂ * d * dθ_dot + (ϕ₂_dot - I₁ * ϕ₁_dot / I₂) * l * l_dot / d ^ 2 + d / (I₂ * ΔR) * (2 * μ * l * l_dot * θ_dot + (μ *
              d ^ 2 + I₁ + I₂) * dθ_dot + Fat1_x * (ry + μ_m1 * l * cos(θ)) - Fat1_y * (rx - μ_m1 * l * cos(θ)) + Fat2_x * (ry + μ_m2 * l * sin(θ)) - Fat2_y *
              (rx + μ_m2 * l * cos(θ)))

    SA[drx, dry, dθ, dl, dϕ₁, dϕ₂, drx_dot, dry_dot, dθ_dot, dl_dot, dϕ₁_dot, dϕ₂_dot]
end

function build_p_normal(m1::Real, m2::Real, l₀::Real, R1::Real, R2::Real, μ₁::Real, μ₂::Real, k::Real, K::Real)
    return [m1, m2, l₀, m1 * m2 / (m1 + m2), m2 / (m1 + m2), m1 / (m1 + m2), m1 + m2,
           R1, R2, R2 - R1, R1 + R2, (2 / 5) * m1 * R1 ^ 2, (2 / 5) * m2 * R2 ^ 2,
           μ₁, μ₂, k, K, 9.8]
end

function build_p_special(m1::Real, m2::Real, l₀::Real, R1::Real, R2::Real, I₁::Real,
                         I₂::Real, μ₁::Real, μ₂::Real, k::Real, K::Real, g::Real)
    return [m1, m2, l₀, m1 * m2 / (m1 + m2), m2 / (m1 + m2), m1 / (m1 + m2), m1 + m2,
            R1, R2, R2 - R1, R1 + R2, I₁, I₂, μ₁, μ₂, k, K, g]
end

function build_p(p::Vector)
    if length(p) == 9
        return build_p_normal(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9])
    elseif length(p) == 12
        return build_p_special(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12])
    end
end

function terminate_affect!(integrator)
    terminate!(integrator)
end

function terminate(u, t, integrator)
    (u[9] < -0.1 && t > 0.1)
end

terminate_cb= DiscreteCallback(terminate, terminate_affect!)

function clean_data(data::DataFrame, p::Vector)
    _, _, _, μ_m1, μ_m2, _, R₁, R₂, ΔR, _, _, _, _, _, _ = p
    rename!(data, [:timestamp, :value1, :value2, :value3, :value4, :value5, :value6,
                   :value7, :value8, :value9, :value10, :value11, :value12] .=> [:t,
                   :rx, :ry, :θ, :l, :ϕ₁, :ϕ₂, :rx_dot, :ry_dot, :θ_dot, :l_dot, :ϕ₁_dot,
                   :ϕ₂_dot])
    data.∂γ₁ = μ_m1 * data.θ_dot - R₁ * data.ϕ₁_dot ./ sqrt.(data.l .^ 2 .+ ΔR ^2)
    data.∂γ₂ = μ_m2 * data.θ_dot - R₂ * data.ϕ₂_dot ./ sqrt.(data.l .^ 2 .+ ΔR ^2)
    data.rx1star = μ_m1 * data.l .* cos.(data.θ)
    data.rx1 = data.rx + μ_m1 * data.rx1star
    data.ry1star = μ_m1 * data.l .* sin.(data.θ)
    data.ry1 = data.ry + data.ry1star
    data.rx2star = -μ_m2 * data.l .* cos.(data.θ)
    data.rx2 = data.rx + data.rx2star
    data.ry2star = -μ_m2 * data.l .* sin.(data.θ)
    data.ry2 = data.ry + data.ry2star
    return data
end

function solution(p::Vector, u0::Vector)
    tspan = (0., 100.)
    ode_problem = ODEProblem{Any}
    p_built = build_p(p)
    ode_problem = ODEProblem(eq_motion, u0, tspan, p_built, callback = terminate_cb)
    sol = solve(ode_problem)
    data = DataFrame(sol)
    return clean_data(data, p)
end

function full_solution(p::Vector, u0::Vector)
    tspan = (0., 100.)
    ode_problem = ODEProblem(eq_motion_dif, u0, tspan, p)
    sol = solve(ode_problem, saveat=0.:0.01:100.)
    data = DataFrame(sol)
    return clean_data(data, p)
end
