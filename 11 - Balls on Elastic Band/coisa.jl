using DifferentialEquations, StaticArrays, LinearAlgebra, BenchmarkTools, Plots

function Fel(l, phi)
    k = 1
    return k * (l - 3)
end

function Tel(l, phi)
    K = 1
    return K * phi
end

function eq_motion_same(u, p, t)
    _, _, θ, l, ϕ, rx_dot, ry_dot, θ_dot, l_dot, ϕ_dot = u
    m, R, I, μ, g = p

    drx = rx_dot; dry = ry_dot
    dθ = θ_dot; dl = l_dot; dϕ = ϕ_dot

    d = sqrt(l ^ 2)
    fel = Fel(d - 2 * R, 2 * ϕ)
    adjusted_fel = -fel * l / d
    tel = Tel(d - ΣR, 2 * ϕ)

    N = m * g
    ∂γ = θ_dot / 2 - R * ϕ_dot / d
    vd = sqrt((l_dot / 2) ^ 2 + (∂γ * l) ^ 2)
    Fat_mod = μ * N / vd
    Fat_x = -Fat_mod * (l_dot / 2 * cos(θ) - ∂γ * l * sin(θ))
    Fat_y = -Fat_mod * (l_dot / 2 * sin(θ) - ∂γ * l * cos(θ))

    drx_dot = 0; dry_dot = 0
    dθ_dot = ((adjusted_fel * sin(θ) + Fat_y) * cos(θ) - (adjusted_fel * cos(θ) + Fat_x) * sin(θ) - 2 * l_dot * θ_dot) / l
    dl_dot = 2 * ((adjusted_fel * cos(θ) + Fat_x) * cos(θ) + (adjusted_fel * sin(θ) + Fat_y) * sin(θ) + l * dθ ^ 2)
    dϕ_dot = ϕ_dot * (l_dot * (-l_dot/ (d^2 * l)) + θ_dot * tan(θ)) + (tel - fel * R * tan(θ) + μ * d ^ 2 * (
              (dl_dot - l * θ_dot ^ 2) * 2 * R / m * tan(θ) -  (l * dθ_dot + 2 * l_dot * θ_dot) * 2 * R / m * tan(θ)) / l) / I

    SA[drx, dry, dθ, dl, dϕ, drx_dot, dry_dot, dθ_dot, dl_dot, dϕ_dot]
end

function eq_motion_dif(u, p, t)
    rx, ry, θ, l, ϕ₁, ϕ₂, rx_dot, ry_dot, θ_dot, l_dot, ϕ₁_dot, ϕ₂_dot = u
    m₁, m₂, μ, μ_m1, μ_m2, M, R₁, R₂, ΔR, ΣR, I₁, I₂, μ₁, μ₂, g = p

    drx = rx_dot; dry = ry_dot
    dθ = θ_dot; dl = l_dot
    dϕ₁ = ϕ₁_dot; dϕ₂ = ϕ₂_dot

    d = sqrt(l ^ 2 + ΔR ^ 2)
    fel = Fel(d - ΣR, ϕ₁ + ϕ₂)
    adjusted_fel = fel * l / d
    tel = -Tel(d - ΣR, ϕ₁ + ϕ₂)

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
    dϕ₁_dot = ϕ₁_dot * (l_dot * (-(ΔR^2 * l_dot)/ (d^2 * l)) + θ_dot * tan(θ)) + (tel - fel * R₁ * tan(θ) - (m₁ * d * R₁ * drx_dot / cos(θ) - μ * d ^ 2 * (
              (dl_dot - l * θ_dot ^ 2) * (R₁ * (tan(θ) / m₁ + 1 / m₂) + R₂ / m₁ * (tan(θ) - 1)) -  (l * dθ_dot + 2 * l_dot * θ_dot) * (R₁ * (tan(θ) / m₂ - 1 / m₁) +
              R₂ / m₁ * (tan(θ) + 1)))) / l) / I₁
    dϕ₂_dot = I₁ / I₂ * dϕ₁_dot - μ * ΔR / I₂ * d * dθ_dot + (ϕ₂_dot - I₁ * ϕ₁_dot / I₂) * l * l_dot / d ^ 2 + d / (I₂ * ΔR) * (2 * μ * l * l_dot * θ_dot + (μ *
              d ^ 2 + I₁ + I₂) * dθ_dot + Fat1_x * (ry + μ_m1 * l * cos(θ)) - Fat1_y * (rx - μ_m1 * l * cos(θ)) + Fat2_x * (ry + μ_m2 * l * sin(θ)) - Fat2_y *
              (rx + μ_m2 * l * cos(θ)))
    SA[drx, dry, dθ, dl, dϕ₁, dϕ₂, drx_dot, dry_dot, dθ_dot, dl_dot, dϕ₁_dot, dϕ₂_dot]
end

function build_p(m1::Float64, m2::Float64, R1::Float64, R2::Float64, μ₁::Float64, μ₂::Float64)
    return [m1, m2, m1 * m2 / (m1 + m2), m2 / (m1 + m2), m1 / (m1 + m2), m1 + m2,
           R1, R2, R2 - R1, R1 + R2, (2 / 5) * m1 * R1 ^ 2, (2 / 5) * m2 * R2 ^ 2,
           μ₁, μ₂, 9.8]
end

function build_p(m::Float64, R::Float64, μ::Float64)
    return [m, R, (2 / 5) * m * R ^ 2, μ, 9.8]
end

p = build_p(1., 2., 1., 2., 1., 1.)
u0 = [0., 0., 0., 4., 100., 100., 0, 0, 0, 0, 0, 0]
tspan = (0., 100.)
ode_problem = ODEProblem{Any}
if length(p) == 15
    ode_problem = ODEProblem(eq_motion_dif, u0, tspan, p)
elseif length(p) == 5
    ode_problem = ODEProblem(eq_motion_same, u0, tspan, p)
end

sol = solve(ode_problem)
plot(sol,tspan=(0.,0.5))