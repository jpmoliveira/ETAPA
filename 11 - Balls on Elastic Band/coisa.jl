using DifferentialEquations, StaticArrays, LinearAlgebra, BenchmarkTools

function Fel(l, phi)
    k = 1
    return k * l
end

function Tel(l, phi)
    K = 1
    return K * phi
end

function coisa(u, p, t)
    rx, ry, θ, l, ϕ₁, ϕ₂, rx_dot, ry_dot, θ_dot, l_dot, ϕ₁_dot, ϕ₂_dot = u
    m1, m2, μ_m1, μ_m2, M, R₁, R₂, ΔR, ΣR, μ₁, μ₂, g = p

    drx = rx_dot; dry = ry_dot
    dθ = θ_dot; dl = l_dot
    dϕ₁ = ϕ₁_dot; dϕ₂ = ϕ₂_dot

    d = sqrt(l ^ 2 + ΔR ^ 2)
    fel = Fel(d - ΣR, ϕ₁ + ϕ₂)
    adjusted_fel = fel * l / d

    N₁ = m1 * g - fel * ΔR / d
    vd₁ = sqrt((μ_m1 * l_dot) ^ 2 + (((μ_m1) * θ_dot - R₁ * ϕ₁) * l) ^ 2)
    ∂γ₁ = μ_m1 * θ_dot - R₁ * ϕ₁_dot / d
    Fat1_mod = μ₁ * N₁ / vd₁
    Fat1_x = -Fat1_mod * (μ_m1 * l_dot * cos(θ) - ∂γ₁ * l * sin(θ))
    Fat1_y = -Fat1_mod * (μ_m1 * l_dot * sin(θ) - ∂γ₁ * l * cos(θ))

    N₂ = m2 * g + fel * ΔR / d
    vd₂ = sqrt((μ_m2 * l_dot) ^ 2 + (((μ_m2) * θ_dot - R₂ * ϕ₂) * l) ^ 2)
    ∂γ₂ = μ_m2 * θ_dot - R₂ * ϕ₂_dot / d
    Fat2_mod = μ₂ * N₂ / vd₂
    Fat2_x = Fat2_mod * (μ_m2 * l_dot * cos(θ) - ∂γ₂ * l * sin(θ))
    Fat2_y = Fat2_mod * (μ_m2 * l_dot * sin(θ) + ∂γ₂ * l * cos(θ))

    drx_dot = (Fat1_x + Fat2_x) / M; dry_dot = (Fat1_y + Fat2_y) / M
    dθ_dot = (adjusted_fel* sin(θ) + Fat1_y - m1 * drx_dot) * cos(θ) - (adjusted_fel + Fat1_x - )
end

dvec = [-l * cos(theta); -l * sin(theta); R2 - R1]
d = norm([-l * cos(theta); -l * sin(theta); R2 - R1])
Dvec = [-l * cos(theta); -l * sin(theta)]
D = norm([-l * cos(theta); -l * sin(theta)])

Fat = Fat1 + Fat2

rxdotdot = Fat[1] / M
rydotdot = Fat[2] / M
thetadotdot = ((Felvec[2] + Fat1[2] - m1 * rydotdot) * cos(theta) - (Felvec[1] + Fat1[1] - m1 * rxdotdot) * sin(theta) - 2 * ldot * thetadot) / l
ldotdot = ((Felvec[1] + Fat1[1] - m1 * rxdotdot) * cos(theta) +(Felvec[2] + Fat1[2] - m1 * rydotdot) + l * thetadot ^ 2) / mu
phi1dotdot = phi1dot * ((ddot / d) + thetadot * tan(theta) - (ldot / l)) + (Tel(l, phi) - Fel(l, phi) * R1 * tan(theta) - (m1 * R1 * rydotdot / cos(theta) - mu ^ 2 * d * ((ldotdot - l * thetadot ^ 2) *
                (R1 * (tan(theta) / m1 + 1 / m2) + R2 / m2 * (tan(theta - 1))) - (l * thetadotdot + 2 * ldot * thetadot) * (R1 * (tan(theta) / m2 - 1 / m1) + R2 / m1 * (tan(theta) + 1)))) / l) / I1
phi2dotdot = I1 / I2 * phi1dotdot + (Fel * (R2 - R1) * cos(theta) - (m1 * R1 + m2 * R2) * d * rxdotdot) / (I2 * sin(theta)) - (phi2dot + (mu * d * (R2 - R1) * thetadot - I1 * phi1dot) / I2) * (l / ldot
                + thetadot / tan(theta)) - (mu * d * (R2 - R1) * thetadotdot - I1 * phi1dotdot + I1 * ddot * phi1dot / d) / I2 - ddot * phi2dot / d
