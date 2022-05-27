using Plots, XLSX
include("coisa.jl")

# p_same = build_p(0.06, 0.02549, 0.3, 7.5, 0.00001)
# u0_same = [0., .1, 50., 0., 0., 0.]

# μ_range = 0.1:0.01:0.5
# I_range = 0.000001:0.000001:0.00004
# k_range = 5.:100.
# K_range = 0.00001:0.000001:0.0001

# ϕ₀_range = 20.:2.:200.

# p_μ = build_p.([(0.06, 0.02549, μ, 7.5, 0.00001) for μ in μ_range])
# p_I = build_p_special.([(0.06, 0.02549, I, 0.3, 7.5, 0.00001, 9.8) for I in I_range])
# p_k = build_p.([(0.06, 0.02549, 0.3, k, 0.00001) for k in k_range])
# p_K = build_p.([(0.06, 0.02549, 0.3, 7.5, K) for K in K_range])

# u0_ϕ₀ = [[0., .1, ϕ₀, 0., 0., 0.] for ϕ₀ in ϕ₀_range]

# data_μ = solution.(p_μ, [u0_same for _ in 1:length(p_μ)])
# times_μ = inv_time.(data_μ)
# df_μ = DataFrame("μ" => μ_range, "t" => times_μ)

# data_I = solution.(p_I, [u0_same for _ in 1:length(p_I)])
# times_I = inv_time.(data_I)
# df_I = DataFrame("I" => I_range, "t" => times_I)

# data_k = solution.(p_k, [u0_same for _ in 1:length(p_k)])
# times_k = inv_time.(data_k)
# df_k = DataFrame("k" => k_range, "t" => times_k)

# data_K = solution.(p_K, [u0_same for _ in 1:length(p_K)])
# times_K = inv_time.(data_K)
# df_K = DataFrame("K" => K_range, "t" => times_K)

# data_ϕ₀ = solution.([p_same for _ in 1:length(u0_ϕ₀)], u0_ϕ₀)
# times_ϕ₀ = inv_time.(data_ϕ₀)
# df_ϕ₀ = DataFrame("ϕ₀" => ϕ₀_range, "t" => times_ϕ₀)

# CSV.write("dados/txmu.csv", df_μ)
# CSV.write("dados/txI.csv", df_I)
# CSV.write("dados/txk.csv", df_k)
# CSV.write("dados/txkappa.csv", df_K)
# CSV.write("dados/txphi0.csv", df_ϕ₀)

data = full_solution(build_p(0.6, 0.02549, 0.3, 7.5, 0.000005),[0., .1, 50., 0., 0., 0.])

# plt1 = plot(μ_range, times_μ, legend = false)
# plt2 = plot(I_range, times_I, legend = false)
# plt3 = plot(k_range, times_k, legend = false)
# plt4 = plot(K_range, times_K, legend = false)
# plt5 = plot(ϕ₀_range, times_ϕ₀, legend = false)

# u0_diff = [0., 0., 0., .1, 50., 50., 0., 0., 0., 0., 0., 0.]

# m₁ = 0.0601; m₂ = m₁ / 2; m₂_range = 0.001:0.001:0.1
# R₁ = 0.02549; R₂ = R₁ / 2; R₂_range = 0.005:0.001:0.1

# p_m = build_p.([(m₁, i, R₁, R₂, 0.3, 0.3, 7.5, 0.00001) for i in m₂_range])
# p_R = build_p.([(m₁, m₂, R₁, i, 0.3, 0.3, 7.5, 0.00001) for i in R₂_range])

# data_m = solution.(p_m, [u0_diff for _ in 1:length(p_m)])
# times_m = inv_time.(data_m)

# data_R = solution.(p_R, [u0_diff for _ in 1:length(p_R)])
# times_R = inv_time.(data_R)

# plt6 = plot(m₂_range / m₁, times_m)
# plt7 = plot(R₂_range / R₁, times_R)
