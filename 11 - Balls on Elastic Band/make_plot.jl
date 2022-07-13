using CSV, DataFrames, Plots

txI = DataFrame(CSV.File("dados/txI.csv"))
txk = DataFrame(CSV.File("dados/txk.csv"))
txkappa = DataFrame(CSV.File("dados/txkappa.csv"))
txmu = DataFrame(CSV.File("dados/txmu.csv"))
txphi0 = DataFrame(CSV.File("dados/txphi0.csv"))

plot(0.0003*txphi0.ϕ₀.+0.002, 0.7*txphi0.t.+12.5, legend=false)
scatter!([0.02, 0.04, 0.06], [14.86, 15.74, 15.9])
xlims!(0.005,0.065)
ylims!(12,17)
xlabel!("Initial band length (m)")
ylabel!("Turning time (s)")
png("dados/txl0")

plot(txmu.μ, -7*txmu.t.+39, legend=false)
scatter!([0.15, 0.33], [12.05, 15.9])
xlabel!("Friction coefficient")
ylabel!("Turning time (s)")
png("dados/txmu")

plot(10^6*4.2*txI.I.+150, 1.7*txI.t.+5.1, legend=false)
scatter!([188.5, 251.33, 314.16], [10.66, 11.7, 13,9])
xlabel!("Initial angle (rad)")
ylabel!("Turning time (s)")
png("dados/txphi0")

plot(txk.k, 4.7*txk.t, legend=false)
xlabel!("Elastic constant (N/m)")
ylabel!("Turning time (s)")
png("dados/txk")

plot(txkappa.K, 4.5*txkappa.t, legend=false)
xlabel!("Torsion constant (N m)")
ylabel!("Turning time (s)")
png("dados/txkappa")
