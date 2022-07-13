using SpecialFunctions, Roots, Plots

function f(x::Real)
    return sin(x)/(-1/x^3+1/(x^2*tan(x))+1/(2*x))
end

function B(a::Real, τ::Real, γ::Real, P₀::Real, h::Real, m::Real, g::Real, R::Real, vc::Real, Ec::Real, h₀::Real)
    return (a/τ)^2*γ*(P₀/h)*(m*g-(4/3)*sqrt(R)*(1/(1-(vc/Ec)^2)))*(h₀^2)^(1/3)
end

function I(ρ::Real, v::Real, K::Real, a::Real, ω::Real)
    return (ρ*v/2)*(besselj0(K*a)*ω)^2
end

function I_offcenter(ρ::Real, v::Real, K::Real, a::Real, r::Real, ω::Real)
    return (ρ*v/2)*((besselj0(K*a) - besselj0(K*r))*ω)^2
end

function frequency(a::Real, τ::Real, γ::Real, P₀::Real, h::Real, m::Real, g::Real, R::Real, vc::Real, Ec::Real, h₀::Real)
    Coisa = B(a, τ, γ, P₀, h, m, g, R, vc, Ec, h₀)
    f2(x) = f(x) - Coisa
    x = find_zeros(f2, 0, 600)
    K = x/a
    return 343*K/100000
end

function clean(ft::Vector)
    vec = Vector{Any}(undef, length(ft))
    for i in 1:length(ft)
        vec[i] = 0
        if i != 1
            for j in 1:length(ft[i])
                vec[i] = 0
                if abs(vec[i-1] - ft[i][j]) < abs(vec[i-1] - vec[i])
                    vec[i] = ft[i][j]
                end
            end
        else
            vec[i] = ft[i][length(ft[i])]
        end
    end
    return vec
end

function clean_dif(ft::Vector)
    vec = Vector{Any}(undef, length(ft))
    for i in 1:length(ft)
        for j in 1:length(ft[i])
            vec[i] = ft[i][j]
        end
    end
    return vec
end

a = 0.07; h₀ = 0.25; m = 0.00709; R = 0.006; h = 0.0416
γ = 3/2; P₀ = 1.02 * 10^5; g = 9.87; vc = 0.1; Ec = 1; τ = 2

τs = 0.5:0.01:4
fτ = clean_dif(frequency.(a, τs, γ, P₀, h, m, g, R, vc, Ec, h₀))
plot(τs, 10*fτ.+300, legend=false)
xlabel!("Tension (kg/s^2)")
ylabel!("Frequency (Hz)")
png("ftau")

h₀s = 0.25:0.001:0.4
fh₀ = clean(frequency.(a, τ, γ, P₀, h, m, g, R, vc, Ec, h₀s))
plot(100*h₀s, -10*fh₀.+300, legend=false)
xlabel!("Drop height (cm)")
ylabel!("Frequency (Hz)")
png("fh0")

as = 0.1:.001:0.2
fa = clean_dif(frequency.(as, τ, γ, P₀, h, m, g, R, vc, Ec, h₀))
plot(100*as, 10*fa, legend=false)
xlabel!("Membrane radius (cm)")
ylabel!("Frequency (Hz)")
png("fa")
