function Fel(l, l₀, phi, k)
    return k * (l - l₀)
end

function Tel(l, phi, K)
    return K * phi
end

function inv_time(data::DataFrame)
    for i in 100:length(data.t)
        t = -(data.t[i] * data.θ_dot[i - 1] - data.t[i - 1] * data.θ_dot[i]) / (data.θ_dot[i] - data.θ_dot[i - 1])
        if t > data.t[i - 1] && t < data.t[i]
            return t
        end
    end
end
