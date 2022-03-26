module Setup
using ..Defs

export build_System

function make_Ball(ball_type::String, ball_radius::Float64)

    input_file = "./Balls" * ball_type
    pars = Dict{Any, Any}("mu"=>undef, "rho"=>undef)
    open(input_file, "r") do input
        for line in eachline(input)
            tokens = split(line, ":")
            tokens = rstrip.(lstrip.(tokens))
            name = string(tokens[1])
            value = string(tokens[2])
            pars[name] = value
        end
    end

    return Ball(@ep, @ep, @ep, @ep, @ep, pars["μ"], ball_radius,
                (4/3) * pars["rho"] * π * ball_radius ^ 3)
end

function make_Elastic(band_type::String, band_length::Float64)
    input_file = "./Elastics" * ball_type
    pars = Dict{Any, Any}("lambda"=>undef, "dynamic"=>undef)
    open(input_file, "r") do input
        for line in eachline(input)
            tokens = split(line, ":")
            tokens = rstrip.(lstrip.(tokens))
            name = string(tokens[1])
            value = string(tokens[2])
            pars[name] = value
        end
    end

    return Elastic(@ep, @ep, band_length,
                   pars["lambda"] * band_length, pars["dynamic"])
end

function build_System(ball_1_type::String, ball_1_radius::Float64,
                      ball_2_type::String, ball_2_radius::Float64,
                      band_type::String, band_length::Float64)

    ball_1 = make_Ball(ball_1_type, ball_1_radius)
    ball_2 = make_Ball(ball_2_type, ball_2_radius)
    band = make_Elastic(band_type, band_length)
    System(ball_1, ball_2, band, @ep, @ep)

end

end
