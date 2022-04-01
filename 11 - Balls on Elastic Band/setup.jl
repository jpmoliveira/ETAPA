module Setup
using ..Defs
using DataFrames, CSV

export build_System

function make_Ball(ball_mass::Float64, ball_mu::Float64, ball_radius::Float64)
    return Ball(ball_mass, ball_mu, ball_radius, (2 / 5) * ball_mass * (ball_radius ^ 2), @ep(6))
end

function make_Elastic(band_mass::String, band_natural_length::Float64, band_type::String)

    input_file = "./Elastics" * band_type * ".csv" # Tipo de arquivo de dados a ser utilizado

    band_data = DataFrame(CSV.File("./Elastics/$(band_type).csv"))

    band_dynamic = 

    return Elastic(band_mass, band_natural_length, band_dynamic)
end

function build_System(ball_1_type::String, ball_1_radius::Float64,
                      ball_2_type::String, ball_2_radius::Float64,
                      band_type::String, band_length::Float64,
                      air_density::Float64, gravity::Float64)

    ball_1 = make_Ball(ball_1_type, ball_1_radius)
    ball_2 = make_Ball(ball_2_type, ball_2_radius)
    band = make_Elastic(band_type, band_length)
    env = Environment(air_density, gravity)
    
    return System(ball_1, ball_2, band, env, @ep(1), @ep(2))
end

end
