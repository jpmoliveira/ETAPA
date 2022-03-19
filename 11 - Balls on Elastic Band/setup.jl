module Setup
using ..Defs

export build_System

function make_Ball(ball::String)

end

function make_Elastic(elastic::String)

end

function build_System(ball_1_type::String, ball_2_type::String, band_type::String)
    ball_1 = make_Ball(ball_1_type)
    ball_2 = make_Ball(ball_2_type)
    band = make_Elastic(band_type)
    u = Vector{Float64}(undef, 100000)
    System(ball_1, ball_2, band, u, u)
end

end
