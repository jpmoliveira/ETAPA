using DifferentialEquations, LinearAlgebra, BenchmarkTools

dvec = [-l * cos(theta); -l * sin(theta); R2 - R1]
d = norm([-l * cos(theta); -l * sin(theta); R2 - R1])
Dvec = [-l * cos(theta); -l * sin(theta)]
D = norm([-l * cos(theta); -l * sin(theta)])

N1 = m1 * g - Fel(l, phi) * (R2 - R1) / d
N2 = m2 * g + Fel(l, phi) * (R2 - R1) / d
vd1 = norm([(mu / m1) * ldot; ((mu / m1) * thetadot - R1 * phi1dot) * l])
Fat1 = -mu1 * N1 / vd1 * [(mu / m1) * ldot * cos(theta) - ((mu / m1) * thetadot - R1 * phi1dot / d) * l * sin(theta);
                          (mu / m1) * ldot * sin(theta) + ((mu / m1) * thetadot - R1 * phi1dot / d) * l * cos(theta)]
vd2 = norm([(mu / m2) * ldot; ((mu / m2) * thetadot - R2 * phi2dot) * l])
Fat2 = -mu2 * N2 / vd2 * [-(mu / m2) * ldot * cos(theta) + ((mu / m2) * thetadot - R2 * phi2dot / d) * l * sin(theta);
                          -(mu / m2) * ldot * sin(theta) - ((mu / m2) * thetadot - R2 * phi2dot / d) * l * cos(theta)]
Fat = Fat1 + Fat2
Felvec = Fel(l, phi1 + phi2) * Dvec

rxdotdot = Fat[1] / M
rydotdot = Fat[2] / M
thetadotdot = ((Felvec[2] + Fat1[2] - m1 * rydotdot) * cos(theta) - (Felvec[1] + Fat1[1] - m1 * rxdotdot) * sin(theta) - 2 * ldot * thetadot) / l
ldotdot = ((Felvec[1] + Fat1[1] - m1 * rxdotdot) * cos(theta) +(Felvec[2] + Fat1[2] - m1 * rydotdot) + l * thetadot ^ 2) / mu
phi1dotdot = 
phi2dotdot = 
