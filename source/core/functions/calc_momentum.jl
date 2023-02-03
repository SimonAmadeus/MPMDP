function calc_totalmomentum(particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity})
    P_tot_x::Float64 = 0
    P_tot_y::Float64 = 0
    P_tot_z::Float64 = 0
    P_tot::Float64 = 0
    for i = 1:length(velocity_vec)
        P_tot_x += velocity_vec[i].velocity[1] * particle_vec[i].mass
        P_tot_y += velocity_vec[i].velocity[2] * particle_vec[i].mass
        P_tot_z += velocity_vec[i].velocity[3] * particle_vec[i].mass
    end
    P_tot = P_tot_x + P_tot_y + P_tot_z
    return P_tot
end
