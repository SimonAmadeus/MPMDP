function calc_inst_temp(particle_vec::Vector{Particle}, K_tot::Float64)
    dim::Int64 = 3
    T::Float64 = K_tot * 2 / (length(particle_vec) * dim - dim)
    return T
end
