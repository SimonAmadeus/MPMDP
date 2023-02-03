export Force_Sum

# Determine properties fieldwise, like force per cell, e.g.

struct Force_Sum <: Field_Sum
    a::Float64
    freq::Int64
end

function initialize_force_cells(args::System, a::Float64)
    nx = Int64(div(args.box[1], a))
    ny = Int64(div(args.box[2], a))
    nz = Int64(div(args.box[3], a))

    force_sum_x = zeros(Float64, nx, ny, nz)
    force_sum_y = zeros(Float64, nx, ny, nz)
    force_sum_z = zeros(Float64, nx, ny, nz)
    return force_sum_x, force_sum_y, force_sum_z
end

function sum_force!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3}, a::Float64)
    for i in 1:length(particle_vec)
        ix = Int64(div(particle_vec[i].pos[1] + args.box[1] / 2, a)) + 1
        iy = Int64(div(particle_vec[i].pos[2] + args.box[2] / 2, a)) + 1
        iz = Int64(div(particle_vec[i].pos[3] + args.box[3] / 2, a)) + 1

        force_sum_x[ix, iy, iz] += force_vec[i].force[1]
        force_sum_y[ix, iy, iz] += force_vec[i].force[2]
        force_sum_z[ix, iy, iz] += force_vec[i].force[3] 
    end
end

function clear_sums!(force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3})
    for i in 1:size(force_sum_x, 1)
        for j in 1:size(force_sum_x, 2)
            for k in 1:size(force_sum_x, 3)
                force_sum_x[i,j,k] = 0
                force_sum_y[i,j,k] = 0
                force_sum_z[i,j,k] = 0
            end
        end
    end
end
