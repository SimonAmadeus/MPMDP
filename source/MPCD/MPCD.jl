export MPCD

struct MPCD <: Collisions
    a::Float64 # Cell size.
    δt::Int64 # Time between collisions.
    θ::Float64 # Angle in degrees.
end

function apply_collisions!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, collisions::MPCD)
    nx = Int64(div(args.box[1], collisions.a))
    ny = Int64(div(args.box[2], collisions.a))
    nz = Int64(div(args.box[3], collisions.a))
    p_counter = zeros(Int64, nx, ny, nz)
    v_cm = zeros(Float64, nx, ny, nz, 3)
    random_axis = Array{Array{Float64, 1}}(undef, nx, ny, nz)
    #random_axis = zeros(Float64, nx, ny, nz, 3)
    step_size = 1e-6
    shift = rand(((- collisions.a / 2) : step_size : (collisions.a / 2)), 3)
    x = similar(particle_vec, Float64) 
    y = similar(particle_vec, Float64)
    z = similar(particle_vec, Float64)
    for i in eachindex(particle_vec)
        x[i] = mod(particle_vec[i].pos[1] + args.box[1] / 2 + shift[1], args.box[1])
        y[i] = mod(particle_vec[i].pos[2] + args.box[2] / 2 + shift[2], args.box[2])
        z[i] = mod(particle_vec[i].pos[3] + args.box[3] / 2 + shift[3], args.box[3])
    end
    for i in eachindex(particle_vec)
        #ix = Int64(div(x[i] + args.box[1] / 2, collisions.a)) + 1
        #iy = Int64(div(y[i] + args.box[2] / 2, collisions.a)) + 1
        #iz = Int64(div(z[i] + args.box[3] / 2, collisions.a)) + 1
        ix = Int64(div(x[i], collisions.a)) + 1
        iy = Int64(div(y[i], collisions.a)) + 1
        iz = Int64(div(z[i], collisions.a)) + 1

        p_counter[ix, iy, iz] += 1

        # Accumulate velocity of particles in cell.
        v_cm[ix, iy, iz, 1] += velocity_vec[i].velocity[1]
        v_cm[ix, iy, iz, 2] += velocity_vec[i].velocity[2]
        v_cm[ix, iy, iz, 3] += velocity_vec[i].velocity[3]
    end

    # Compute center of mass velocity for each cell.
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                if p_counter[i, j, k] > 0
                    v_cm[i, j, k, 1] /= p_counter[i, j, k]
                    v_cm[i, j, k, 2] /= p_counter[i, j, k]
                    v_cm[i, j, k, 3] /= p_counter[i, j, k]
                end
            end
        end
    end

    # Create random axis for each cell.
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                random_axis[i, j, k] = normalize(Float64.(rand(3)))
            end
        end
    end

    # velocity update in collision step
    for i in eachindex(particle_vec)
        ix = Int64(div(x[i], collisions.a)) + 1
        iy = Int64(div(y[i], collisions.a)) + 1
        iz = Int64(div(z[i], collisions.a)) + 1

        # update particle velocity
        v_cm_cell = [v_cm[ix, iy, iz, 1], v_cm[ix, iy, iz, 2], v_cm[ix, iy, iz, 3]]
        v_rel = velocity_vec[i].velocity .- v_cm_cell
        R = generate_rotation_matrix(random_axis[ix, iy, iz], collisions.θ)

        velocity_vec[i].velocity[1] = v_cm_cell[1] + R[1, 1] * v_rel[1] + R[1, 2] * v_rel[2] + R[1, 3] * v_rel[3]
        velocity_vec[i].velocity[2] = v_cm_cell[2] + R[2, 1] * v_rel[1] + R[2, 2] * v_rel[2] + R[2, 3] * v_rel[3]
        velocity_vec[i].velocity[3] = v_cm_cell[3] + R[3, 1] * v_rel[1] + R[3, 2] * v_rel[2] + R[3, 3] * v_rel[3]
    end

    # Thermostat:
    #apply_MBS_thermostat!(args, particle_vec, velocity_vec, v_cm, collisions)

end

function apply_MBS_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, v_cm::Array{Float64, 4}, collisions::MPCD)
    for i in eachindex(particle_vec)
        ix = Int64(div(particle_vec[i].pos[1] + args.box[1] / 2, collisions.a)) + 1
        iy = Int64(div(particle_vec[i].pos[2] + args.box[2] / 2, collisions.a)) + 1
        iz = Int64(div(particle_vec[i].pos[3] + args.box[3] / 2, collisions.a)) + 1

        #v_cm_cell = [v_cm[ix, iy, iz, 1], v_cm[ix, iy, iz, 2], v_cm[ix, iy, iz, 3]]
        #v_rel = velocity_vec[i].velocity .- v_cm_cell
        #E_k = 0.5 * sum(v_rel .^ 2)
        E_k = 0.5 * particle_vec[i].mass * sum(velocity_vec[i].velocity .^2)

        # Create a normal distribution with mean 0 and standard deviation sqrt(k_B*T/m)
        d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))

        # Draw velocities from Maxwell-Boltzmann distribution and calculate kinetic energy
        v_sample = [rand(d) for _ in 1:3]
        E_k_target = 0.5 * particle_vec[i].mass * sum(v_sample .^ 2)

        ξ = sqrt(E_k_target / E_k)

        # Update particle velocity according to MBS thermostat
        velocity_vec[i].velocity *= ξ 
    end
end

# Generate a rotation matrix around a random axis with a given angle. 
function generate_rotation_matrix(axis::Vector{Float64}, angle::Float64)
    u = normalize(axis)
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    one_minus_cos_angle = 1 - cos_angle

    R = Matrix{Float64}(LinearAlgebra.I, 3, 3)
    R[1, 1] = cos_angle + u[1]^2 * one_minus_cos_angle
    R[1, 2] = u[1] * u[2] * one_minus_cos_angle - u[3] * sin_angle
    R[1, 3] = u[1] * u[3] * one_minus_cos_angle + u[2] * sin_angle
    R[2, 1] = u[2] * u[1] * one_minus_cos_angle + u[3] * sin_angle
    R[2, 2] = cos_angle + u[2]^2 * one_minus_cos_angle
    R[2, 3] = u[2] * u[3] * one_minus_cos_angle - u[1] * sin_angle
    R[3, 1] = u[3] * u[1] * one_minus_cos_angle - u[2] * sin_angle
    R[3, 2] = u[3] * u[2] * one_minus_cos_angle + u[1] * sin_angle
    R[3, 3] = cos_angle + u[3]^2 * one_minus_cos_angle

    return R
end