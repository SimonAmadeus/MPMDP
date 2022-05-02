export Berendsen_Thermostat, Andersen_Thermostat

function calc_totalkineticenergy(particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity})
    K_tot::Float64 = 0
    for i = 1:length(velocity_vec)
        for j = 1:3
            K_tot += sum(velocity_vec[i].velocity[j]^2 * 0.5 * particle_vec[i].mass)
        end
    end
    return K_tot
end

function calc_inst_temp(particle_vec::Vector{Particle}, K_tot::Float64)
    dim::Int64 = 3
    T::Float64 = K_tot * 2 / (length(particle_vec) * dim - dim)
    return T
end

function calc_totalmomentum(particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity})
    P_tot::Float64 = 0
    for i = 1:length(velocity_vec)
        P_tot += sum(velocity_vec[i].velocity .* particle_vec[i].mass)
    end
    return P_tot
end

function apply_1st_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, ::No_Thermostat) end

function apply_2nd_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, ::No_Thermostat) end

struct Berendsen_Thermostat <: Thermostat 
    τ::Float64
end

function apply_1st_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, Thermostat::Berendsen_Thermostat) 
    # Update positions.
    for i = 1:length(particle_vec)
        particle_vec[i].pos[1] += (velocity_vec[i].velocity[1] + 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t) * args.∆t
        particle_vec[i].pos[2] += (velocity_vec[i].velocity[2] + 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t) * args.∆t
        particle_vec[i].pos[3] += (velocity_vec[i].velocity[3] + 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t) * args.∆t
    end 

    # Update velocities.
    for i = 1:length(particle_vec)
        velocity_vec[i].velocity[1] += 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[2] += 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[3] += 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t
    end 

    # Update images.
    for i = 1:length(particle_vec)
        particle_vec[i].image[1] += Int(floor((particle_vec[i].pos[1] + args.box[1] / 2) / args.box[1]))
        particle_vec[i].image[2] += Int(floor((particle_vec[i].pos[2] + args.box[2] / 2) / args.box[2]))
        particle_vec[i].image[3] += Int(floor((particle_vec[i].pos[3] + args.box[3] / 2) / args.box[3]))
        particle_vec[i].pos[1]   -= args.box[1] * floor((particle_vec[i].pos[1] + args.box[1] / 2) / args.box[1])
        particle_vec[i].pos[2]   -= args.box[2] * floor((particle_vec[i].pos[2] + args.box[2] / 2) / args.box[2])
        particle_vec[i].pos[3]   -= args.box[3] * floor((particle_vec[i].pos[3] + args.box[3] / 2) / args.box[3])
    end 
end

function apply_2nd_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, Thermostat::Berendsen_Thermostat) 

    totalkineticenergy = calc_totalkineticenergy(particle_vec, velocity_vec)
    curr_T = calc_inst_temp(particle_vec, totalkineticenergy)
    scalT  = sqrt(1.0 + args.∆t * (args.T / curr_T - 1.0) / Thermostat.τ)
    
    # update regulated velocities
    for i = 1:length(particle_vec)
        velocity_vec[i].velocity[1] = (velocity_vec[i].velocity[1] + 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t) * scalT
        velocity_vec[i].velocity[2] = (velocity_vec[i].velocity[2] + 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t) * scalT
        velocity_vec[i].velocity[3] = (velocity_vec[i].velocity[3] + 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t) * scalT
    end 
end

struct Andersen_Thermostat <: Thermostat 
    σ::Float64
end

function apply_1st_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, Thermostat::Andersen_Thermostat) 
    # Update positions.
    for i = 1:length(particle_vec)
        particle_vec[i].pos[1] += (velocity_vec[i].velocity[1] + 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t) * args.∆t
        particle_vec[i].pos[2] += (velocity_vec[i].velocity[2] + 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t) * args.∆t
        particle_vec[i].pos[3] += (velocity_vec[i].velocity[3] + 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t) * args.∆t
    end 

    # Update velocities.
    for i = 1:length(particle_vec)
        velocity_vec[i].velocity[1] += 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[2] += 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[3] += 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t
    end 

    # Update images.
    for i = 1:length(particle_vec)
        particle_vec[i].image[1] += fld(particle_vec[i].pos[1], args.box[1])
        particle_vec[i].image[2] += fld(particle_vec[i].pos[2], args.box[2])
        particle_vec[i].image[3] += fld(particle_vec[i].pos[3], args.box[3])
        particle_vec[i].pos[1]   -= args.box[1] * fld(particle_vec[i].pos[1], args.box[1])
        particle_vec[i].pos[2]   -= args.box[2] * fld(particle_vec[i].pos[2], args.box[2])
        particle_vec[i].pos[3]   -= args.box[3] * fld(particle_vec[i].pos[3], args.box[3])
    end 
end

function apply_2nd_integration!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force}, Thermostat::Andersen_Thermostat) 
    # Update regulated velocities.
    for i = 1:length(particle_vec)
        d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
        rn = rand(1)
        if rn[1] < Thermostat.σ * args.∆t
            velocity_vec[i].velocity = rand(d,3) # randn(rng,Float64) * sqrt(args.temp/atoms[velii].mass)
        else
            velocity_vec[i].velocity[1] += 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t
            velocity_vec[i].velocity[2] += 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t
            velocity_vec[i].velocity[3] += 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t
        end
    end
end