export Berendsen_Thermostat, Andersen_Thermostat, Lowe_Andersen_Thermostat

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, ::No_Thermostat) end

struct Berendsen_Thermostat <: Thermostat 
    τ::Float64
end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, Thermostat::Berendsen_Thermostat)

    K_tot = calc_total_kinetic_energy(particle_vec, velocity_vec)
    curr_T = calc_inst_temp(particle_vec, K_tot)
    scalT  = sqrt(1.0 + args.∆t * (args.T / curr_T - 1.0) / Thermostat.τ)
    
    # Update velocities.
    for i = 1:length(particle_vec)
        velocity_vec[i].velocity[1] *= scalT
        velocity_vec[i].velocity[2] *= scalT
        velocity_vec[i].velocity[3] *= scalT
    end

end

struct Andersen_Thermostat <: Thermostat 
    σ::Float64
end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, Thermostat::Andersen_Thermostat)

    # Update velocities.
    for i = 1:length(particle_vec)
        rn = rand(1)
        if rn[1] < Thermostat.σ * args.∆t
            d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
            velocity_vec[i].velocity = rand(d,3) # randn(rng,Float64) * sqrt(args.temp/atoms[velii].mass)
        #else
        #    velocity_vec[i].velocity[1] += 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t
        #    velocity_vec[i].velocity[2] += 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t
        #    velocity_vec[i].velocity[3] += 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t
        end
    end

end

struct Lowe_Andersen_Thermostat
    Γ::Float64
end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, Thermostat::Lowe_Andersen_Thermostat)

    # Update velocities.
    for i = 1:length(neighbor_list)
        for j = 1:length(neighbor_list[i])
            rn = rand(1)
            if rn[1] < Thermostat.Γ * args.∆t
                d = Normal(0.0, sqrt(2 * args.T / particle_vec[i].mass))
                ∆d = particle_vec[neighbor_list[i][j]].pos .- particle_vec[i].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                r_unit = ∆d / r
                ∆v = velocity_vec[neighbor_list[i][j]].velocity .- velocity_vec[i].velocity 
                ∆v_relative = rand(d, 3)
                ∆ij = 0.5 * r_unit * (∆v_relative - ∆v) * r_unit
                velocity_vec[i].velocity += ∆ij
                velocity_vec[neighbor_list[i][j]].velocity -= ∆ij
            end    
        end
    end
end