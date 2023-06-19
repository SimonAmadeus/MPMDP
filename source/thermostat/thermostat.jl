export Berendsen_Thermostat, Andersen_Thermostat, Lowe_Andersen_Thermostat, Langevin_Thermostat

# Berendsen thermostat.
struct Berendsen_Thermostat <: Thermostat 
    τ::Float64
end

# Andersen thermostat.
struct Andersen_Thermostat <: Thermostat 
    σ::Float64
end

# Lowe-Andersen thermostat.
struct Lowe_Andersen_Thermostat <: Thermostat
    Γ::Float64
    r_c::Float64
end

struct Langevin_Thermostat <: Thermostat
    γ::Float64
end

# No Thermostat:
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, ::No_Thermostat) end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, Thermostat::Berendsen_Thermostat)
    # Determine scaling factor.
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

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neihbor_list::Any, Thermostat::Andersen_Thermostat)
    # Update velocities.
    for i = 1:length(particle_vec)
        rn = rand(1)
        if rn[1] < Thermostat.σ * args.∆t
            d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
            velocity_vec[i].velocity = rand(d, 3)
        end
    end

end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neighbor_list::Any, Thermostat::Lowe_Andersen_Thermostat)
    # Update velocities in a pairwise fashion.
    for i = 1:length(neighbor_list)
        for j = 1:length(neighbor_list[i])
            rn = rand(1)
            if rn[1] < Thermostat.Γ * args.∆t
                d = Normal(0.0, sqrt(2 * args.T / particle_vec[i].mass))
                ∆d = particle_vec[i].pos .- particle_vec[neighbor_list[i][j]].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r < Thermostat.r_c
                    r_unit = ∆d / r
                    ∆v = velocity_vec[i].velocity .- velocity_vec[neighbor_list[i][j]].velocity 
                    v_random = rand(d, 3)
                    ∆v_relative = v_random - ∆v
                    ∆ij = 0.5 * r_unit * dot(∆v_relative, r_unit)
                    velocity_vec[i].velocity += ∆ij
                    velocity_vec[neighbor_list[i][j]].velocity -= ∆ij
                end
            end    
        end
    end
end

function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neighbor_list::Any, thermostat::Langevin_Thermostat)
    Random.seed!(1234)
    for i = 1:length(particle_vec)
        mass = particle_vec[i].mass
        σ = sqrt(2 * thermostat.γ * args.T / mass)
        
        f_r = σ * randn(3)
        f_f = - thermostat.γ * mass * velocity_vec[i].velocity
        f = f_f + f_r
        a = f / mass

        # Update velocities using the Langevin equation
        velocity_vec[i].velocity .+= a * args.∆t
    end
end