export run_MD

include("system.jl")

# Main simulation.
function run_MD(args::System)

    # Read data
    particle_vec, l_chain, velocity_vec, bond_vec, angle_vec, args.box = read_lammps_x(args.data_file, args.system_type)
    #particle_vec, args.box = read_lammps_x(args.data_file, args.system_type)
    
    mesh = init_mesh(args.box, 10, 10, 10)
    if typeof(args.non_bonded_interactions) == original_hPF_Interactions
        args.non_bonded_interactions = original_hPF_Interactions(mesh, 0.2)
    elseif typeof(args.non_bonded_interactions) == spectral_hPF_Interactions
        args.non_bonded_interactions = spectral_hPF_Interactions(mesh, 0.2, 0.5)
    end

    first_step = args.first_step

    log_file = open(args.log_file, "w")
    traj_file = open(args.traj_file, "w")

    if typeof(args.shear) != No_Shear
        shear_file = open(args.shear_file, "w")
        momentum_file = open(args.momentum_file, "w")
    end
    if args.restart == true
        restart_file = open(args.restart_file, "w")
    end

    # First simulation step.
    if first_step == true

        # initialize veclocities!
        if args.init_velocity == "Z"
            for i = 1:length(particle_vec)
                velocity_vec[i] = Velocity(zeros(3))
            end
        elseif args.init_velocity == "G"
            for i = 1:length(particle_vec)
                d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
                velocity_vec[i] = Velocity(rand(d,3))
            end
        elseif args.init_velocity == "R"
            println("Initial velocities are read from restart file.") 
        else
            println("Please set a mode for the initial velocity!")
        end
        
        force_vec = Vector{Force}(undef, length(particle_vec))
        for i = 1:length(particle_vec)
            force_vec[i] = Force(zeros(3))
        end

        energy_vec = zeros(Float64, length(particle_vec))

        apply_1st_integration!(args, particle_vec, velocity_vec, force_vec, args.thermostat)
        
        apply_logger!(log_file, first_step, 1, velocity_vec, energy_vec, particle_vec, args.logger)
        #apply_dump(args.traj_file, 1, args, particles, forces, velocities, args.trajectorydump)
        first_step = false

    end

    @showprogress for step_i = 2:args.n_steps
        
        if typeof(args.shear) != No_Shear && step_i%args.shear_freq == 0.0
            apply_shear(args, particle_vec, velocity_vec, args.shear)
        end
        
        # Apply interactions.
        apply_bonded_interactions!(args, particle_vec, force_vec, energy_vec, bond_vec, args.bonded_interactions)
        apply_non_bonded_interactions!(args, particle_vec, force_vec, energy_vec, args.non_bonded_interactions)
        
        # Apply integration.
        apply_1st_integration!(args, particle_vec, velocity_vec, force_vec, args.thermostat)
        apply_2nd_integration!(args, particle_vec, velocity_vec, force_vec, args.thermostat)
        
        # Clear forces and energies after integration. 
        clear_force!(force_vec)
        clear_energy!(energy_vec)
        
        if step_i%1000 == 0
            VAR = 1
        end

        # Output.
        if step_i%args.log_freq == 0.0
            apply_logger!(log_file, first_step, step_i, velocity_vec, energy_vec, particle_vec, args.logger)
        end

        if step_i%args.traj_freq == 0.0
            apply_dump(traj_file, step_i, args, particle_vec, force_vec, velocity_vec, args.trajectorydump)
        end

        if step_i%args.n_steps == 0.0 && typeof(args.shear) != No_Shear
            sample_velocity_profile(shear_file, momentum_file, step_i, args, particle_vec, velocity_vec, args.shear)
        end

        if step_i%args.n_steps == 0.0 && args.restart == true
            apply_restart!(restart_file, l_chain, args, particle_vec, velocity_vec, bond_vec, angle_vec, args.system_type)
        end
        
    end
end