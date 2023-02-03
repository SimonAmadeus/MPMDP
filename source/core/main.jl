export run_MD

include("system.jl")

# Main simulation.

function run_MD(args::System)
    
    println("You are currently using MPMDP, written by Simon Alberti (Version Feb23).")
    println(args.non_bonded_interactions, " interactions are used.")

    energy = 0.0e0
    virial = 0.0e0
    curr_p = 0.0e0
    #momentum = 0.0e0
    av_std_sum = 0.0e0
    av_sum = 0.0e0
    
    # Stress tensor: xx, yy, zz, xy, xz, yz
    stress = Vector(zeros(6))
    # Kinetic energy: xx, yy, zz
    E_kin = Vector(zeros(6)) 

    mesh = init_mesh(args.box, [10, 10, 10]) # Dummy mesh, used in hPF. 

    if typeof(args.bonded_interactions) != No_Bonds
        args.system_type = Molecular()
    end
    
    if typeof(args.non_bonded_interactions) == DPD
        if args.thermostat != No_Thermostat
            println("While using DPD, there is no additional thermostat needed!")
            args.thermostat = No_Thermostat()
        end
    elseif args.thermostat != No_Thermostat
        println(args.thermostat, " is used.")
    end
    
    # Read data
    particle_vec, c_l, velocity_vec, bond_vec, angle_vec, args.box = read_lammps_x(args.data_file, args.system_type)
    N = length(particle_vec)

    if typeof(args.sum) != No_Sum
        force_sum_x, force_sum_y, force_sum_z = initialize_force_cells(args, args.sum.a)
    end

    neighbor_list = []
    if typeof(args.nl) != No_NL
        for i in 1:N
            push!(neighbor_list, [])
        end
    end

    #c_l = Int(N / particle_vec[N].i_mol) # Chain length.
    
    # Slip-spring initialization.
    if typeof(args.t_bonds) != No_T_Bonds
        t_bond_vec = init_t_bonds!(args, c_l, particle_vec, bond_vec, args.t_bonds)
    end

    # Shear initialization.
    if typeof(args.shear) != No_Shear
        d_slab, z_bins = initialize_shear!(args, args.shear)
    end

    # hPF initialization.
    if typeof(args.non_bonded_interactions) == original_hPF_Interactions
        mesh = init_mesh(args.box, args.non_bonded_interactions.mesh_points)
        args.non_bonded_interactions = original_hPF_Interactions(args.non_bonded_interactions.κ, args.non_bonded_interactions.mesh_points)
    elseif typeof(args.non_bonded_interactions) == spectral_hPF_Interactions
        mesh = init_mesh(args.box, args.non_bonded_interactions.mesh_points)
        args.non_bonded_interactions = spectral_hPF_Interactions(args.non_bonded_interactions.κ, args.non_bonded_interactions.σ, args.non_bonded_interactions.mesh_points)
    end

    first_step = args.first_step

    # Open output files.

    log_file = open(args.log_file, "w")
    traj_file = open(args.traj_file, "w")

    if args.restart == true
        restart_file = open(args.restart_file, "w")
    end

    if typeof(args.shear) != No_Shear
        shear_file = open(args.shear_file, "w")
        momentum_file = open(args.momentum_file, "w")
    end

    if typeof(args.t_bonds) != No_T_Bonds
        entanglement_file = open("entanglement.out", "w")
    end

################################################################################

    # First simulation step.

################################################################################

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
        
        # Initialize force vector.
        force_vec = Vector{Force}(undef, length(particle_vec))
        for i = 1:length(particle_vec)
            force_vec[i] = Force(zeros(3))
        end

        #energy_vec = zeros(Float64, length(particle_vec))

        virial = stress[1] + stress[2] + stress[3]
        E_kin = calc_kinetic_energy(args, particle_vec, velocity_vec, E_kin)
        V = args.box[1] * args.box[2] * args.box[3]
        curr_p = 2 / V * (E_kin[1] + E_kin[2] + E_kin[3] + virial) / 3

        update_positions!(args, particle_vec, velocity_vec, force_vec)

        create_neighbor_list!(args, particle_vec, neighbor_list, args.nl)

        apply_logger(log_file, first_step, 0, velocity_vec, energy, curr_p, particle_vec, args.logger)
        apply_dump(traj_file, 0, args, particle_vec, force_vec, velocity_vec, args.trajectorydump)

        first_step = false

    end

################################################################################

    # Main simulation loop.

################################################################################

    @showprogress for step_i = 1:args.n_steps
    #for step_i = 1:args.n_steps
        # Shear flow (Momentum exchange). 
        if typeof(args.shear) != No_Shear && step_i%args.shear.freq == 0.0
            apply_shear!(args, particle_vec, velocity_vec, args.shear)
        end

        # Slip-spring chain end relocation and migration.
        if typeof(args.t_bonds) != No_T_Bonds && step_i%args.t_bonds.freq_mc == 0.0
            migration!(first_step, args, c_l, particle_vec, t_bond_vec, args.t_bonds)
        end
        
        # Apply first integration.
        update_velocities!(args, particle_vec, velocity_vec, force_vec)
        update_positions!(args, particle_vec, velocity_vec, force_vec)
        clear_force!(force_vec)

        # Apply interactions.
        energy, stress = apply_bonded_interactions!(args, particle_vec, force_vec, energy, stress, bond_vec, args.bonded_interactions)
        energy1 = energy
        energy, stress = apply_angle_interactions!(args, particle_vec, force_vec, energy, stress, angle_vec, args.angle_interactions)
        energy2 = energy - energy1 
        energy, stress = apply_non_bonded_interactions!(args, args.system_type, particle_vec, velocity_vec, c_l, force_vec, energy, stress, mesh, args.non_bonded_interactions)
        energy3 = energy - energy1 - energy2
        # Slip-spring interactions, if active.
        if typeof(args.t_bonds) != No_T_Bonds
            energy, stress = apply_bonded_interactions!(args, particle_vec, force_vec, energy, stress, t_bond_vec, args.t_bonds)
            energy4 = energy - energy1 - energy2 - energy3
        end
        # energy1: bonded
        # energy2: anngles
        # energy3: non-bonded
        # energy4: slip-springs
        
        # Apply second integration.
        update_velocities!(args, particle_vec, velocity_vec, force_vec)

        # Calcylate virial, kinetic energy and pressure.

        virial = stress[1] + stress[2] + stress[3]
        E_kin = calc_kinetic_energy(args, particle_vec, velocity_vec, E_kin)
        V = args.box[1] * args.box[2] * args.box[3]
        curr_p = 2 / V * (E_kin[1] + E_kin[2] + E_kin[3] + virial) / 3

        # Neighbor list is only created, if args.nl != No_NL.
        create_neighbor_list!(args, particle_vec, neighbor_list, args.nl)

        apply_thermostat!(args, particle_vec, velocity_vec, neighbor_list, args.thermostat)   
        apply_barostat!(args, particle_vec, curr_p, args.barostat)

        # Output.
        if step_i%args.log_freq == 0.0
            apply_logger(log_file, first_step, step_i, velocity_vec, energy, curr_p, particle_vec, args.logger)
        end
        if step_i%args.traj_freq == 0.0
            apply_dump(traj_file, step_i, args, particle_vec, force_vec, velocity_vec, args.trajectorydump)
        end
        if args.shear_sample_freq != 0 && step_i%args.shear_sample_freq == 0.0 && typeof(args.shear) != No_Shear
            z_bins = sample_velocity_profile!(args, particle_vec, velocity_vec, d_slab, z_bins, args.shear)
        end
        if step_i%args.n_steps == 0.0 && typeof(args.shear) != No_Shear
            apply_shear_log(shear_file, momentum_file, step_i, args, d_slab, z_bins, args.shear)
        end
        if step_i%args.n_steps == 0.0 && args.restart == true
            apply_restart(restart_file, c_l, args, particle_vec, velocity_vec, bond_vec, angle_vec, args.system_type)
        end
        if step_i%args.n_steps == 0.0 && typeof(args.t_bonds) != No_T_Bonds
            apply_t_bond_dump(entanglement_file, t_bond_vec)
        end
        if typeof(args.sum) != No_Sum && step_i%args.sum.freq == 0.0
            sum_force!(args, particle_vec, force_vec, force_sum_x, force_sum_y, force_sum_z, args.sum.a)
            println("Step $(step_i)")
            println(force_sum_x)
            println(std(force_sum_x))
            println(mean(force_sum_x))
            println(sum(force_sum_x))
            println(force_sum_y)
            println(std(force_sum_y))
            println(mean(force_sum_x))
            println(sum(force_sum_y))
            println(force_sum_z)
            println(std(force_sum_z))
            println(mean(force_sum_x))
            println(sum(force_sum_z))
            println("Average: ")
            av_std_sum += (std(force_sum_x) + std(force_sum_x) + std(force_sum_z)) / (3 * (args.n_steps / args.sum.freq))
            av_sum += (sum(force_sum_x) + sum(force_sum_x) + sum(force_sum_z)) / (3 * (args.n_steps / args.sum.freq))
            println(av_std_sum)
            println(av_sum)
        end

        # Clear energies after simulation step ends. 
	    energy = 0.0e0
        #clear_energy!(energy_vec)
        clear_stress_tensor!(stress)
        clear_kin_energy_tensor!(E_kin)
        if typeof(args.sum) != No_Sum
            clear_sums!(force_sum_x, force_sum_y, force_sum_z)
        end    
    end
end
