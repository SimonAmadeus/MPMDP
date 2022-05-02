# Initialize coordinates in lammps format. Written by Simon Alberti in January 2022.

using LinearAlgebra
using ProgressMeter

mutable struct Particle
    spec::Int64
    mass::Float64
    charge::Float64
    pos::Vector{Float64}
    image::Vector{Int64}
end

mutable struct Bond
    type::Int64
    p_1::Int64
    p_2::Int64
end

mutable struct Angle
    type::Int64
    p_1::Int64
    p_2::Int64
    p_3::Int64
end

function correct_distance_pbc(box, d)
    if d[1] > box[1] / 2
        d[1] -= box[1]
    elseif d[1] < -box[1] / 2
        d[1] += box[1]
    end
    if d[2] > box[2] / 2
        d[2] -= box[2]
    elseif d[2] < -box[2] / 2
        d[2] += box[2]
    end
    if d[3] > box[3] / 2
        d[3] -= box[3]
    elseif d[3] < -box[3] / 2
        d[3] += box[3]
    end
    return d
end

#mutable struct system
#    box::Vector{Float64}
#end

#function sys_init()
#    box = [0.0, 0.0, 0.0]
#end

box = [0.0, 0.0, 0.0]

println("Enter the box size in x-direction.")
box[1] = parse(Float64, chomp(readline()))
println("Enter the box size in y-direction.")
box[2] = parse(Float64, chomp(readline()))
println("Enter the box size in z-direction.")
box[3] = parse(Float64, chomp(readline()))

V = box[1] * box[2] * box[3]

#println("$box[1] $box[2] $box[2]")

chains = ""

println("Polymer chain? (y/n)")
chains = chomp(readline())
if chains != "y" && chains != "n"
    println("Only y and n dumbass!")
    exit()
end

# NO MOLECULES.
if chains == "n"
    n_spec = 1
    #println("Enter number of different species.")
    #n_spec = chomp(readline())
    println("Enter total number of particles.")
    n_particles = parse(Int64, chomp(readline()))

    println("Enter tolerance distance σ.")
    σ = parse(Float64, chomp(readline()))

    #ρ = n_particle_vec / V
    spec = 1
    mass = 100
    charge = 0
    coords = [0.0, 0.0 , 0.0]
    image = [0, 0, 0]

    particle_vec = Vector{Particle}(undef, n_particles) # Initialize particle vector.

    check = false

    # Generate data.

    for i in 1:n_particles
        global check = false
        if i == 1 # Place first particle randomly within the box.
            global coords = rand(Float64, 3) .* box .- box / 2 
            particle_i = Particle(spec, mass, charge, coords, image)
            particle_vec[i] = particle_i
        else
            while check == false # Create random coordinates unti r > σ.
                global check = true
                global coords = rand(Float64, 3) .* box .- box / 2 
                for j in 1:(i - 1) 
                    d = particle_vec[j].pos .- coords
                    d = correct_distance_pbc(box, d)
                    r = norm(d)
	                if r < σ
                        global check = false
		            end
                end
            end
            particle_i = Particle(spec, mass, charge, coords, image)
            particle_vec[i] = particle_i
        end
    end

    # Write output.
    open("input.data", "w") do file
        write(file, "LAMMPS data file generated via init_lammps.jl\n \n")
        write(file, "$n_particles atoms\n")
        write(file, "$n_spec atom types\n")
        write(file, "\n")
        write(file, "$(-box[1]/2) $(box[1]/2) xlo xhi\n")
        write(file, "$(-box[2]/2) $(box[2]/2) ylo yhi\n")
        write(file, "$(-box[3]/2) $(box[3]/2) zlo zhi\n")
        write(file, "\n")
        write(file, "Masses\n")
        write(file, "\n")
        for i in 1:n_spec
            write(file, "$i $mass\n")
        end
        write(file, "\n")
        write(file, "Atoms\n")
        write(file, "\n")
        for i in 1:n_particles
            write(file, "$i $(particle_vec[i].spec) $(particle_vec[i].spec) $(particle_vec[i].charge) $(particle_vec[i].pos[1]) $(particle_vec[i].pos[2]) $(particle_vec[i].pos[3]) $(particle_vec[i].image[1]) $(particle_vec[i].image[2]) $(particle_vec[i].image[3])\n")
        end
    end

# Polymer chains.
elseif chains == "y"
    println("Enter total number of particles.")
    n_particles = parse(Int64, chomp(readline()))
    println("Enter chain length.")
    l_chain = parse(Int64, chomp(readline()))
    n_mol = n_particles / l_chain
    n_mol = Int(n_mol)
    if isinteger(n_mol) == false
        println("Check you numbers, dumbass!")
        exit()
    end
    n_bonds = (l_chain - 1) * n_mol
    n_bonds = Int(n_bonds)
    n_angles = (l_chain - 2) * n_mol
    n_angles = Int(n_angles)
    println("Enter tolerance distance σ.")
    σ = parse(Float64, chomp(readline()))

    l_bond = 1
    spec = 1
    mass = 100
    charge = 0
    coords = [0.0, 0.0 , 0.0]
    image = [0, 0, 0]

    particle_vec = Vector{Particle}(undef, n_particles)
    #bond_vec = Vector{Bond}(undef, n_bonds)
    #angle_vec = Vector{Angle}(undef, n_angles)

    check = false
    i_counter = 0

    for i in 1:n_mol
        for j in 1:l_chain
            global check = false
            global i_counter += 1
            if i_counter == 1 # Place first particle randomly.
                global coords = rand(Float64, 3) .* box .- box / 2 
                particle_i = Particle(spec, mass, charge, coords, image)
                particle_vec[i_counter] = particle_i
            elseif j == 1 && i_counter != 1 # At new chain, check tolerance distance from all other particles.
                while check == false
                    global check = true
                    global coords = rand(Float64, 3) .* box .- box / 2 
                    for k in 1:(i_counter - 1) 
                        d = particle_vec[k].pos .- coords
                        d = correct_distance_pbc(box, d)
                        r = norm(d)
	                    if r < σ
                            global check = false
		                end
                    end
                end
                particle_i = Particle(spec, mass, charge, coords, image)
                particle_vec[i_counter] = particle_i
            else
                while check == false # Create positions with a distance l_bond from previous particle and check tolrecance.
                    global check = true
                    global coords = rand(Float64, 3) .* box .- box / 2
                    len = norm(coords)
                    coords = coords ./ len
                    coords = particle_vec[i_counter - 1].pos .+ coords
                    for k in 1:(i_counter - 1) 
                        d = particle_vec[k].pos .- coords
                        d = correct_distance_pbc(box, d)
                        r = norm(d)
	                    if r < σ
                            global check = false
		                end
                    end
                end
                particle_i = Particle(spec, mass, charge, coords, image)
                particle_vec[i_counter] = particle_i
            end
        end
    end

    # Write output.
    open("input.data", "w") do file
        write(file, "LAMMPS data file generated via init_lammps_x.jl. Written by Simon Alberti.\n \n")
        write(file, "$n_particles atoms\n")
        n_spec = 1
        write(file, "$n_spec atom types\n")
        write(file, "$n_bonds bonds\n")
        n_spec = 1
        write(file, "$n_spec bond types\n")
        write(file, "$n_angles angles\n")
        n_spec = 1
        write(file, "$n_spec angle types\n")
        write(file, "\n")
        
        write(file, "$(-box[1]/2) $(box[1]/2) xlo xhi\n")
        write(file, "$(-box[2]/2) $(box[2]/2) ylo yhi\n")
        write(file, "$(-box[3]/2) $(box[3]/2) zlo zhi\n")
        
        write(file, "\n")
        write(file, "Masses\n")
        write(file, "\n")
        for i in 1:n_spec
            write(file, "$i $mass\n")
        end
        write(file, "\n")

        write(file, "Atoms\n")
        write(file, "\n")
        
        i_counter = 0

        for i in 1:n_mol
            for j in 1:l_chain
                i_counter += 1
                write(file, "$i_counter $i $(particle_vec[i_counter].spec) $(particle_vec[i_counter].charge) $(particle_vec[i_counter].pos[1]) $(particle_vec[i_counter].pos[2]) $(particle_vec[i_counter].pos[3]) $(particle_vec[i_counter].image[1]) $(particle_vec[i_counter].image[2]) $(particle_vec[i_counter].image[3])\n")
            end
        end
        
        write(file, "\n")
        write(file, "Bonds\n")
        write(file, "\n")
        
        cl_counter = 1 # Chain length counter.
        global i_counter = 0 # Index counter.
        global i_counter2 = 0 # Bond/Angle counter.
        bond_type = 1

        for i in 1:n_mol
            global cl_counter = 1
            for j in 1:l_chain
                global i_counter += 1
                if cl_counter != l_chain
                    global cl_counter += 1
                    global i_counter2 += 1
                    write(file, "$i_counter2 $bond_type $i_counter $(i_counter + 1)\n")
                end
            end
        end

        write(file, "\n")
        write(file, "Angles\n")
        write(file, "\n")
        
        global i_counter = 0 # Index counter.
        global i_counter2 = 0 # Bond/Angle counter.

        for i in 1:n_mol
            global cl_counter = 1
            for j in 1:l_chain
                global i_counter += 1
                if cl_counter < l_chain - 1
                    global cl_counter += 1
                    global i_counter2 += 1
                    write(file, "$i_counter2 $bond_type $i_counter $(i_counter + 1) $(i_counter + 2)\n")
                end
            end
        end
    end

else
    println("Only y or n you dumbass.")
end