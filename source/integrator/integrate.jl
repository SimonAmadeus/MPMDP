# Velocity Verlet intgratior.
function update_positions!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force})
    # Update positions.
    @inbounds @simd for i = 1:length(particle_vec)
        particle_vec[i].pos .+= velocity_vec[i].velocity .* args.∆t
    end 
    # Update images.
    @inbounds @simd for i = 1:length(particle_vec)
        particle_vec[i].image .+= floor.((particle_vec[i].pos .+ args.box ./ 2) ./ args.box)
        particle_vec[i].pos   .-= args.box .* floor.((particle_vec[i].pos .+ args.box ./ 2) ./ args.box)
    end
end

function update_velocities!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force})
    # Update velocities.
    @inbounds @simd for i = 1:length(particle_vec)
        velocity_vec[i].velocity .+= 0.5 .* (force_vec[i].force ./ particle_vec[i].mass) .* args.∆t
    end 
end
#=
function update_positions!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force})

    # Update positions.
    for i = 1:length(particle_vec)
        particle_vec[i].pos[1] += velocity_vec[i].velocity[1] * args.∆t
        particle_vec[i].pos[2] += velocity_vec[i].velocity[2] * args.∆t
        particle_vec[i].pos[3] += velocity_vec[i].velocity[3] * args.∆t
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

function update_velocities!(args::System, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, force_vec::Vector{Force})

    # Update velocities.
    for i = 1:length(particle_vec)
        velocity_vec[i].velocity[1] += 0.5 * force_vec[i].force[1] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[2] += 0.5 * force_vec[i].force[2] / particle_vec[i].mass * args.∆t
        velocity_vec[i].velocity[3] += 0.5 * force_vec[i].force[3] / particle_vec[i].mass * args.∆t
    end 

end=#