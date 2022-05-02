export XYZTrajectoryDump

function apply_dump(fileO, step::Int64, args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, velocity_vec::Vector{Velocity}, ::No_TrajectoryDump) end

struct XYZTrajectoryDump <: TrajectoryDump
end

function apply_dump(fileO, step::Int64, args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, velocity_vec::Vector{Velocity}, ::XYZTrajectoryDump) 
    
    x = 0 
    y = 0
    z = 0
    
    string_out=""
    string_out=string(length(particle_vec))
    string_out*="\n"
    write(fileO, string_out)
    string_out=""
    string_out*=string(step)
    string_out*="\n"
    write(fileO, string_out)
    for i = 1:length(particle_vec)
        string_out=""
        type_ = particle_vec[i].spec
        x = particle_vec[i].pos[1] + (particle_vec[i].image[1] * args.box[1])
        y = particle_vec[i].pos[2] + (particle_vec[i].image[2] * args.box[2])
        z = particle_vec[i].pos[3] + (particle_vec[i].image[3] * args.box[3])
        string_out = string(type_)*" "*string(x)*" "*string(y)*" "*string(z)
        string_out*="\n"
        write(fileO, string_out)
    end
end