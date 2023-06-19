# This code calculates the Radius of gyration of a polymer chain.

using Plots
using LinearAlgebra
using ProgressMeter
using Statistics

#println("Enter input file (XYZ format).")
#traj_file = chomp(readline())
#println("Enter out-file (e.g. rog.txt).")
#out_file = chomp(readline())
#println("How many frames would you like to analyze?")
#n_frames = parse(Int64, chomp(readline()))
println("Enter chain length?")
len_chn = parse(Int64, chomp(readline()))

traj_file = "trj.xyz"
out_file = "rog.txt"
n_frames = 100

counter = 0
t = 0
fh = 2 # File header.
N = 0 # Number of particles.

#box = [xlen, ylen, zlen]

open(traj_file) do file
    global N = parse(Int64, split(readline(file))[1])
    global n_chains = Int(N / len_chn)
    println("$n_chains chains in the system.")
end

x = zeros(n_frames, N)
y = zeros(n_frames, N)
z = zeros(n_frames, N)
xcom = zeros(n_frames, n_chains)
ycom = zeros(n_frames, n_chains)
zcom = zeros(n_frames, n_chains)

# Read x, y and z coordinates of all frames.
open(traj_file) do file
    for (line_number, line) in enumerate(eachline(file))
        m = mod(line_number, fh + N)
        frame = (line_number ÷ (fh + N)) + 1
        if (m != 1 && m != 2 && frame <= n_frames)
            a = split(line)
            ind = m - fh
            if(m == 0)
                ind = N
                frame -= 1
            end
            x[frame, ind] = parse(Float64, a[2])
            y[frame, ind] = parse(Float64, a[3])
            z[frame, ind] = parse(Float64, a[4])                  
        end        
    end
end

for i in 1:n_chains
    b = (i - 1) * len_chn + 1 # First bead if the chain.
    c = i * len_chn # Last bead of the chain.
    for frame in 1:n_frames
        xcom[frame, i] = mean(x[frame, b:c])
        ycom[frame, i] = mean(y[frame, b:c])
        zcom[frame, i] = mean(z[frame, b:c])
    end
end

rog = zeros(n_frames)
rogx = zeros(n_frames)
rogy = zeros(n_frames)
rogz = zeros(n_frames)

for frame in 1:n_frames
    for i in 1:n_chains
        b = (i - 1) * len_chn + 1
        c = i * len_chn    
        for j in b:c
            rog[frame] += ((x[frame, j] - xcom[frame, i])^2 + (y[frame, j] - ycom[frame, i])^2 + (z[frame, j] - zcom[frame, i])^2) / (n_chains * len_chn)
            rogx[frame] += ((x[frame, j] - xcom[frame, i])^2) / (n_chains * len_chn)
            rogy[frame] += ((y[frame, j] - ycom[frame, i])^2) / (n_chains * len_chn)
            rogz[frame] += ((z[frame, j] - zcom[frame, i])^2) / (n_chains * len_chn)
        end
    end
end

out = open(out_file, "w")
string_out = string(mean(rog)) * string(var(rog)) * "\n"
write(out, string_out)
string_out = string(mean(rogx)) * string(var(rogx)) * "\n"
write(out, string_out)
string_out = string(mean(rogy)) * string(var(rogy)) * "\n"
write(out, string_out)
string_out = string(mean(rogz)) * string(var(rogz)) * "\n"
write(out, string_out)

close(out)

println("ROG^2 (x) = ", mean(rogx), " std = ", var(rogx))
println("ROG^2 (y) = ", mean(rogy), " std = ", var(rogy))
println("ROG^2 (z) = ", mean(rogz), " std = ", var(rogz))
println("ROG^2 (total) = ", mean(rog), " std = ", var(rog))
