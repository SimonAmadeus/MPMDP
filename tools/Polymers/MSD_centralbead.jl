#!/usr/bin/env julia
# This code calculates the MSD of the central bead of a polymer chain.
ENV["GKSwstype"] = "100"
using Plots
using LinearAlgebra
using ProgressMeter
using Statistics

#println("Enter input file (XYZ format).")
#traj_file = chomp(readline())
#println("Enter out-file (e.g. msd.txt).")
#out_file = chomp(readline())
#println("How many frames would you like to analyze?")
#n_frames = parse(Int64, chomp(readline()))
println("Enter chain length?")
len_chn = parse(Int64, chomp(readline()))

traj_file = "trj.xyz"
out_file = "msd.txt"
n_frames = 1000

dt = 10
fh = 2 # File header.

open(traj_file) do file
    global N = parse(Int64, split(readline(file))[1])
    global n_chains = Int(N / len_chn)
    println("$n_chains chains in the system.")
end

x = zeros(n_frames, N)
y = zeros(n_frames, N)
z = zeros(n_frames, N)

xcb = zeros(n_frames, n_chains)
ycb = zeros(n_frames, n_chains)
zcb = zeros(n_frames, n_chains)

open(traj_file) do file
    t = 0
    for line in eachline(file)
        t += 1
        m = mod(t, fh + N)
        frame = (t ÷ (fh + N)) + 1
        if (m != 1 && m != 2 )
            a = split(line)
            ind = m - fh
            if(m == 0)
                ind = N
                frame -= 1
            end
            if frame <= n_frames
                x[frame, ind] = parse(Float64, a[2])
                y[frame, ind] = parse(Float64, a[3])
                z[frame, ind] = parse(Float64, a[4])
            else
                break
            end                  
        end        
    end
end

for i in 1:n_chains
    b = (i - 1) * len_chn + 1 # First bead if the chain.
    c = i * len_chn # Last bead of the chain.
    for frame in 1:n_frames
        xcb[frame, i] = x[frame, Int(ceil(0.5 * (c - b)))]
        ycb[frame, i] = y[frame, Int(ceil(0.5 * (c - b)))]
        zcb[frame, i] = z[frame, Int(ceil(0.5 * (c - b)))]
    end
end


l_msd = n_frames ÷ 2
msd = zeros(l_msd)
for i in 1:(n_frames-l_msd)
    l = 0
    for j in i + 1:(i + l_msd)
        l += 1
        for k in 1:n_chains
            msd[l] += (xcb[j, k] - xcb[i, k])^2 + (ycb[j, k] - ycb[i, k])^2 + (zcb[j, k] - zcb[i, k])^2
        end
    end         
end

time = collect(1:l_msd)
time *= dt

msd /= n_chains * l_msd

out_file = open(out_file, "w")
string_out = "t msd\n"
write(out_file, string_out)
for i in 1:length(msd)
    global string_out = string(time[i]) * " " * string(msd[i]) * "\n"
    write(out_file, string_out)
end

# Plotting:
default(size=(800, 600), linewidth=2, markersize=4, legendfontsize=12, guidefontsize=14, tickfontsize=12)
gr()

# Create the plot
p = plot(time, msd, 
    line=:solid, 
    marker=:circle, 
    markercolor=:blue, 
    markeralpha=0.5, 
    markersize=4, 
    xlabel="timesteps", 
    ylabel="MSD", 
    title="MSD", 
    legend=:none, 
    grid=:off, 
    framestyle=:box,
    xaxis=:log,
    yaxis=:log)

# Save the plot
savefig(p, "msd_cb.png")
