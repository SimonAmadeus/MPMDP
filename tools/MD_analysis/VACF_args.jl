# This code calculates the velocity autocorrelation function of a given trajectory.

ENV["GKSwstype"] = "100"
using Plots
using LinearAlgebra
using ProgressMeter
using Statistics

n_frames = parse(Int64, ARGS[1])
dt = parse(Int64, ARGS[2])

traj_file = "trj.xyz"
out_file = "vacf.txt"

fh = 2 # File header.
n_particles = 0 # Number of particles.

function autocorrelation(v, t)
    mean_v = mean(v)
    var_v = var(v)
    function C(τ)
        N = length(v) - τ
        sum = 0.0
        for t = 1:N
            sum += (v[t] - mean_v) * (v[t+τ] - mean_v)
        end
        return sum / ( N * var_v)
    end
    return [C(τ) for τ = 0:t-1]
end

function average_autocorrelation(vx_data, vy_data, vz_data, t)
    autocorrs = [(autocorrelation(vx, t) .+ autocorrelation(vy, t) .+ autocorrelation(vz, t)) ./ 3 for (vx, vy, vz) in zip(eachcol(vx_data), eachcol(vy_data), eachcol(vz_data))]
    return mean(autocorrs)
end

function diffusion_coefficient(autocorr, Δt)
    integral = sum(autocorr) * Δt
    return integral / 3  # Divide by 3 for 3 dimensions.
end

open(traj_file) do file
    global N = parse(Int64, readline(file))
    println("$N particles in the system.")
end

vx = zeros(n_frames, N)
vy = zeros(n_frames, N)
vz = zeros(n_frames, N)

open(traj_file) do file
    for (t, line) in enumerate(eachline(file))
        m = mod(t, fh + N)
        frame = (t ÷ (fh + N)) + 1
        if (m != 1 && m != 2 ) && frame <= n_frames
            a = split(line)
            ind = m - fh
            if m == 0
                ind = N
            end
            vx[frame, ind] = parse(Float64, a[5])
            vy[frame, ind] = parse(Float64, a[6])
            vz[frame, ind] = parse(Float64, a[7])
        end
    end
end

# Substract center of mass velocity.
for i in 1:n_frames
    vx[i, :] .-= mean(vx[i, :])
    vy[i, :] .-= mean(vy[i, :])
    vz[i, :] .-= mean(vz[i, :])
end

# Calculate the autocorrelation and diffusion coefficient.
autocorr = average_autocorrelation(vx, vy, vz, n_frames)
D = diffusion_coefficient(autocorr, dt)

println("Diffusion coefficient: ", D)

out_file = open(out_file, "w")
try
    for (i, C) in enumerate(autocorr)
        write(out_file, "$((i-1)*dt) $C\n")
    end
finally
    close(out_file)
end

# plotting:
default(size=(800, 600), linewidth=2, markersize=4, legendfontsize=12, guidefontsize=14, tickfontsize=12)
gr()

time = collect(0:length(autocorr) - 1)
time *= dt

# Create the VACF plot
p = plot(time, autocorr,
    line=:solid,
    marker=:circle,
    markercolor=:blue,
    markeralpha=0.5,
    markersize=4,
    xlabel="Timesteps",
    ylabel="VACF",
    title="Velocity Autocorrelation Function",
    legend=:none,
    grid=:off,
    framestyle=:box)

# Save the VACF plot
savefig(p, "VACF.png")