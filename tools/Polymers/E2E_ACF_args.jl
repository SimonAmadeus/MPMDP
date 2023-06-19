ENV["GKSwstype"] = "100"
using Plots
using LinearAlgebra
using ProgressMeter
using Statistics

len_chn = parse(Int64, ARGS[1])
n_frames = parse(Int64, ARGS[2])
dt = parse(Int64, ARGS[3])

traj_file = "trj.xyz"
out_file = "e2e_acf.txt"

fh = 2 # File header.

open(traj_file) do file
    global N = parse(Int64, split(readline(file))[1])
    global n_chains = Int(N / len_chn)
    println("$n_chains chains in the system.")
end

x = zeros(n_frames, N)
y = zeros(n_frames, N)
z = zeros(n_frames, N)

open(traj_file) do file
    for (t, line) in enumerate(eachline(file))
        m = mod(t, fh + N)
        frame = (t ÷ (fh + N)) + 1
        if m != 1 && m != 2 && frame <= n_frames
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

x_e2e = zeros(n_frames, n_chains)
y_e2e = zeros(n_frames, n_chains)
z_e2e = zeros(n_frames, n_chains)

for i in 1:n_chains
    b = (i - 1) * len_chn + 1 # First bead of the chain.
    c = i * len_chn # Last bead of the chain.
    for frame in 1:n_frames
        e2e = [x[frame, c] - x[frame, b], y[frame, c] - y[frame, b], z[frame, c] - z[frame, b]]
        mag = norm(e2e) # magnitude.
        e2e_normalized = e2e ./ mag # normalising.
        
        x_e2e[frame, i], y_e2e[frame, i], z_e2e[frame, i] = e2e_normalized
    end
end

l_acf = n_frames ÷ 2
e2e_acf = zeros(l_acf)
for i in 1:(n_frames - l_acf)
    l = 0
    for j in i + 1:(i + l_acf)
        l += 1
        for k in 1:n_chains
            e2e_acf[l] += (x_e2e[j, k] * x_e2e[i, k]) + (y_e2e[j, k] * y_e2e[i, k]) + (z_e2e[j, k] * z_e2e[i, k])
        end
    end         
end

e2e_acf /= (n_chains * (n_frames - l_acf))

time = collect(1:l_acf)
time *= dt

out_file = open(out_file, "w")
string_out = "t acf\n"
write(out_file, string_out)
for i in 1:length(e2e_acf)
    global string_out = string(time[i]) * " " * string(e2e_acf[i]) * "\n"
    write(out_file, string_out)
end

# Calculate the integral of the E2E ACF using the trapezoidal rule
e2e_acf_integral = cumsum(0.5 * (e2e_acf[1:end-1] + e2e_acf[2:end]) * dt)

# Find the index when the integral reaches 1/e or 0.368
threshold = 1 / exp(1)
index_relaxation = findfirst(e2e_acf_integral .>= threshold)

# Get the relaxation time from the time array
if index_relaxation !== nothing
    relaxation_time = time[index_relaxation]
    println("Relaxation time: ", relaxation_time)
else
    println("Relaxation time not found within the given time range.")
end

# Create the ACF plot
p = plot(time, e2e_acf,
    line=:solid,
    marker=:circle,
    markercolor=:blue,
    markeralpha=0.5,
    markersize=4,
    xlabel="Timesteps",
    ylabel="ACF",
    title="End-to-end Autocorrelation Function",
    legend=:none,
    grid=:off,
    framestyle=:box)

# Save the ACF plot
savefig(p, "E2E_ACF.png")