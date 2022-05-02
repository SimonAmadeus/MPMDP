using LinearAlgebra
using ProgressMeter

mutable struct Positions
    pos::Vector{Float64}
end

#println("Enter input file (XYZ format).")
#traj_file = chomp(readline())
traj_file = "trj.xyz"
#println("Enter out-file (e.g. rdf.txt).")
#out_file = chomp(readline())
out_file = "rdf.txt"
#println("How many frames would you like to analyze?")
#n_frames = parse(Int64, chomp(readline()))
n_frames = 100

counter = 0
fh = 2 # File header.

N = 0 # Number of particles.
box = [10, 10, 30]

n_bins = 100 # Number of bins.
bins = zeros(Float64, 100) # Storage for particle numbers.
d_bins = 5 # Distrance of RDF.
s_bin = d_bins / n_bins # Size of a bin.

open(traj_file) do file
    for line in eachline(file)
        global counter += 1
        a = split(line)
        if counter == 1
            global N = parse(Int64, a[1])
            println("$N particles in the system.")
        end
        break
    end
end

ρ = N / (box[1] * box[2] * box[3])
coords = Vector{Positions}(undef, N)

counter = 0
counter2 = 0 # Within frame -> particle counter.

@showprogress for i in 1:n_frames
    
    # Read data, 1 frame at once.
    open(traj_file) do file
        counter = 0
        counter2 = 0
        for line in eachline(file)
            global counter += 1
            a = split(line)
            start_i = fh + (N + fh) * (i - 1)
            end_i = (N + fh) + (N + fh) * (i - 1)
            if counter > start_i && counter < (end_i + 1)
                counter2 += 1
                pos = [parse(Float64, a[2]), parse(Float64, a[3]), parse(Float64, a[4])]
                pos -= box .* floor.((pos .+ box / 2) ./ box)
                global pos_i = Positions(pos)
                coords[counter2] = pos_i
            end
        end
    end
    
    # RDF loop for every frame.
    for i in 2:N
        for j in 1:i-1
            ∆d = coords[i].pos .- coords[j].pos
            ∆d .-= box .* floor.((∆d .+ box / 2) ./ box)
            r = norm(∆d)
            ii = Int(floor(r / s_bin) )
            if ii <= n_bins - 1
                bins[ii + 1] += 2 / (n_frames * N)
            end
        end
    end

end
       
for i in 1:n_bins
    v_bin = ((i + 1)^3 - i^3)*s_bin^3 # Volume of a bin.
    nid = (4 / 3) * π * v_bin * ρ # Average number of particles per volume.
    bins[i] = bins[i] / nid # Normalize.
end

out_file = open(out_file, "w")

string_out = "r g\n"
write(out_file, string_out)

for i in 1:n_bins
    global string_out = string(i * s_bin) * " " * string(bins[i]) * "\n"
    write(out_file, string_out)
end
