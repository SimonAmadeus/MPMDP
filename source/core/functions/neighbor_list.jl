export standard_NL

struct standard_NL <: NL end

function create_neighbor_list!(args::System, particle_vec::Vector{Particle}, neighbor_list::Vector{Any}, ::No_NL) end

function create_neighbor_list!(args::System, particle_vec::Vector{Particle}, neighbor_list::Vector{Any}, NL::standard_NL)
    # First clear neighbor list, then recalculate.
    for i in 1:length(neighbor_list)
        neighbor_list[i] = []
    end
    for i = 2:length(neighbor_list)
        for j = 1: i - 1
            ∆d = particle_vec[i].pos .- particle_vec[j].pos
            ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
            r = norm(∆d)
            r_cutoff = 3
            if r < r_cutoff
                push!(neighbor_list[i], j)
                push!(neighbor_list[j], i)
            end
        end
    end
end