export LJ_12_6

struct LJ_12_6 <: Non_Bonded_Interactions
    ε::Float64
    σ::Float64
end

# Lennard-Jones interactions for monoatomic systems.
function apply_non_bonded_interactions!(args::System, ::No_Mol, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, c_l::Int64, force_vec::Vector{Force}, energy::Float64, stress::Vector{Float64}, mesh::Mesh, lj126::LJ_12_6)

    for i = 2:length(particle_vec)
        for j = 1:i-1
            ∆d = particle_vec[i].pos .- particle_vec[j].pos
            ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
            r = norm(∆d)
            r6i = 0
            #r_cutoff = norm(args.box) / 2
            r_cutoff = 3
            if r < r_cutoff
                r6i = 1.0 / r^6
                f = 48 * lj126.ε * (lj126.σ * r6i^2 - 0.5 * lj126.σ * r6i)
                force_vec[i].force[1] += ∆d[1] * f / (r^2)
                force_vec[j].force[1] -= ∆d[1] * f / (r^2)
                force_vec[i].force[2] += ∆d[2] * f / (r^2)
                force_vec[j].force[2] -= ∆d[2] * f / (r^2)
                force_vec[i].force[3] += ∆d[3] * f / (r^2)
                force_vec[j].force[3] -= ∆d[3] * f / (r^2)
                energy += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i)
                # xx, yy, zz, xy, xz, yz
                stress[1] += force_vec[i].force[1] * ∆d[1]
                stress[2] += force_vec[i].force[2] * ∆d[2]
                stress[3] += force_vec[i].force[3] * ∆d[3]
                stress[4] += force_vec[i].force[1] * ∆d[2]
                stress[5] += force_vec[i].force[1] * ∆d[3]
                stress[6] += force_vec[i].force[2] * ∆d[3]
            end
        end
    end
    return energy, stress
end

# 1-2 and 1-3 exclusion for bonded atoms.
function apply_non_bonded_interactions!(args::System, ::Molecular, particle_vec::Vector{Particle}, c_l::Int64, force_vec::Vector{Force}, energy::Float64, stress::Vector{Float64}, mesh::Mesh, lj126::LJ_12_6)

    for i = 2:length(particle_vec)
        for j = 1:i-1
            if abs(i - j) > 2 || ( mod(i, c_l) == 1 && (j - i == -1 || j - i == - 2)) || ( mod(j, c_l) == 1 && (i - j == -1 || i - j == - 2)) || ( mod(i, c_l) == 2 && j - i == - 2 ) || ( mod(j, c_l) == 2 && i - j == - 2 ) || ( mod(i, c_l) == 0 && (j - i == 1 || j - i == 2)) || ( mod(j, c_l) == 0 && (i - j == 1 || i - j == 2)) || ( mod(i, c_l) == c_l - 1 && j - i == 2) || ( mod(j, c_l) == c_l - 1 && i - j == 2)
                ∆d = particle_vec[i].pos .- particle_vec[j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                r6i = 0
                #r_cutoff = norm(args.box) / 2
                r_cutoff = 3
                if r < r_cutoff
                    r6i = 1.0 / r^6
                    f = 48 * lj126.ε * (lj126.σ * r6i^2 - 0.5 * lj126.σ * r6i)
                    force_vec[i].force[1] += ∆d[1] * f / (r^2)
                    force_vec[j].force[1] -= ∆d[1] * f / (r^2)
                    force_vec[i].force[2] += ∆d[2] * f / (r^2)
                    force_vec[j].force[2] -= ∆d[2] * f / (r^2)
                    force_vec[i].force[3] += ∆d[3] * f / (r^2)
                    force_vec[j].force[3] -= ∆d[3] * f / (r^2)
                    energy += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i)
                    # xx, yy, zz, xy, xz, yz
                    stress[1] += force_vec[i].force[1] * ∆d[1]
                    stress[2] += force_vec[i].force[2] * ∆d[2]
                    stress[3] += force_vec[i].force[3] * ∆d[3]
                    stress[4] += force_vec[i].force[1] * ∆d[2]
                    stress[5] += force_vec[i].force[1] * ∆d[3]
                    stress[6] += force_vec[i].force[2] * ∆d[3]
                end
            end
        end
    end
    return energy, stress
end
