export DPD

struct DPD <: Non_Bonded_Interactions
    a::Float64 # Repulstion parameter.
    γ::Float64 # Friction parametter.
end

function apply_non_bonded_interactions!(args::System, ::System_Type, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, c_l::Int64, force_vec::Vector{Force}, energy::Float64, stress::Vector{Float64}, mesh::Mesh, DPD::DPD)
    σ = sqrt(2 * args.T * DPD.γ)
    for i = 2:length(particle_vec)
        for j = 1:i-1
            ∆d = particle_vec[i].pos .- particle_vec[j].pos
            ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
            r = norm(∆d)
            r_c = 1
            if r < r_c

                ∆v = velocity_vec[i].velocity .- velocity_vec[j].velocity
                ∆v_ij = ∆v[1]*∆d[1] + ∆v[2]*∆d[2] + ∆v[3]*∆d[3]
                ω = 1 - r / r_c 
                f_c = DPD.a * ω # Conservative force.
                f_d = - DPD.γ * ω^2 * ∆v_ij / r # Dissipative force.
                f_r = (2 * rand() - 1) * σ * ω / sqrt(args.∆t) # Random force.

                f = f_c + f_d + f_r # Total force.

                force_vec[i].force[1] += ∆d[1] * f / (r^2)
                force_vec[j].force[1] -= ∆d[1] * f / (r^2)
                force_vec[i].force[2] += ∆d[2] * f / (r^2)
                force_vec[j].force[2] -= ∆d[2] * f / (r^2)
                force_vec[i].force[3] += ∆d[3] * f / (r^2)
                force_vec[j].force[3] -= ∆d[3] * f / (r^2)
                energy += 0.5 * DPD.a *(1 - r / r_c)^2
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