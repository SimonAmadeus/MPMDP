export Harmonic_Bonded, FENE_Bonded, LJ_12_6, original_hPF_Interactions, spectral_hPF_Interactions

function clear_force!(force_vec::Vector{Force})
    for i = 1:length(force_vec)
        force_vec[i].force[1] = 0.0
        force_vec[i].force[2] = 0.0
        force_vec[i].force[3] = 0.0
    end
end

function clear_energy!(energy_vec::Vector{Float64})
    for i = 1:length(energy_vec)
        energy_vec[i] = 0.0
    end
end


# Bonded interactions.
function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, bond_vec::Vector{Bond}, ::No_Bonds) end

struct Harmonic_Bonded <: Bonded_Interactions 
    κ::Float64
    r0::Float64
end

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, bond_vec::Vector{Bond}, harmonic_bond::Harmonic_Bonded) 
    
    # 0.5×k×(r-r0)^2
    for i = 1:length(bond_vec)
        pos_1 = particle_vec[bond_vec[i].p_2].pos .+ particle_vec[bond_vec[i].p_2].image .* args.box
        pos_2 = particle_vec[bond_vec[i].p_1].pos .+ particle_vec[bond_vec[i].p_1].image .* args.box
        ∆d = pos_2 .- pos_1
        #∆d .-= args.box .* floor.((∆d + args.box / 2)./ args.box) # PBCs.
        r = norm(∆d)
        f = - harmonic_bond.κ * (harmonic_bond.r0 / r - 1) 
        force_divr = f / r
        force_vec[bond_vec[i].p_1].force .-= ∆d .* force_divr 
        force_vec[bond_vec[i].p_2].force .+= ∆d .* force_divr
        energy_vec[bond_vec[i].p_1] += 0.5 * 0.5 * harmonic_bond.κ * (r - harmonic_bond.r0)^2
        energy_vec[bond_vec[i].p_2] += 0.5 * 0.5 * harmonic_bond.κ * (r - harmonic_bond.r0)^2
    end

end

struct FENE_Bonded <: Bonded_Interactions
    K::Float64 #    (energy/distance^2)
    R0::Float64 #   (distance)
    ϵ::Float64 #    (energy)
    σ::Float64 #    (distance)
end

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, bond_vec::Vector{Bond}, fenebond::FENE_Bonded) 

    for i = 1:length(bond_vec)

        pos_1 = particle_vec[bond_vec[i].p_2].pos .+ particle_vec[bond_vec[i].p_2].image .* args.box
        pos_2 = particle_vec[bond_vec[i].p_1].pos .+ particle_vec[bond_vec[i].p_1].image .* args.box
        ∆d = pos_2 .- pos_1

        r = norm(∆d)
        rsq = r^2

        r0sq = fenebond.R0^2
        rlogarg = 1 - rsq / r0sq

        fbond = - fenebond.K / rlogarg

        if rsq < fenebond.σ^2
            sr2 = fenebond.σ^2 / rsq
            sr6 = sr2*sr2*sr2
            force_divr += 48*fenebond.ϵ * sr6 * (sr6 - 0.5) / rsq
        end

        E = -0.5*fenebond.K*fenebond.R0^2*log(1-(r/fenebond.R0)^2)+4*fenebond.ϵ*((fenebond.σ/r)^12-(fenebond.σ/r)^6)+fenebond.ϵ

        energy[bond_vec[i].p_1+1]+=E*0.5
        energy[bond_vec[i].p_2+1]+=E*0.5

        # calculate forces
        force_vec[bond_vec[i].p_1+1].force .+= delta_d .* force_divr
        force_vec[bond_vec[i].p_2+1].force .+= -delta_d .* force_divr
    end
end

# Non-bonded interactions.

function apply_non_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, ::No_Non_Bonded) end

struct LJ_12_6 <: Non_Bonded_Interactions
    ε::Float64
    σ::Float64
end

struct original_hPF_Interactions <: Non_Bonded_Interactions 
    mesh::Mesh
    κ::Float64
end

struct spectral_hPF_Interactions <: Non_Bonded_Interactions 
    mesh::Mesh
    κ::Float64
    σ::Float64
end

function apply_non_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, lj126::LJ_12_6)
    
    for i = 2:length(particle_vec)
        for j = 1:i-1
            ∆d = particle_vec[i].pos .- particle_vec[j].pos
            ∆d .-= args.box .* floor.((∆d .+ args.box / 2)./args.box)
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
                energy_vec[i] += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i) 
                energy_vec[j] += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i)
            end 
        end
    end
    
end

function Particle_Field_Interactions!(particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, mesh::Mesh, hPF::original_hPF_Interactions)
    avg_den = length(particle_vec) / prod(mesh.boxsize)
    pic = zeros(8)
    
    # force on the grids: -1 * derivative of the potential on grids
    grid_grad_x, grid_grad_y, grid_grad_z = Grad_DensVertex(mesh)

    for i = 1:length(particle_vec)

        ix::Int64 = floor(particle_vec[i].pos[1] / mesh.edge_size[1])
        iy::Int64 = floor(particle_vec[i].pos[2] / mesh.edge_size[2])
        iz::Int64 = floor(particle_vec[i].pos[3] / mesh.edge_size[3])

        δx = particle_vec[i].pos[1] - ix * mesh.edge_size[1]
        δy = particle_vec[i].pos[2] - iy * mesh.edge_size[2]
        δz = particle_vec[i].pos[3] - iz * mesh.edge_size[3]

        pic[1] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[2] = (mesh.edge_size[1] - δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                 
        pic[3] = (mesh.edge_size[1] - δx) * (δy) * (δz) / prod(mesh.edge_size)
        pic[4] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
        pic[5] = (δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                  
        pic[6] = (δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)                                   
        pic[7] = (δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                                   
        pic[8] = (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      

        icell = ix
        jcell = iy
        kcell = iz

        cell_index_x = pbc_mesh(icell, mesh.N[1])
        cell_index_y = pbc_mesh(jcell, mesh.N[2])
        cell_index_z = pbc_mesh(kcell, mesh.N[3])

        cell_index_x_plus = pbc_mesh(cell_index_x + 1, mesh.N[1])
        cell_index_y_plus = pbc_mesh(cell_index_y + 1, mesh.N[2])
        cell_index_z_plus = pbc_mesh(cell_index_z + 1, mesh.N[3])

        cell_index_x = cell_index_x + 1
        cell_index_y = cell_index_y + 1
        cell_index_z = cell_index_z + 1

        cell_index_x_plus = cell_index_x_plus + 1
        cell_index_y_plus = cell_index_y_plus + 1
        cell_index_z_plus = cell_index_z_plus + 1

        energy_vec[i] += mesh.grids[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += mesh.grids[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * 1 / (avg_den * hPF.κ)

        energy_vec[i] -= 1 / hPF.κ

        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * - 1 / hPF.κ / avg_den
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * - 1 / hPF.κ / avg_den

        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * - 1 / hPF.κ / avg_den
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * - 1 / hPF.κ / avg_den

        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * - 1 / hPF.κ / avg_den 
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * - 1 / hPF.κ / avg_den
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * - 1 / hPF.κ / avg_den

    end
   
end

function Particle_Field_Interactions!(particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, mesh::Mesh, hPF::spectral_hPF_Interactions)
    
    avg_den = length(particle_vec) / prod(mesh.boxsize)
    pic=zeros(8)

    # apply gaussian filter

    # forward transform
    complex_field = fft(mesh.grids) 
    
    #convolution with gaussian
    conv_field = complex_field .* exp.(-0.5 .* hPF.σ^2 .* mesh.Kmag)
    
    # backward transform
    grids = real(ifft(conv_field)) 
    
    # force on the grids: -1 * derivative of the potential on grids
    grid_grad_x = real(ifft( conv_field .* 1 / (avg_den*hPF.κ) .* -1 .* mesh.Ls .* 1im)) 
    grid_grad_y = real(ifft( conv_field .* 1 / (avg_den*hPF.κ) .* -1 .* mesh.Ks .* 1im)) 
    grid_grad_z = real(ifft( conv_field .* 1 / (avg_den*hPF.κ) .* -1 .* mesh.Zs .* 1im)) 

    for i = 1:length(particle_vec)

        ix::Int64 = floor(particle_vec[i].pos[1] / mesh.edge_size[1])
        iy::Int64 = floor(particle_vec[i].pos[2] / mesh.edge_size[2])
        iz::Int64 = floor(particle_vec[i].pos[3] / mesh.edge_size[3])

        δx = particle_vec[i].pos[1] - ix * mesh.edge_size[1]
        δy = particle_vec[i].pos[2] - iy * mesh.edge_size[2]
        δz = particle_vec[i].pos[3] - iz * mesh.edge_size[3]

        pic[1] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[2] = (mesh.edge_size[1] - δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                 
        pic[3] = (mesh.edge_size[1] - δx) * (δy) * (δz) / prod(mesh.edge_size)                                   
        pic[4] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)                 
        pic[5] = (δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                  
        pic[6] = (δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                                   
        pic[7] = (δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                                   
        pic[8] = (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      

        icell = ix
        jcell = iy
        kcell = iz

        cell_index_x = pbc_mesh(icell, mesh.N[1])
        cell_index_y = pbc_mesh(jcell, mesh.N[2])
        cell_index_z = pbc_mesh(kcell, mesh.N[3])

        cell_index_x_plus = pbc_mesh(cell_index_x + 1, mesh.N[1])
        cell_index_y_plus = pbc_mesh(cell_index_y + 1, mesh.N[2])
        cell_index_z_plus = pbc_mesh(cell_index_z + 1, mesh.N[3])

        cell_index_x = cell_index_x + 1
        cell_index_y = cell_index_y + 1
        cell_index_z = cell_index_z + 1

        cell_index_x_plus = cell_index_x_plus + 1
        cell_index_y_plus = cell_index_y_plus + 1
        cell_index_z_plus = cell_index_z_plus + 1

        energy_vec[i] += grids[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * 1 / (avg_den * hPF.κ)
        energy_vec[i] += grids[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * 1 / (avg_den * hPF.κ)
        energy_vec[i] -= 1 / hPF.κ

        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y, cell_index_z]                * pic[1]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7]
        force_vec[i].force[1] +=  grid_grad_x[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8]

        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y, cell_index_z]                * pic[1]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7]
        force_vec[i].force[2] +=  grid_grad_y[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8]

        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y, cell_index_z]                * pic[1]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7]
        force_vec[i].force[3] +=  grid_grad_z[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8]

    end
   
end


function apply_non_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy_vec::Vector{Float64}, hPF::Non_Bonded_Interactions)

    clear_mesh!(hPF.mesh)
    for i = 1:length(particle_vec)
        cloudincell!(particle_vec[i].pos, hPF.mesh)
    end
    Particle_Field_Interactions!(particle_vec, force_vec, energy_vec, hPF.mesh, hPF)
end