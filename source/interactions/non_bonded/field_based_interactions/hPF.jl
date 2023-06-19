export original_hPF_Interactions, spectral_hPF_Interactions

struct original_hPF_Interactions <: Non_Bonded_Interactions 
    κ::Float64
    mesh_points::Vector{Int64}
end

struct spectral_hPF_Interactions <: Non_Bonded_Interactions 
    κ::Float64
    σ::Float64
    mesh_points::Vector{Int64}
end

function apply_non_bonded_interactions!(args::System, ::System_Type, particle_vec::Vector{Particle}, velocity_vec::Vector{Velocity}, neighbor_list::Any, c_l::Int64, force_vec::Vector{Force}, energy::Float64, stress::Vector{Float64}, mesh::Mesh, hPF::Non_Bonded_Interactions)

    clear_mesh!(mesh)
    for i = 1:length(particle_vec)
        cloudincell!(particle_vec[i].pos, mesh, i)
    end

    Particle_Field_Interactions!(particle_vec, force_vec, energy, mesh, hPF)

    return energy, stress
end

function Particle_Field_Interactions!(particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy::Float64, mesh::Mesh, hPF::original_hPF_Interactions)
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

        energy += mesh.particle_grid[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * 1 / (avg_den * hPF.κ)
        energy += mesh.particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * 1 / (avg_den * hPF.κ)

        energy -= 1 / hPF.κ

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

    return energy

end

function Particle_Field_Interactions!(particle_vec::Vector{Particle}, force_vec::Vector{Force}, energy::Float64, mesh::Mesh, hPF::spectral_hPF_Interactions)

    avg_den = length(particle_vec) / prod(mesh.boxsize)
    pic=zeros(8)

    # apply gaussian filter

    # forward transform
    complex_field = fft(mesh.particle_grid)

    #convolution with gaussian
    conv_field = complex_field .* exp.(-0.5 .* hPF.σ^2 .* mesh.Kmag)

    # backward transform
    particle_grid = real(ifft(conv_field))

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

        energy += particle_grid[cell_index_x, cell_index_y, cell_index_z]                * pic[1] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x, cell_index_y_plus, cell_index_z]           * pic[2] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x, cell_index_y_plus, cell_index_z_plus]      * pic[3] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x, cell_index_y, cell_index_z_plus]           * pic[4] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x_plus, cell_index_y, cell_index_z]           * pic[5] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x_plus, cell_index_y, cell_index_z_plus]      * pic[6] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z]      * pic[7] * 1 / (avg_den * hPF.κ)
        energy += particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] * pic[8] * 1 / (avg_den * hPF.κ)
        energy -= 1 / hPF.κ

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

    return energy

end
