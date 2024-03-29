function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
    vz::AbstractVector{T}) where T
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

function init_mesh(boxsize::Vector{Float64}, N::Vector{Int64})
    egde_x = boxsize[1] / N[1]
    egde_y = boxsize[2] / N[2]
    egde_z = boxsize[3] / N[3]
    edge_size = [egde_x, egde_y, egde_z]
    N_total::Int64 = N[1] * N[2] * N[3]
    particle_grid = zeros(N[1], N[2], N[3])
    boxmatrix = [boxsize[1] 0 0; 0 boxsize[2] 0; 0 0 boxsize[3]]

    ks = FFTW.fftfreq(N[1], 2π / (boxsize[1] / N[1]))
    ls = FFTW.fftfreq(N[2], 2π / (boxsize[2] / N[2])) # N[1]
    zs = FFTW.fftfreq(N[3], 2π / (boxsize[3] / N[3])) # N[1]

    Ks, Ls, Zs = meshgrid(ks, ls, zs)
    #print(Ks)
    
    Kmag = (Ks.^2 .+ Ls.^2 .+ Zs.^2)

    return Mesh(edge_size, N, particle_grid, boxsize, boxmatrix, Kmag, Ks, Ls, Zs)
end

function clear_mesh!(mesh::Mesh)
    for k = 1:mesh.N[1]
        for j = 1:mesh.N[2]
            for i = 1:mesh.N[3]
                mesh.particle_grid[k, j, i] = 0.0
            end
        end
    end
end

function getCellIndex(ix::Int64, iy::Int64, iz::Int64, N::Vector{Int64})
    cellindex::Int64 = ix + (iy - 1) * N[1] + (iz - 1) * N[1] * N[2]
    return cellindex
end

function getiXYZfromCellIndex(ixyz::Vector{Int64}, cellindex::Int64, N::Vector{Int64})
    ixyz[3] = floor((cellindex - 1) / (N[1] * N[2]))
    ixyz[2] = floor(((cellindex - 1) - N[1] * N[2] * ixyz[3]) / N[1])
    ixyz[1] = ((cellindex - 1) - N[1] * N[2] * ixyz[3]) % N[1]
    return ixyz
end

function pbc_mesh(x_index::Int64, N_mesh::Int64)
    x_index -= N_mesh * floor(x_index / N_mesh)
    x_index = Int(x_index)
    return x_index
end

function pbc_particle(position::Vector{Float64}, mesh::Mesh)
    position .-= mesh.boxsize .* floor.( position ./ mesh.boxsize)
    return position
end

function pbc_particle(position::Float64, boxsize::Float64)
    position -= boxsize * floor( position / boxsize)
    return position
end

function cloudincell(wrapped_pos::Array{Float64,1}, gridedge::Array{Float64,1}, meshnumber::Array{Int64,1})
    
    δx_lower = wrapped_pos[1] - floor(wrapped_pos[1] / edge_size[1]) * edge_size[1]
    δy_lower = wrapped_pos[2] - floor(wrapped_pos[2] / edge_size[2]) * edge_size[2]
    δz_lower = wrapped_pos[3] - floor(wrapped_pos[3] / edge_size[3]) * edge_size[3]

    δx_upper = edge_size[1] - δx_lower
    δy_upper = edge_size[1] - δy_lower
    δz_upper = edge_size[1] - δz_lower

end

function cloudincell!(position::Vector{Float64}, mesh::Mesh, i)
    # cell index = 1 + floor(x/cell_sizex) + floor(y/cell_sizey)*Nx + floor(z/cell_sizez)*Nx*Ny 

    # 3           4
    # -------------
    # |  |  |  |  | 
    # |  |  |  |  | 
    # |  |  |  |  |
    # -------------
    # 1           2 
    #if i < 15
    #    println("Nr: ", i)
    #    println("Position before pbc: ", position)
    #end

    #position_x = pbc_particle(position[1], mesh.boxsize[1])
    #position_y = pbc_particle(position[2], mesh.boxsize[2])
    #position_z = pbc_particle(position[3], mesh.boxsize[3])

    position_x = position[1] + 0.5 * mesh.boxsize[1]
    position_y = position[2] + 0.5 * mesh.boxsize[2]
    position_z = position[3] + 0.5 * mesh.boxsize[3]
    
    #if i < 15
    #    position = [position_x, position_y, position_z]
    #    println("Position after pbc: ", position)
    #end

    δx = position[1] - floor(position[1] / mesh.edge_size[1]) * mesh.edge_size[1]
    δy = position[2] - floor(position[2] / mesh.edge_size[2]) * mesh.edge_size[2]
    δz = position[3] - floor(position[3] / mesh.edge_size[3]) * mesh.edge_size[3]

    cell_index_x = floor(Int64, position_x / mesh.edge_size[1]) + 1
    cell_index_y = floor(Int64, position_y / mesh.edge_size[2]) + 1
    cell_index_z = floor(Int64, position_z / mesh.edge_size[3]) + 1

    cell_index_x_plus = pbc_mesh(cell_index_x, mesh.N[1]) + 1
    cell_index_y_plus = pbc_mesh(cell_index_y, mesh.N[2]) + 1
    cell_index_z_plus = pbc_mesh(cell_index_z, mesh.N[3]) + 1


    #if i < 15
    #    println("Cell_index x: ", cell_index_x)
    #    println("Cell_index y: ", cell_index_y)
    #    println("Cell_index z: ", cell_index_z)
    #end 

    mesh.particle_grid[cell_index_x, cell_index_y, cell_index_z]                += (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
    mesh.particle_grid[cell_index_x, cell_index_y_plus, cell_index_z]           += (mesh.edge_size[1] - δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                 
    mesh.particle_grid[cell_index_x, cell_index_y_plus, cell_index_z_plus]      += (mesh.edge_size[1] - δx) * (δy) * (δz) / prod(mesh.edge_size)
    mesh.particle_grid[cell_index_x, cell_index_y, cell_index_z_plus]           += (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)                 
    mesh.particle_grid[cell_index_x_plus, cell_index_y, cell_index_z]           += (δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                  
    mesh.particle_grid[cell_index_x_plus, cell_index_y, cell_index_z_plus]      += (δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
    mesh.particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z]      += (δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
    mesh.particle_grid[cell_index_x_plus, cell_index_y_plus, cell_index_z_plus] += (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      
    
end

function Grad_DensVertex(mesh::Mesh)
    
    grid_grad_x = zeros(mesh.N[1], mesh.N[2], mesh.N[3])
    grid_grad_y = zeros(mesh.N[1], mesh.N[2], mesh.N[3])
    grid_grad_z = zeros(mesh.N[1], mesh.N[2], mesh.N[3])

    for kcell = 1:mesh.N[3]
        for jcell = 1:mesh.N[2]
            for icell = 1:mesh.N[1]


                icell_plus = pbc_mesh(icell, mesh.N[1]) + 1
                jcell_plus = pbc_mesh(jcell, mesh.N[2]) + 1
                kcell_plus = pbc_mesh(kcell, mesh.N[3]) + 1

                icell_minus = pbc_mesh(icell - 2, mesh.N[1]) + 1
                jcell_minus = pbc_mesh(jcell - 2, mesh.N[2]) + 1
                kcell_minus = pbc_mesh(kcell - 2, mesh.N[3]) + 1
                
                grad_x = 0.5 * (mesh.particle_grid[icell_plus, jcell, kcell] - mesh.particle_grid[icell_minus, jcell, kcell])
                grad_y = 0.5 * (mesh.particle_grid[icell, jcell_plus, kcell] - mesh.particle_grid[icell, jcell_minus, kcell])
                grad_z = 0.5 * (mesh.particle_grid[icell, jcell, kcell_plus] - mesh.particle_grid[icell, jcell, kcell_minus])
                grid_grad_x[icell, jcell, kcell] += grad_x
                grid_grad_y[icell, jcell, kcell] += grad_y
                grid_grad_z[icell, jcell, kcell] += grad_z
            end
        end
    end
    return grid_grad_x, grid_grad_y, grid_grad_z
end
