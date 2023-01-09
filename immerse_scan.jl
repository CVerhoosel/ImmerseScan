include("./defs.jl") # Constants and structs
include("./utils.jl") # Utility functions
include("./io.jl") # Input/output functions

function main(name, p, n, depth)

    # Read the scan data
    shape, spacing, voxels = read_scan("data/"*name)
    ndims = (shape[3]!=1 ? 3 : 2)

    # B-spline polynomial degrees per direction
    degrees  = Tuple(min(shape[dim]-1,p) for dim=1:3)
    n        = ndims==3 ? n : (n[1], n[2], 1)
    meshsize = (shape.*spacing./n)

    # Perform the convolution operation to obtain the cardinal B-spline coefficients
    stencils = Tuple(cardinal_stencil(degree) for degree in degrees)
    lcoeffs  = zeros((shape.+degrees)...)
    weights  = zeros((shape.+degrees)...)
    for v in CartesianIndices(voxels)
        lcoeffs_indices = Tuple(cardinal_indices(v[d], degrees[d]) for d in 1:3)
        for i in CartesianIndices(degrees.+1)
            lcoeffs_index = Tuple(lcoeffs_indices[d][i[d]] for d in 1:3)
            w = prod(stencils[d][i[d]] for d in 1:3)
            lcoeffs[lcoeffs_index...] += w*voxels[v]
            weights[lcoeffs_index...] += w
        end
    end
    lcoeffs = lcoeffs./weights

    # Evaluate the levelset function
    levelset = zeros((shape.+(shape.!=1))...)
    B = Tuple(cardinal_basis(degree, [0.,1.]) for degree in degrees)
    for v in CartesianIndices(voxels)

        lc = lcoeffs[Tuple(cardinal_indices(v[d], degrees[d]) for d in 1:3)...]
        levelset_indices = Tuple(cardinal_indices(v[d], shape[d]!=1) for d in 1:3)

        cell_levelset = zeros(((shape.!=1).+1)...)
        for c in CartesianIndices(cell_levelset)
            for i in CartesianIndices(size(lc))
                cell_levelset[c] += B[1][i[1],c[1]]*B[2][i[2],c[2]]*B[3][i[3],c[3]]*lc[i]
            end
        end
        levelset[levelset_indices...] = cell_levelset
    end

    # Write the image to a vti file
    vtk = WriteVTK.vtk_grid("voxeldata", (shape.+(shape.!=1))..., spacing=(spacing...,))
    vtk["scan"] = voxels
    vtk["levelset"] = levelset
    WriteVTK.vtk_save(vtk)

    ############
    # Trimming #
    ############

    mesh = Array{Union{Cell{Array{Float64, 3}, ndims, Float64, 2^ndims}, Nothing}}(undef, n...)
    npoints = 2^depth+1
    unit_points = Vector(LinRange(0.0, 1.0, npoints))
    V = Array{Float64}(undef, (size(mesh,d)!=1 ? n[d]*(npoints-1)+1 : 1 for d in 1:3)...)
    for c in CartesianIndices(size(mesh))
        cell_levelset = zeros(((size(mesh,d)!=1 ? npoints : 1) for d in 1:3)...)
        for p in CartesianIndices(cell_levelset)
            x  = [(c[d]+unit_points[p[d]]-1)*meshsize[d] for d in 1:3]
            v  = [(size(mesh,d)!=1 ? min(max(ceil(Int64, x[d]/spacing[d]),1),size(voxels,d)) : 1) for d in 1:3]
            ξ  = [(x[d]-(v[d]-1)*spacing[d])/spacing[d] for d in 1:3]
            lc = lcoeffs[Tuple(cardinal_indices(v[d], degrees[d]) for d in 1:3)...]
            Φ = Tuple(cardinal_basis(degrees[d], [ξ[d]]) for d in 1:3)
            for i in CartesianIndices(size(lc))
                cell_levelset[p] += Φ[1][i[1],1]*Φ[2][i[2],1]*Φ[3][i[3],1]*lc[i]
            end
        end

        V[(size(mesh,d)!=1 ? ((c[d]-1)*(npoints-1)+1:c[d]*(npoints-1)+1) : 1 for d in 1:3)...] = cell_levelset

        root = Cell(SVector(((Tuple(c).-1).*meshsize)[1:ndims]...),
                    SVector(meshsize[1:ndims]...), cell_levelset)

        if !all(cell_levelset.>0.)
            adaptivesampling!(root, LevelsetRefinery())
        end
        mesh[c] = root

    end

    # Write the levelset data to a vti file
    vtk = WriteVTK.vtk_grid("levelset", size(V)..., spacing=Tuple(meshsize./(npoints-1)))
    vtk["levelset"] = V
    WriteVTK.vtk_save(vtk)

    # Gather all subcells
    vtkcells = Vector{VTKCellWithCoordinates}(undef, 0)
    for element in filter(!isnothing, mesh)
        for subcell in allleaves(element)
            if all(subcell.data.<=0.)
                continue
            elseif all(subcell.data.>0.)
                push!(vtkcells, VTKCellWithCoordinates(vtktypes[(ndims,2^ndims)], hcat(collect(vertices(subcell.boundary))...)))
            else
                push!(vtkcells, get_mosaic(subcell)...)
            end
        end
    end

    vtk = write_vtk("trimmed", vtkcells)
end

# Read command line and callt the main script
main(read_cli()...)