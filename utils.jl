using StaticArrays: SVector

# Local cardinal B-spline filtering stencil
"""
    cardinal_stencil(p)

Compute the local filtering stencil of the cardinal B-spline.

# Arguments
- `p::Int64`: B-spline polynomial order.
"""
function cardinal_stencil(p::Int64)::Vector{Float64}
    n = p+1
    s = zeros(n)
    for i=0:n-1
        for k=0:n
            s[i+1]+=((-1)^k)*binomial(n,k)*(i<k ? 0 : (i+1-k)^n-(i-k)^n)
        end
    end
    s ./= factorial(n)
end

# Local cardinal B-spline basis
"""
    cardinal_basis(p, η)

Compute the `p+1` basis functions at the local coordinates `η∈[0,1]`.

# Arguments
- `p::Int64`: B-spline polynomial order.
- `η::Vector{Float64}`: array of local coordinates.
"""
function cardinal_basis(p::Int64, η::Vector{Float64})::Matrix{Float64}

    if p==0
        return ones(1, length(η))
    end

    ξ = Vector(p:-1:0) .+ reverse(η)'
    B = zeros(Float64, size(ξ))
    n = p+1
    for k=0:n
        ξ_minus_k = (ξ.-k).*(ξ.>k)
        B+=((-1)^k)*binomial(n,k)*(ξ_minus_k).^p
    end
    B ./= factorial(p)
end

# Get the global cardianal B-spline indices
function cardinal_indices(i::Int64, p::Integer)::UnitRange{Int64}
    i:i+p
end

# Octree subdivision criterion
function RegionTrees.needs_refinement(r::LevelsetRefinery, cell)
    if size(cell.data,1)==2 || all(cell.data.>0.) || all(cell.data.<=0.)
        return false
    end
    size(cell.data,1)>2
end

# Octree subdivision data split
function RegionTrees.refine_data(r::LevelsetRefinery, cell::Cell, indices)
    npoints_child = div(size(cell.data,1)-1, 2)+1
    if length(indices)==2
        cell.data[((indices[d]-1)*(npoints_child-1)+1:indices[d]*(npoints_child-1)+1 for d in 1:2)...,:]
    elseif length(indices)==3
        cell.data[((indices[d]-1)*(npoints_child-1)+1:indices[d]*(npoints_child-1)+1 for d in 1:3)...]
    end
end

# Mosaic function for lowest level of subdivision
function get_mosaic(cell::Cell)::Vector{VTKCellWithCoordinates}

    ndims = length(cell.boundary.origin)
    # println("Call to get_mosaic with ndims=", ndims)

    shape  = size(cell.data)
    coords = ndims>0 ? hcat(collect(vertices(cell.boundary))...) : Array{Float64}(undef, 0, 1)
    data   = collect(Iterators.flatten(cell.data))

    @assert size(coords)==(ndims,2^ndims)
    @assert size(data)==(2^ndims,)

    # All level set values are positive
    if all(cell.data.>0.) # Hyperrectangle of dimension `ndims`
        # println("All positive")
        return [VTKCellWithCoordinates(vtktypes[(ndims, 2^ndims)], coords)]

    # All level set values are non-positive
    elseif all(cell.data.<=0.)
        # println("All non-positive")
        return []

    # Both positive and non-positive level set values
    else
        # println("Both positive and negative")

        # Compute the value in the centroid based on linear interpolation
        ccoord = reshape(sum(Matrix(coords),dims=2),ndims)./size(coords,2)
        cdata  = sum(data)/length(data)

        # Compute the midpoint
        if cdata≈0.
            midpoint = ccoord
        else
            number_of_zeros = 0
            midpoint = zeros(size(ccoord))
            for (vcoord, vdata) in zip(eachcol(coords), data)
                # Zero point when the level set changes sign
                if cdata*vdata<=0.
                    ξ = cdata/(cdata-vdata)
                    midpoint += ccoord.+ξ*(vcoord.-ccoord)
                    number_of_zeros+=1
                end
            end
            midpoint ./= number_of_zeros
            # println("midpoint = ", midpoint)
        end

        # Loop over the faces (edges in 2D) of the cell
        mosaic = Vector{VTKCellWithCoordinates}(undef, 0)
        for offset_direction in reverse(1:ndims)
            for offset_index in 1:2
                bindices = boundary_indices[(ndims,offset_direction,offset_index)]
                bscale   = Array{Float64}(reshape(hcat(collect([((j>=offset_direction ? i-1 : i)==j) for i=1:ndims] for j=1:ndims-1)...), ndims, ndims-1))
                boffset  = Array{Float64}([((i==offset_direction)&&(offset_index==2)) for i=1:ndims])

                bdata   = reshape(data[bindices], shape[1:ndims-1])

                # Create the lower-dimensional unit cell
                bcell = Cell(SVector(zeros(ndims-1)...), SVector(ones(ndims-1)...), bdata)

                # Get the mosaic of the lower-dimensional unit cell
                bmosaic = get_mosaic(bcell)

                # Extrude boundaries toward the midpoint
                for btile in bmosaic
                    btile_coords = (bscale*btile.coords.+boffset).*cell.boundary.widths.+cell.boundary.origin
                    if btile.vtktype==vtktypes[(0,1)] # Vertex to line segment
                        # println("vertex -> line")
                        push!(mosaic, VTKCellWithCoordinates(vtktypes[(1,2)], hcat(btile_coords, midpoint)))
                    elseif btile.vtktype==vtktypes[(1,2)] # Line segment to triangle
                        # println("line -> triangle")
                        push!(mosaic, VTKCellWithCoordinates(vtktypes[(2,3)], hcat(btile_coords, midpoint)))
                    elseif btile.vtktype==vtktypes[(2,4)] # Quadrilateral to pyramid
                        push!(mosaic, VTKCellWithCoordinates(vtktypes[(3,5)], hcat(btile_coords, midpoint)))
                    elseif btile.vtktype==vtktypes[(2,3)] # Triangle to tetrahedron
                        push!(mosaic, VTKCellWithCoordinates(vtktypes[(3,4)], hcat(btile_coords, midpoint)))
                    else
                        error("Extrusion not implemented for VTK type ", btile.vtktype)
                    end
                end

            end
        end

        return mosaic
    end
end