import WriteVTK
using RegionTrees

vtktypes = Dict((0,1)=>WriteVTK.VTKCellTypes.VTK_VERTEX,
                (1,2)=>WriteVTK.VTKCellTypes.VTK_LINE,
                (2,3)=>WriteVTK.VTKCellTypes.VTK_TRIANGLE,
                (2,4)=>WriteVTK.VTKCellTypes.VTK_QUAD,
                (3,4)=>WriteVTK.VTKCellTypes.VTK_TETRA,
                (3,5)=>WriteVTK.VTKCellTypes.VTK_PYRAMID,
                (3,8)=>WriteVTK.VTKCellTypes.VTK_HEXAHEDRON)

vtkordering = Dict(WriteVTK.VTKCellTypes.VTK_QUAD      =>[1,2,4,3],
                   WriteVTK.VTKCellTypes.VTK_HEXAHEDRON=>[1,2,4,3,5,6,8,7])

mutable struct VTKCellWithCoordinates
    vtktype::WriteVTK.VTKCellType
    coords::Array{Float64}

    VTKCellWithCoordinates(vtktype, coords) = new(vtktype, vtktype ∈ keys(vtkordering) ? coords[:,vtkordering[vtktype]] : coords)
end

vtktypes = Dict((0,1)=>WriteVTK.VTKCellTypes.VTK_VERTEX,
                (1,2)=>WriteVTK.VTKCellTypes.VTK_LINE,
                (2,3)=>WriteVTK.VTKCellTypes.VTK_TRIANGLE,
                (2,4)=>WriteVTK.VTKCellTypes.VTK_QUAD,
                (3,4)=>WriteVTK.VTKCellTypes.VTK_TETRA,
                (3,5)=>WriteVTK.VTKCellTypes.VTK_PYRAMID,
                (3,8)=>WriteVTK.VTKCellTypes.VTK_HEXAHEDRON)

boundary_indices = Dict()
for ndims in 1:3
    for offset_direction in 1:ndims
        for offset_index in 1:2
            boundary_indices[(ndims, offset_direction, offset_index)] = []
            for (i, index) in enumerate(CartesianIndices(Tuple(repeat([2],ndims))))
                if index[offset_direction]==offset_index
                    push!(boundary_indices[(ndims, offset_direction, offset_index)], i)
                end
            end
        end
    end
end

struct LevelsetRefinery <: AbstractRefinery
end