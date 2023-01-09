import JSON

function read_scan(jsonfile::String)::Tuple{Vector{Int64}, Vector{Float64}, Array{Int16, 3}}

    # Read the json file
    info = JSON.parsefile(jsonfile)
    @assert info["dtype"]=="<i2" "Unsupported data type. Must be 2 byte little endian (<i2)."

    # Append a singleton dimension for two-dimensional data
    shape   = info["shape"]
    spacing = info["spacing"]
    if length(shape)==2
        push!(shape,1)
        push!(spacing,0.)
    end

    # Read the data
    intensity = Array{Int16}(undef, shape...)
    if info["order"]=="F" # first index changes fastest
        read!("data/" * info["fname"], intensity)
    elseif info["order"]=="C" # last index changes fastest
        Cordered = Array{Int16}(undef, reverse(shape)...)
        read!("data/" * info["fname"], Cordered)
        for i in CartesianIndices(intensity)
            intensity[i] = Cordered[i[3],i[2],i[1]]
        end
    end

    shape, spacing, intensity
end

# Parsing the cli
using ArgParse
function read_cli()::Tuple{String,Int64,Tuple{Int64,Int64,Int64},Int64}

    s = ArgParseSettings()
    @add_arg_table s begin
        "--order", "-p"
            help = "cardinal B-spline degree"
            arg_type = Int64
            default = 3
        "--depth", "-d"
            help = "sub-cell refinement depth"
            arg_type = Int64
            default = 2
        "--elements", "-n"
            help = "number of elements per direction (comma separated)"
            arg_type = Union{Int64,Vector{Int64}}
            default = [10,10,10]
        "name"
            help = "scan data json file"
            arg_type = String
            default = "synthetic_50x50.json"
    end

    parsed_args = parse_args(s)

    parsed_args["name"], parsed_args["order"], Tuple(parsed_args["elements"]), parsed_args["depth"]
end

function ArgParse.parse_item(::Type{Union{Int64,Vector{Int64}}}, x::AbstractString)
    array = [parse(Int64, String(strip(s))) for s ∈ split(x,',')]
    if length(array)==1
        return repeat(array, 3)
    elseif length(array)==2
        push!(array,1)
    end
    return array
end

# VTK unstructured grid plot
function write_vtk(name::String, cells::Vector{VTKCellWithCoordinates})

    ncells    = length(cells)
    vtkcells  = Array{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, Vector{Int64}}}(undef, ncells)
    vtkpoints = Array{Array{Float64}}(undef, ncells)

    nnodes = 0
    for (cid, cell) in enumerate(cells)
        vtkcells[cid]  = WriteVTK.MeshCell(cell.vtktype, nnodes.+collect(1:cell.vtktype.nodes))
        vtkpoints[cid] = cell.coords

        nnodes += cell.vtktype.nodes
    end

    vtk = WriteVTK.vtk_grid(name, hcat(vtkpoints...), vtkcells)
    WriteVTK.vtk_save(vtk)

    return vtk
end