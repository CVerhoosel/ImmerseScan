# Unit tests
using Test

include("./immerse_scan.jl")

@testset verbose = true "Mosaic" begin

    # Line segment from -1 to 3
    #
    # ( 1 )---(-1 )
    mosaic = get_mosaic( Cell(SVector(-1.), SVector(4.), [1., -1.]) )
    @test hcat([tile.coords for tile in mosaic]...)≈[-1 1]

    # (-1 )---( 1 )
    mosaic = get_mosaic( Cell(SVector(-1.), SVector(4.), [-1., 1.]) )
    @test hcat([tile.coords for tile in mosaic]...)≈[3 1]

    # ( 3 )---( -1 )
    mosaic = get_mosaic( Cell(SVector(-1.), SVector(4.), [3., -1.]) )
    @test hcat([tile.coords for tile in mosaic]...)≈[-1 2]

    # Unit square with:
    #
    # j ^
    #   |
    # ( 1 )---(-1 )
    #   |       |
    # ( 3 )---( 1 )-> i
    mosaic = get_mosaic( Cell(SVector(0., 0.), SVector(1., 1.), [3. 1. ; 1. -1.]) )
    result = [0  1  3/4  0  1/2  3/4  0  0  3/4  1  1    3/4
              0  0  3/4  1  1    3/4  0  1  3/4  0  1/2  3/4]
    @test hcat([tile.coords for tile in mosaic]...)≈result

    # ( 1 )---(-1 )
    #   |       |
    # ( 1 )---( 1 )
    mosaic = get_mosaic( Cell(SVector(0., 0.), SVector(1., 1.), [1. 1. ; 1. -1.]) )
    result = [0  1  2/3  0  1/2  2/3  0  0  2/3  1  1    2/3
              0  0  2/3  1  1    2/3  0  1  2/3  0  1/2  2/3]
    @test hcat([tile.coords for tile in mosaic]...)≈result

    # ( 1 )---(-1 )
    #   |       |
    # ( 1 )---(-1 )
    mosaic = get_mosaic( Cell(SVector(0., 0.), SVector(1., 1.), [1. 1. ; -1. -1.]) )
    result = [ 0  1/2  1/2  0  1/2  1/2  0  0  1/2
               0  0    1/2  1  1    1/2  0  1  1/2]
    @test hcat([tile.coords for tile in mosaic]...)≈result

    # 3×2 Rectanlge with origin at (-2,4):
    #
    # ( 1 )---(-1 )
    #   |       |
    # ( 3 )---( 1 )
    mosaic = get_mosaic( Cell(SVector(-2., 4.), SVector(3., 2.), [3. 1. ; 1. -1.]) )
    result = [-2  1  1/4   -2  -1/2  1/4   -2  -2  1/4   1  1  1/4
               4  4  11/2   6   6    11/2   4   6  11/2  4  5  11/2]
    @test hcat([tile.coords for tile in mosaic]...)≈result

    # Unit cube with:
    #
    #    j ^
    #      |
    #    ( 3 )----( 1 )
    #     /        /|
    #    / |      / |
    # ( 1 )----(-1 )|
    #   |  |     |  |
    #   |( 5 ) - |( 3 )-> i
    #   | /      | /
    #   |/       |/
    # ( 3 )----( 1 )
    #  /
    # k
    mosaic = get_mosaic( Cell(SVector(0., 0., 0.), SVector(1., 1.,1.), reshape([5 3 3 1 3 1 1 -1], 2, 2, 2)) )
    @test length(mosaic)==15 && sum([vtktypes[(3,4)]==tile.vtktype for tile in mosaic])==12 && sum([vtktypes[(3,5)]==tile.vtktype for tile in mosaic])==3 && size(hcat([tile.coords for tile in mosaic]...))==(3,63)

    vtk = write_vtk("mosaic", mosaic)

end # testset

@testset verbose = true "Cardinal B-spline" begin
    @test cardinal_stencil(0)≈[1.]
    @test cardinal_stencil(1)≈[1/2,1/2]
    @test cardinal_stencil(2)≈[1/6,2/3,1/6]
    @test cardinal_stencil(3)≈[1/24,11/24,11/24,1/24]
    @test cardinal_indices(2, 4)==(2:6)
end # testset