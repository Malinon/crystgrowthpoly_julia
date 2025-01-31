import Pkg

Pkg.develop(path="../../crystgrowthpoly_julia")

using Test
import crystgrowthpoly

function get_all_files_from_directory(directory)
    files = readdir(directory)
    return [joinpath(directory, file) for file in files]
end

INPUT_FILES = get_all_files_from_directory("tessellations")

function is_valid_orphic_diagram(orph)
    @test typeof(orph) == crystgrowthpoly.OrphicDiagram
    points_num = length(orph.points)
    @test points_num > 0

    @test length(orph.lines) == 2 * points_num
    @test length(orph.rectangles) == points_num
    rectangle_poly = orph.rectangles[1].growth_f[1]
    @test rectangle_poly[2] == 0
    @test rectangle_poly[3] == 0
    @test rectangle_poly[4] == 0
    for rectangle in orph.rectangles
        @test rectangle.growth_f[1] == rectangle_poly
    end
end

@testset "Generic Test Orphic Diagram" begin
    for path in INPUT_FILES
        @testset "Testing file: $path" begin
            tes = crystgrowthpoly.read_tessellation_from_file(path)
            reduced_tes = crystgrowthpoly.reduce_polygon(tes.vertices, tes.edges, tes.faces)
            orph = crystgrowthpoly.create_orphic_diagram(reduced_tes)
            is_valid_orphic_diagram(orph)
        end
    end
end

@testset "4444 test" begin
    tes = crystgrowthpoly.read_tessellation_from_file("tessellations/4444.txt")
    reduced_tes = crystgrowthpoly.reduce_polygon(tes.vertices, tes.edges, tes.faces)
    orph = crystgrowthpoly.create_orphic_diagram(reduced_tes)
    @test length(orph.rectangles) == 1
    @test orph.points[1].growth_f == ([1,1,1,1], [2, 1, 1, 0], [1, 0, 0, 0])
    @test orph.points[1].point == (0,0)
    @test orph.rectangles[1].growth_f == ([1,0,0,0], [2, -1, -1, 0], [1, -1, -1, 1])
    @test orph.rectangles[1].higher_right_corner == (1,1)
    @test orph.rectangles[1].lower_left_corner == (0,0)
    @test orph.lines[1].growth_f == ([1,1,0,0], [2, 1, -1, -1], [1, 0, -1, 0])
    @test orph.lines[1].start == (0,0)
    @test orph.lines[1].end_point == (1,0)
    @test orph.lines[2].growth_f == ([1,0,1,0], [2, -1, 1, -1], [1, -1, 0, 0])
    @test orph.lines[2].start == (0,0)
    @test orph.lines[2].end_point == (0,1)
end

@testset "884 test" begin
    tes = crystgrowthpoly.read_tessellation_from_file("tessellations/884.txt")
    reduced_tes = crystgrowthpoly.reduce_polygon(tes.vertices, tes.edges, tes.faces)
    orph = crystgrowthpoly.create_orphic_diagram(reduced_tes)
    @test length(orph.rectangles) == 4
    # Check the growth functions of the points
    @test orph.points[1].growth_f == ([2,1,1,1], [6, 1, 1, 0], [4, 0, 0, 0])
    @test orph.points[1].point == (0,0)
    @test orph.points[2].growth_f == ([2,1,1,0], [6, 1, -1, -1], [4, 0, -2, 0])
    @test orph.points[2].point == (0.5, 0)
    @test orph.points[3].growth_f == ([2,1,1,0], [6, -1, 1, -1], [4, -2, 0, 0])
    @test orph.points[3].point == (0, 0.5)
    @test orph.points[4].growth_f == ([2,1,1,1], [6, -1, -1, 0], [4, -2, -2, 0])
    @test orph.points[4].point == (0.5, 0.5)
    # Check the growth functions of the rectangles
    @test orph.rectangles[1].growth_f == ([2,0,0,0], [6, -3, -3, 1], [4, -3, -3, 2])
    @test orph.rectangles[1].higher_right_corner == (0.5,0.5)
    @test orph.rectangles[1].lower_left_corner == (0,0)
    @test orph.rectangles[2].growth_f == orph.rectangles[1].growth_f # Same growth function
    @test orph.rectangles[2].higher_right_corner == (1,0.5)
    @test orph.rectangles[2].lower_left_corner == (0.5,0)
    @test orph.rectangles[3].growth_f == orph.rectangles[1].growth_f # Same growth function
    @test orph.rectangles[3].higher_right_corner == (0.5,1)
    @test orph.rectangles[3].lower_left_corner == (0,0.5)
    @test orph.rectangles[4].growth_f == orph.rectangles[1].growth_f # Same growth function
    @test orph.rectangles[4].higher_right_corner == (1,1)
    @test orph.rectangles[4].lower_left_corner == (0.5,0.5)
    # Check the growth functions of the lines
    @test orph.lines[1].growth_f == ([2,1,0,0], [6, 1, -3, -1], [4, 0, -3, 0])
    @test orph.lines[1].start == (0,0)
    @test orph.lines[1].end_point == (0.5,0)
    @test orph.lines[2].growth_f == ([2,0,1,0], [6, -3, 1, -1], [4, -3, 0, 0])
    @test orph.lines[2].start == (0,0)
    @test orph.lines[2].end_point == (0,0.5)
    @test orph.lines[3].growth_f == orph.lines[1].growth_f # Same growth function
    @test orph.lines[3].start == (0.5,0)
    @test orph.lines[3].end_point == (1,0)
    @test orph.lines[4].growth_f == ([2,0,1,0], [6, -3, -1, 0], [4, -3, -2, 1])
    @test orph.lines[4].start == (0.5,0)
    @test orph.lines[4].end_point == (0.5,0.5)
    @test orph.lines[5].growth_f == ([2,1,0,0], [6, -1, -3, 0], [4, -2, -3, 1])
    @test orph.lines[5].start == (0,0.5)
    @test orph.lines[5].end_point == (0.5,0.5)
    @test orph.lines[6].growth_f == orph.lines[2].growth_f
    @test orph.lines[6].start == (0,0.5)
    @test orph.lines[6].end_point == (0,1)
    @test orph.lines[7].growth_f == orph.lines[5].growth_f
    @test orph.lines[7].start == (0.5,0.5)
    @test orph.lines[7].end_point == (1,0.5)
    @test orph.lines[8].growth_f == orph.lines[4].growth_f 
    @test orph.lines[8].start == (0.5,0.5)
    @test orph.lines[8].end_point == (0.5,1)
end
