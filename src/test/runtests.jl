include("test_file.jl")

using Test




@testset "Test Suite" begin
    @test 1 + 1 == 2
    @test 2 * 3 == 6
    @test 10 - 5 == 5
end




