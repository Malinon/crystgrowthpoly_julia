using Nemo
include("tessellation.jl")
include("visualisation.jl")
include("preprocessor.jl")
# TODO: Choose better type for A


function get_topological_growth_polynomials(cells::AbstractVector, translation_vectors::AbstractVector, symmetric_growth::Bool=true)
    """Function finding topological growth functions"""
    dim = length(cells) - 1
    if symmetric_growth
        args = ntuple(i -> (i,i), dim+1)
        variables_num = 1
    else
        # TODO: Implement generating arguments for higher dimmensions.
        if dim == 2
            args = ((1,1), (1,2), (2,1), (2,2))
        elseif dim == 3
            args = ((1,1,1), (1,1,2), (1,2,1), (2,1,1), (2,2,1), (2,1,2), (1,2,2), (2,2,2))
        end
        variables_num = dim
    end
    polynomial_k_cells_coeff_lists = Vector{fmpq_mat}(undef, dim + 1)
    # Calculate coefficients of growth function for cells with dimmension from 0 to dim - 2
    gen_val = arg -> cryst_private.get_k_cells_num(arg, cells[1], translation_vectors, 0)
    polynomial_k_cells_coeff_lists[1] = cryst_private.find_poly(gen_val, args, symmetric_growth)
    for i in 2:(length(cells) - 2)
        gen_val = arg -> cryst_private.get_k_cells_num(arg, cells[i], translation_vectors, 1)
        polynomial_k_cells_coeff_lists[i] = cryst_private.find_poly(gen_val, args, symmetric_growth)
    end 
    # Calculate coefficients of growth function related to cells of highest dimmensions.
    highest_dim_function = zeros(Integer, length(args))
    highest_dim_function[1] = length(cells[length(cells)])
    polynomial_k_cells_coeff_lists[dim + 1] = matrix(QQ, 1,  length(args), highest_dim_function)#highest_dim_function[0]
    # Calculate coefficients of growth functions related to dim - 1 cells based on Euler Characteristic
    polynomial_k_cells_coeff_lists[dim] = (
        cryst_private.calculate_coefficients_based_ne_euler_charcteristics(
            polynomial_k_cells_coeff_lists, highest_dim_function, dim))
    return Tuple(gf.growth_function((polynomial_k_cells_coeff_lists[i],), variables_num, i) for i in 1:dim+1)
end

function get_crystallographic_growth_functions(points::AbstractVector, edges::AbstractVector, faces::AbstractVector, x0, symmetric_growth::Bool=true)
        return (gf.growth_function(cryst_private.get_vertices_polynomial(points, x0, symmetric_growth), 2, 1),
        gf.growth_function(cryst_private.get_faces_polynomial(edges, x0, symmetric_growth), 2, 2),
    gf.growth_function(cryst_private.get_faces_polynomial(faces, x0, symmetric_growth), 2, 3))
end

function get_growth_polynomials_parallelogram(tes::Tessellation, scaler_v1, scaler_v2, x0, symmetric_growth::Bool=false)
    return get_crystallographic_growth_functions(tes.polygon.cells[1], tes.polygon.cells[2], tes.polygon.cells[3], x0, symmetric_growth)
end

function read_tessellation_from_file(file_path)
    f = open(file_path, "r")
    nums_of_cells = split(readline(f), " ", keepempty=false)
    num_vertices = parse(Int64, (nums_of_cells[1]))
    vertices = Vector{Tuple{fmpq, fmpq}}(undef, num_vertices)
    for i in 1:num_vertices
        vertices[i] = line_to_point(readline(f))
    end
    num_edges = parse(Int64, nums_of_cells[2])
    edges = Vector{Tuple{Tuple{fmpq, fmpq}, Tuple{fmpq, fmpq}}}(undef, num_edges)
    for i in 1:num_edges
        vertices_ids = split(readline(f), " ", keepempty=false)
        edges[i] = ids_to_faces(vertices, vertices_ids)
    end
    num_faces = parse(Int64,nums_of_cells[3])
    faces = Vector{Tuple}(undef, num_faces)
    for i in 1:num_faces
        vertices_ids = split(readline(f), " ", keepempty=false)
        faces[i] = ids_to_faces(vertices, vertices_ids)
    end
    close(f)
    return Tessellation(Polygon(vertices, edges, faces))
end


