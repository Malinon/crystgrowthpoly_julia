include("growth_function.jl")
include("regions.jl")
include("orphic_utils.jl")

function __init__()
    global C = CalciumField()
end

struct OrphicDiagram
    """Structure describing orphic diagram"""
    rectangles::Vector{OpenRectangle}
    lines::Vector{OpenLine}
    points::Vector{SpecialPoint}
end

function get_growth_for_edge(reference_polynomials::Vector{Int64}, number_of_vertices::Int64, points_at_line::Int64,
    modifiers_edge::Vector{ModifierChangeDescription}, modifiers_face::ModifierChangeDescription, direction::Direction)
    """Compute growth functions for horizontal edge"""
    vertice_polynomial = [number_of_vertices, 0, 0, 0]
    vertice_polynomial[3 - direction] = points_at_line

    edge_polynomial = calculate_function_for_edge(reference_polynomials[1], modifiers_edge, direction)
    face_polynomial = calculate_function_for_edge(reference_polynomials[2], modifiers_face, direction)

    return (vertice_polynomial, edge_polynomial, face_polynomial)
end

function get_growth_for_rectangle(reference_polynomials, number_of_vertices, modifiers_edge, modifiers_face, direction)
    """Compute growth functions for horizontal edge"""
    vertice_polynomial = [number_of_vertices, 0, 0, 0]

    edge_polynomial = calculate_function_for_rectangle(reference_polynomials[1], modifiers_edge, direction)
    face_polynomial = calculate_function_for_rectangle(reference_polynomials[2], modifiers_face, direction)

    return (vertice_polynomial, edge_polynomial, face_polynomial)
end

function create_orphic_diagram(tessellation::Tessellation)
    """Find regions within which growth functions are identical."""
    x_endpoints, x_vert_counts = find_boundary_lines(tessellation.vertices, X)
    y_endpoints, y_vert_counts = find_boundary_lines(tessellation.vertices, Y)

    # Const values
    number_of_vertices = length(tessellation.vertices)
    x_len = length(x_endpoints)
    y_len = length(y_endpoints)
    lengths = (x_len, y_len)

    # Assign integers to coordinates to eliminate operations on rational/real numbers
    x_dict = Dict((x_endpoints[i], i) for i in 1:x_len)
    y_dict = Dict((y_endpoints[i], i) for i in 1:y_len)
    dicts = (x_dict, y_dict)

    # Generating descriptions and modifiers for edges
    modifiers_edge, descriptions_edges = describe_change_points_and_cells(tessellation.edges, lengths, dicts)
    modifiers_face, descriptions_faces = describe_change_points_and_cells(tessellation.faces, lengths, dicts)

    # Buffers for regions
    found_rectangles = Vector{OpenRectangle}(undef, x_len * y_len)
    points = Vector{SpecialPoint}(undef, x_len * y_len)
    lines = Vector{OpenLine}(undef, 2 * x_len * y_len)

    point_counter = 1
    for y_iter in 1:y_len
        for x_iter in 1:x_len
            # Find growth functions for 0-dimensional region
            if (x_iter != 1) || (y_iter != 1)
                # It is not first vertice, so we will use growth function from neighbour vertice
                # Determine from which direction we should take growth functions
                if x_iter == 1
                    # We are in first column, so we should take growth functions from below
                    id_of_reference_point = (y_iter - 1) * x_len + x_iter
                    shared_coord = y_iter
                    diff_coord = x_iter
                    direction_of_the_neighbour = Y
                    perpendicular_direction = X
                else
                    # We are not in first column, so we should take growth functions from left
                    id_of_reference_point = point_counter - 1
                    shared_coord = x_iter
                    diff_coord = y_iter
                    direction_of_the_neighbour = X
                    perpendicular_direction = Y
                end

                new_poly_edge = calculate_polynomial_in_next_point(points[id_of_reference_point].growth_f[2],
                                                                   modifiers_edge[direction_of_the_neighbour][shared_coord - 1],
                                                                   modifiers_edge[direction_of_the_neighbour][shared_coord],
                                                                   descriptions_edges[perpendicular_direction],
                                                                   direction_of_the_neighbour,
                                                                   diff_coord)
                new_poly_faces = calculate_polynomial_in_next_point(points[id_of_reference_point].growth_f[3],
                                                                    modifiers_face[direction_of_the_neighbour][shared_coord - 1],
                                                                    modifiers_face[direction_of_the_neighbour][shared_coord],
                                                                    descriptions_faces[perpendicular_direction],
                                                                    direction_of_the_neighbour,
                                                                    diff_coord)
            else
                """Calculate polynomials for first 0-dimensional region """
                new_poly_edge = get_coefs_poly_point(descriptions_edges, x_iter, y_iter)
                new_poly_faces = get_coefs_poly_point(descriptions_faces, x_iter, y_iter)
            end

            new_poly_points = [number_of_vertices, y_vert_counts[y_iter], x_vert_counts[x_iter], 0]
            points[point_counter] = SpecialPoint((x_endpoints[x_iter], y_endpoints[y_iter]),
                                                 (new_poly_points, new_poly_edge, new_poly_faces))

            # Find regions for 1/2-dimensional adjacent regions
            # Find polynomials for horizontal edge
            next_x = (x_iter == length(x_endpoints)) ? x_endpoints[1] + 1 : x_endpoints[x_iter + 1]
            horizontal_edge_functions = get_growth_for_edge(points[point_counter].growth_f,
                                                            points_at_line,
                                                            modifiers_edge[X],
                                                            modifiers_face[X],
                                                            X)
            lines[point_counter * 2 - 1] = OpenLine((x_endpoints[x_iter], y_endpoints[y_iter]),
                                                    (next_x, y_endpoints[y_iter]),
                                                    horizontal_edge_functions)

            # Find polynomials for vertical edge
            next_y = (y_iter == length(y_endpoints)) ? y_endpoints[1] + 1 : y_endpoints[y_iter + 1]
            vertical_edge_functions = get_growth_for_edge(points[point_counter].growth_f,
                                                          points_at_line,
                                                          modifiers_edge[Y],
                                                          modifiers_face[Y],
                                                          Y)
            lines[point_counter * 2] = OpenLine((x_endpoints[x_iter], y_endpoints[y_iter]),
                                                (x_endpoints[x_iter], next_y),
                                                vertical_edge_functions)

            # Find regions for 2-dimensional adjacent region
            rectangle_functions = get_growth_for_rectangle(horizontal_edge_functions,
                                                           modifiers_edge[Y],
                                                           modifiers_face[Y],
                                                           X)
            found_rectangles[point_counter] = OpenRectangle((next_x, next_y),
                                                            (x_endpoints[x_iter], y_endpoints[y_iter]),
                                                            rectangle_functions)
            point_counter += 1
        end
    end

    # Add +1 to growth function for 0d regions with vertices
    for v in tessellation.vertices
        points[x_dict[v[1]] + (y_dict[v[2]] - 1) * x_len].growth_f[1][4] = 1
    end

    return OrphicDiagram(found_rectangles, lines, points)
end