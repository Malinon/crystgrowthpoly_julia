include("growth_function.jl")
include("regions.jl")

function __init__()
  global C = CalciumField()
end


""" This function finds breakpoints defining regions in orphic diagrams. It also counts number of unique  vertices on each line.
It returns coordinates describing one-dimmension of orphic diagrams and number of translation-unique vertices with these coordinates"""
function find_boundary_lines(tes::Vector{Tuple{A, A}}, index) where A
    """Finds lines perpendicular to chosen axis defining domains"""
    unique_points_coord = Dict{A, Int64}()
    # Find unique coordinates and count vertices on these lines
    """ Every unique x/y coordinate in vertex's description defines line defining orphic diagrams"""
    for p in tes
        unique_points_coord[p[index]] = get(unique_points_coord, p[index], 0) + 1
    end
    """ x- and y- endpoints combained describe shapes of orphic diagrams """
    endpoints = collect(keys(unique_points_coord))
    sort!(endpoints)
    modifiers = [unique_points_coord[c] for c in endpoints]::Array{Int64}
    return (endpoints, modifiers)
end

function describe(self::OpenRectangle, v1, v2)
    p2 = ((self.higher_right_corner[1], self.lower_left_corner[2]))
    p4 = ((self.lower_left_corner[1], self.higher_right_corner[2]))
    println("Parallelogram ", (self.lower_left_corner, p2, self.higher_right_corner, p4))
end

function get_polynomials(self)
        return(self.growth_f[1].polynomials_coefficients, self.growth_f[2].polynomials_coefficients, self.growth_f[3].polynomials_coefficients)
end

function describe(self::OpenLine, v1, v2)
    println("Line ", self.start, self.end_point)
end

function describe(self::SpecialPoint, v1, v2)
    println("Point ", self.point)
end

struct RegionsDescription
    """This class describes all regions of orphic diagram."""
    tessellation
    description_dict
end

"""This struct contains all usefull informations from edge/face description"""
struct MinimalFigDesc
    """ Minimal coordinate of figure"""
    zeroer_id::Int64
    """ Maximal coordinate of figure"""
    edger_id::Int64
    """ Maximal modifier of the figure. X-modifier of figure is equal  for given x0 is minimal k for which, there exist equivalent Fig,
    that for all p in Fig,   k > p.x - x_0 >= 0.  """
    max_mod_val::Int64
end

struct EndpointModifier
    """ Id of corresponding MinimalFigDesc """
    id::Int64
    """ Is coordinate equal to maximal coordinate of figure"""
    is_edger::Bool
end

function generate_modifier_and_description(x_dict, fig_raw, x_modifiers, index, id)
    min_max = cryst_private.get_min_max(fig_raw, index)# Find extreme coordinates of figure
    floored_min = floor(min_max[1])
    floored_max = floor(min_max[2])
    mantisee_max = min_max[2] - floored_max
    mantisee_min = min_max[1] - floored_min
    # Replace coordinates with id numbers
    zeroer_id = x_dict[mantisee_min]
    edger_id = x_dict[mantisee_max]
    modifier = Int(ZZ(floor(min_max[2] - min_max[1]))) # It is minimal modifier
    # Save information. Impact/status of the figur changes in mantisee_min and mantisee_max lines
    push!(x_modifiers[edger_id], EndpointModifier(id, true))
    push!(x_modifiers[zeroer_id], EndpointModifier(id, false))
    if edger_id != zeroer_id
        # Maximal modifier is equal to minimal modifier + 1
        return MinimalFigDesc(zeroer_id, edger_id, modifier + 1)
    else
        # If x/y-lenth of figure is integer, then modifier is constant.
        return MinimalFigDesc(zeroer_id, edger_id, modifier)

    end
end

function generate_modifiers_array_and_fig_description(x_len, y_len, x_dict, y_dict, figs_raw, x_modifiers, y_modifiers, fig_desriptions)
    """ Generate modifiers and minimal figure descriptions for edges or faces. """
    for i in 1:x_len
        x_modifiers[i] = Vector{EndpointModifier}()
    end
    for i in 1:y_len
        y_modifiers[i] = Vector{EndpointModifier}()
    end
    id = 1
    for fig_raw in figs_raw
        x_desc = generate_modifier_and_description(x_dict, fig_raw, x_modifiers, 1, id)
        y_desc = generate_modifier_and_description(y_dict, fig_raw, y_modifiers, 2, id)
        fig_desriptions[id] = (x_desc, y_desc)
        id += 1
    end
end



""" Update polynomial if impact of figure has changed. Version for 2-dimmensional regions"""
function simple_update!(coeffs, modifier, fig_desc,index, mod_iter)
    if modifier.is_edger
        """ If it is edger, then for points on the right we get new_mod_x = mod_x -1
        Impact of fig before = (x -mod_x)(y -mod_y) +edge_modifiers_x + edge_modifiers_y + edge_modifiers_both
        Impact of fig after = (x -mod_x + 1)(y -mod_y)+ edge_modifiers_y
        Impact of fig_after = Impact of fig before - mod_y + y
        """
        """ Remove edge influence
        Single edge modificator before = y - mod_y
        Both edge modfdificatot before = y - mod_y + x - mod_x +1
        Single edge modificator before = 0
        Both edge modfdificatot before = x - mod_x + 1  (Because now it is egdger on y)

        So in this case nothing happens
        """
    else
        """ If it is zeroer, then for points on the right we get new_mod_x = mod_x + 1
        Impact of fig before = (x -mod_x)(y -mod_y) + edge_modifiers_y
        Impact of fig after = (x -mod_x - 1)(y -mod_y)+ edge_modifiers_y
        Impact of fig_after = Impact of fig before + mod_y - y
        """
        coeffs[4 - index] -= 1# -y
        """ Remove edge influence
        Single edge modificator before = x - mod_x
        Single edge modificator after = x - mod_x - 1
        """
        if is_edgy(mod_iter, fig_desc)
            coeffs[4] += (get_modifier_0d(mod_iter, fig_desc) - 1) # mod_y -1
        else
            coeffs[4] += get_modifier_0d(mod_iter, fig_desc)# Mod_y
        end
    end
end



""" Update polynomial if impact of figure has changed. Version for 2-dimmensional regions"""
function simple_update_rect!(coeffs, modifier, fig_desc,index, mod_iter)
    if modifier.is_edger
        """ If it is edger, then for points on the right we get new_mod_x = mod_x -1
        Impact of fig before = (x -mod_x)(y -mod_y) +edge_modifiers_x + edge_modifiers_y + edge_modifiers_both
        Impact of fig after = (x -mod_x + 1)(y -mod_y)+ edge_modifiers_y
        Impact of fig_after = Impact of fig before - mod_y + y
        """
        """ Remove edge influence
        Single edge modificator before = y - mod_y
        Both edge modfdificatot before = y - mod_y + x - mod_x +1
        Single edge modificator before = 0
        Both edge modfdificatot before = x - mod_x + 1  (Because now it is now egdger on y)

        So in this case nothing happens
        """
    else
        """ If it is zeroer, then for points on the right we get new_mod_x = mod_x + 1
        Impact of fig before = (x -mod_x)(y -mod_y) + edge_modifiers_y
        Impact of fig after = (x -mod_x - 1)(y -mod_y)+ edge_modifiers_y
        Impact of fig_after = Impact of fig before + mod_y - y
        """

        """
        It is 2-dimmenstional region so no figures are on edges.
        """
        coeffs[2] -= 1
        """We use another function for calculating modifier, because we already have moved from point (dict_x[x_iter], dict_y[y_iter]) to 
        (dict_x[x_iter] + epsilion , dict_y[y_iter]). """ 
        coeffs[4] += get_modifier(mod_iter, fig_desc)
    end
end

function is_edgy(index::Int64, fig::MinimalFigDesc)
    return index == fig.edger_id
end

# TODO: Check if we can eliminate <,= >=, here 
""" Get modifier for figure """
function get_modifier_0d(index, fig::MinimalFigDesc)
    if fig.zeroer_id > fig.edger_id
        if index > fig.zeroer_id || index <= fig.edger_id
            return fig.max_mod_val
        else
            return fig.max_mod_val - 1
        end
    elseif fig.zeroer_id == fig.edger_id
        return fig.max_mod_val
    else
        if index > fig.edger_id || index <= fig.zeroer_id
            return fig.max_mod_val - 1
        else
            return fig.max_mod_val
        end
    end
end

""" Get modifier for figure. Strength of inequalities are changes because we already have moved by espilon in x-direction """
function get_modifier(index, fig::MinimalFigDesc)
    if fig.zeroer_id > fig.edger_id
        if index >= fig.zeroer_id || index < fig.edger_id
            return fig.max_mod_val
        else
            return fig.max_mod_val - 1
        end
    elseif fig.zeroer_id == fig.edger_id
        return fig.max_mod_val
    else
        if index >= fig.edger_id || index < fig.zeroer_id
            return fig.max_mod_val - 1
        else
            return fig.max_mod_val
        end
    end
end


"""Calculate growth function for 0-dimmensional region based on functions for another 0-dimmensional regions"""
function update_polynomial_in_next_point!(coefs_last, modifiers, minimal_fig_desc_list, mod_iter, index, row_iter)
    new_poly = copy(coefs_last)
    # Last point modification
    """ Apply modification defined in previous node. In previous node some figures changed their status.
    If the last node wasn't edger for them their change of their impact starts here"""
    for mod in modifiers[mod_iter - 1]
        simple_update!(new_poly, mod, minimal_fig_desc_list[mod.id][3 - index], index, row_iter)
    end
    # Edgy modifiers
    """ Apply edgy modifiers. Some figures are now inside closuer of half-opend frame, so their impact has changers.
    Main modiers are not changed but we need to change edge impact/modifiers """
    for mod in modifiers[mod_iter]
        if mod.is_edger
            """ fig is now on edge.
            Old edge_modifier sinlge: 0
            New edge_modifier sinlge: y - mod_y
            """
            fig = minimal_fig_desc_list[mod.id][3-index]
            new_poly[4-index] +=  1
            new_poly[4] -= get_modifier_0d(row_iter, fig)
            if is_edgy(row_iter, fig)
                new_poly[4] += 1
            end
        end
    end
    return new_poly
end


function common_update!(coefficients, mod_k, mod_l)
    coefficients[2]  = coefficients[2] - mod_l
    coefficients[3]  = coefficients[3] - mod_k
    coefficients[4] =  coefficients[4] + mod_k * mod_l
end


""" Add impact of figure to growth polynomial """
function update_polynomial_point(edge_desc::Tuple{MinimalFigDesc, MinimalFigDesc}, index_k, index_l, coefficients)
    mod_k = get_modifier_0d(index_k, edge_desc[1])
    mod_l = get_modifier_0d(index_l, edge_desc[2])
    common_update!(coefficients, mod_k, mod_l)
    if is_edgy(index_k, edge_desc[1])
         coefficients[3] = coefficients[3] + 1
         coefficients[4] = coefficients[4] - mod_l
        if is_edgy(index_l, edge_desc[2])
             coefficients[2] = coefficients[2] + 1
             coefficients[4] = coefficients[4] + 1 - mod_k
       end
    elseif is_edgy(index_l, edge_desc[2])
         coefficients[2] = coefficients[2] + 1
         coefficients[4] = coefficients[4] - mod_k
    end
end


""" Find growth function (1-cells, 2-cells) for 0 dimmensional region """
function get_coefs_poly_point(face_descs, index_x, index_y)
    coefficients = [length(face_descs), 0, 0, 0]::Array{Int64}
    for desc in face_descs
        update_polynomial_point(desc, index_x, index_y,  coefficients)
    end
    return coefficients
end

""" Based on growth functions in 0-dimmensional region find polynomials for adjacent regions (verical and horsional line starting in this point
and rectangle for which the point is lower-left corner)"""
function update_lines_and_rect_poly(poly_horisontal_line, poly_vertical_line, x_iter, y_iter, x_modifier, y_modifier, fig_desc)
    # Horisontal line
    for mod in x_modifier
        simple_update!(poly_horisontal_line, mod, fig_desc[mod.id][2], 1, y_iter)
    end
    # Vertical line and rectangle
    poly_rect = copy(poly_horisontal_line)
    for mod in y_modifier
        fig = fig_desc[mod.id][1]
        simple_update!(poly_vertical_line, mod, fig, 2, x_iter)
        simple_update_rect!(poly_rect, mod, fig, 2, x_iter)
    end
    return poly_rect
end

""" Constract domains and insert into corresponfing arrays"""
function add_regions(x_iter, y_iter, lines, found_rectangles, counter, x_endpoints, y_endpoints, rect_funs, horisontal_line_funs, verical_line_funs)
    next_x = (x_iter == length(x_endpoints[1])) ? x_endpoints[1][1] + 1 : x_endpoints[1][x_iter + 1]
    next_y = (y_iter == length(y_endpoints[1])) ? y_endpoints[1][1] + 1 : y_endpoints[1][y_iter + 1]
    found_rectangles[counter] = OpenRectangle((next_x, next_y), (x_endpoints[1][x_iter], y_endpoints[1][y_iter]), rect_funs)
    lines[counter * 2 - 1] = OpenLine((x_endpoints[1][x_iter], y_endpoints[1][y_iter]), (next_x, y_endpoints[1][y_iter]), horisontal_line_funs)
    lines[counter * 2] = OpenLine((x_endpoints[1][x_iter], y_endpoints[1][y_iter] ), (x_endpoints[1][x_iter], next_y), verical_line_funs)
end

""" Constract domains and insert into corresponfing arrays"""
function add_regions_0d(x_iter, y_iter, lines, found_rectangles, counter, x_endpoints, y_endpoints, rect_funs, horisontal_line_funs, verical_line_funs)
    next_x = (x_iter == length(x_endpoints[1])) ? x_endpoints[1][1] + 1 : x_endpoints[1][x_iter + 1]
    next_y = (y_iter == length(y_endpoints[1])) ? y_endpoints[1][1] + 1 : y_endpoints[1][y_iter + 1]
    found_rectangles[counter] = OpenRectanglePrimitive((next_x, next_y), (x_endpoints[1][x_iter], y_endpoints[1][y_iter]), rect_funs)
    lines[counter * 2 - 1] = OpenLinePrimitive((x_endpoints[1][x_iter], y_endpoints[1][y_iter]), (next_x, y_endpoints[1][y_iter]), horisontal_line_funs)
    lines[counter * 2] = OpenLinePrimitive((x_endpoints[1][x_iter], y_endpoints[1][y_iter] ), (x_endpoints[1][x_iter], next_y), verical_line_funs)
end

function find_regions(vertices::Vector{Tuple{A,A}}, ts_edges::Vector{Tuple{Tuple{A,A}, Tuple{A,A}}},
        ts_faces ,symmetric_growth::Bool=true, full_plot::Bool=false) where A
    """Find regions within which growth functions are identical."""
    x_endpoints = find_boundary_lines(vertices, 1)
    y_endpoints = find_boundary_lines(vertices, 2)
    println(x_endpoints)
    println(y_endpoints)

    # Const values
    number_of_vertices = length(vertices)
    x_len = length(x_endpoints[1])
    y_len = length(y_endpoints[1])
    # Assign integers to coordinates to eliminate operations on rational/real numbers
    x_dict = Dict((x_endpoints[1][i], i) for i in 1:length(x_endpoints[1]))
    y_dict = Dict((y_endpoints[1][i], i) for i in 1:length(y_endpoints[1]))
    # Generating desciptions and modifiers for edges
    edge_desriptions = Vector{Tuple{MinimalFigDesc, MinimalFigDesc}}(undef, length(ts_edges))
    x_modifiers_edges = Vector{Vector{EndpointModifier}}(undef, length(x_endpoints[1]))
    y_modifiers_edges = Vector{Vector{EndpointModifier}}(undef, length(y_endpoints[1]))
    generate_modifiers_array_and_fig_description(x_len, y_len, x_dict, y_dict, ts_edges, x_modifiers_edges, y_modifiers_edges, edge_desriptions)
    # Generating desciptions and modifiers for faces
    face_desriptions = Vector{Tuple{MinimalFigDesc, MinimalFigDesc}}(undef, length(ts_faces))
    x_modifiers_faces = Vector{Vector{EndpointModifier}}(undef, length(x_endpoints[1]))
    y_modifiers_faces = Vector{Vector{EndpointModifier}}(undef, length(y_endpoints[1]))
    generate_modifiers_array_and_fig_description(x_len, y_len, x_dict, y_dict, ts_faces, x_modifiers_faces, y_modifiers_faces, face_desriptions)
    # Buffers for regions
    found_rectangles = Vector{OpenRectangle}(undef, x_len * y_len)
    points = Vector{SpecialPoint}(undef, x_len * y_len)
    lines = Vector{OpenLine}(undef, 2 * x_len * y_len)

    counter = 1
    for y_iter in 1:y_len
        for x_iter in 1:x_len
            if x_iter != 1
                # Set next points
                """Calculate polynomials for 0-dimmensional region based on polynomials for corresponding point on left"""
                new_poly_points = [number_of_vertices, y_endpoints[2][y_iter], x_endpoints[2][x_iter], 0]
                new_poly_edge = update_polynomial_in_next_point!(points[counter - 1].growth_f[2], x_modifiers_edges, edge_desriptions, x_iter ,1, y_iter)
                new_poly_faces = update_polynomial_in_next_point!(points[counter - 1].growth_f[3], x_modifiers_faces, face_desriptions, x_iter ,1, y_iter)
            else
                if y_iter != 1
                    """Calculate polynomials for 0-dimmensional region based on polynomials for corresponding point 1 row below"""
                    new_poly_points = [number_of_vertices, y_endpoints[2][y_iter], x_endpoints[2][x_iter], 0]
                    new_poly_edge = update_polynomial_in_next_point!(points[(y_iter - 2) * x_len + 1].growth_f[2], y_modifiers_edges, edge_desriptions, y_iter , 2, x_iter)
                    new_poly_faces = update_polynomial_in_next_point!(points[(y_iter - 2) * x_len + 1].growth_f[3], y_modifiers_faces, face_desriptions, y_iter , 2, x_iter)
                else
                    """Calculate polynomials for first 0-dimmensional region """
                    new_poly_points = [number_of_vertices, y_endpoints[2][y_iter], x_endpoints[2][x_iter], 0]
                    new_poly_edge =  get_coefs_poly_point(edge_desriptions, x_iter, y_iter)
                    new_poly_faces = get_coefs_poly_point(face_desriptions, x_iter, y_iter)
                end
            end
            points[counter] = SpecialPoint((x_endpoints[1][x_iter], y_endpoints[1][y_iter]), (new_poly_points, new_poly_edge, new_poly_faces))
            """Start finding growth functions for adjacent regions """
            # Vertices polynomials
            poly_points_horisontal_line = [number_of_vertices, y_endpoints[2][y_iter], 0, 0]
            poly_points_vertical_line = [number_of_vertices, 0, x_endpoints[2][x_iter], 0]
            poly_points_rect = [number_of_vertices, 0, 0, 0]
            """ Based on growth functions in 0-dimmensional region find polynomials for adjacent regions (verical and horsiontal line starting in this point
            and rectangle for which the point is lower-left corner)"""
            # Edges polynomials
            poly_edges_horisontal_line = copy(points[counter].growth_f[2])
            poly_edges_vertical_line = copy(points[counter].growth_f[2])
            poly_edges_rect = update_lines_and_rect_poly( poly_edges_horisontal_line, poly_edges_vertical_line, x_iter, y_iter,
                x_modifiers_edges[x_iter], y_modifiers_edges[y_iter], edge_desriptions)
            # Faces polynomials
            poly_faces_horisontal_line = copy(points[counter].growth_f[3])
            poly_faces_vertical_line = copy(points[counter].growth_f[3])
            poly_faces_rect = update_lines_and_rect_poly( poly_faces_horisontal_line, poly_faces_vertical_line, x_iter, y_iter,
                x_modifiers_faces[x_iter], y_modifiers_faces[y_iter], face_desriptions)
            # Construct regions and assign found growth functions
            add_regions(x_iter, y_iter, lines, found_rectangles, counter, x_endpoints, y_endpoints,
                (poly_points_rect, poly_edges_rect, poly_faces_rect),
                (poly_points_horisontal_line, poly_edges_horisontal_line, poly_faces_horisontal_line),
                (poly_points_vertical_line, poly_edges_vertical_line, poly_faces_vertical_line))
            counter += 1
        end
    end

    # Add corner modifier to vertice polynomials (0d regions only). (If frames starts at vertice we need to add +1 to point polynomial)
    for v in vertices
        points[x_dict[v[1]] + (y_dict[v[2]] - 1) * x_len].growth_f[1][4] = 1
    end

    return [found_rectangles, lines, points]
end

function find_regions(vertices::Vector{Tuple{A,A}},symmetric_growth::Bool=true, full_plot::Bool=false) where A
    """Find regions within which growth functions are identical."""
    x_endpoints = find_boundary_lines(vertices, 1)
    y_endpoints = find_boundary_lines(vertices, 2)
    # Const values
    number_of_vertices = length(vertices)
    x_len = length(x_endpoints[1])
    y_len = length(y_endpoints[1])
    # Assign integers to coordinates to eliminate operations on rational/real numbers
    x_dict = Dict((x_endpoints[1][i], i) for i in 1:length(x_endpoints[1]))
    y_dict = Dict((y_endpoints[1][i], i) for i in 1:length(y_endpoints[1]))
    # Buffers for regions
    found_rectangles = Vector{OpenRectanglePrimitive}(undef, x_len * y_len)
    points = Vector{SpecialPointPrimitive}(undef, x_len * y_len)
    lines = Vector{OpenLinePrimitive}(undef, 2 * x_len * y_len)

    counter = 1
    for y_iter in 1:y_len
        for x_iter in 1:x_len
            new_poly_points = [number_of_vertices, y_endpoints[2][y_iter], x_endpoints[2][x_iter], 0]
            points[counter] = SpecialPointPrimitive((x_endpoints[1][x_iter], y_endpoints[1][y_iter]), new_poly_points)
            """Start finding growth functions for adjacent regions """
            # Vertices polynomials
            poly_points_horisontal_line = [number_of_vertices, y_endpoints[2][y_iter], 0, 0]
            poly_points_vertical_line = [number_of_vertices, 0, x_endpoints[2][x_iter], 0]
            poly_points_rect = [number_of_vertices, 0, 0, 0]
            """ Based on growth functions in 0-dimmensional region find polynomials for adjacent regions (verical and horsiontal line starting in this point
            and rectangle for which the point is lower-left corner)"""
            # Construct regions and assign found growth functions
            add_regions_0d(x_iter, y_iter, lines, found_rectangles, counter, x_endpoints, y_endpoints,
                poly_points_rect,
                poly_points_horisontal_line,
                poly_points_vertical_line)
            counter += 1
        end
    end

    # Add corner modifier to vertice polynomials (0d regions only). (If frames starts at vertice we need to add +1 to point polynomial)
    for v in vertices
        points[x_dict[v[1]] + (y_dict[v[2]] - 1) * y_len].growth_f[4] = 1
    end

    return [found_rectangles, lines, points]
end

function save_primitive_diagram(file_name, rectangles, lines, points)
    open(file_name, "w") do file
        write(file, "h-r corner;l-l corner;growth function\n")
        for rect in rectangles
            write(file, to_string(rect))
            write(file, "\n")
        end
        write(file, "\n")
        write(file, "start;end;growth function\n")
        for l in lines
            write(file, to_string(l))
            write(file, "\n")
        end
        write(file, "\n")
        write(file, "point;growth function\n")
        for p in points
            write(file, to_string(p))
            write(file, "\n")
        end
    end
end