import Base:-

struct CellDescription
    """Structure describing cell"""
    length_in_directions::Int64
    min_coordinate_id::Int64
    max_coordinate_id::Int64
end

struct ModifierChangeDescription
    """Structure describing change of modifier"""
    modifier_change::Int64
    is_change_temporary::Bool # True if new value is valid only on single line
    is_change_immediate::Bool # True if new value is valid also at this line - closed interval
    cell_description::CellDescription
end

@enum Direction X=1 Y=2
Base.to_index(d::Direction) = Int(d)
-(v::Int, d::Direction) = v - Int(d)

function get_modifier_at_point(description::CellDescription, coord_id::Int64)
    modifier = description.length_in_directions
    if description.min_coordinate_id < coord_id
        modifier = modifier + 1
    end

    if description.max_coordinate_id <= coord_id
        modifier = modifier - 1
    end
    return modifier
end

function get_modifier_in_domain(description::CellDescription, coord_id::Int64)
    modifier = description.length_in_directions
    if description.min_coordinate_id < coord_id
        modifier = modifier + 1
    end

    if description.max_coordinate_id < coord_id
        modifier = modifier - 1
    end
    return modifier
end

""" Find growth function (1-cells, 2-cells) for 0 dimmensional region """
function get_coefs_poly_point(descriptions_x::Vector{ModifierChangeDescription}, descriptions_y::Vector{ModifierChangeDescription},
        index_x::Int64, index_y::Int64)
    number_of_cells = length(descriptions_x)
    coefficients = [number_of_cells, 0, 0, 0]::Array{Int64}
    for i in 1:number_of_cells
        x_modifier = get_modifier_at_point(descriptions_x[i], index_x)
        y_modifier = get_modifier_at_point(descriptions_y[i], index_y)
        coefficients[2] = coefficients[2] - y_modifier
        coefficients[3] = coefficients[3] - x_modifier
        coefficients[4] = coefficients[4] + x_modifier * y_modifier
    end
    return coefficients
end

function update_polynomial!(polynomial::Vector{Int64}, change::Int64, perpendicular_modifier::Int64, direction::Direction)
    polynomial[4] = polynomial[4] + perpendicular_modifier * change
    polynomial[4 - direction] -= change
end

function calculate_polynomial_in_next_point(reference_polynomial::Vector{Int64}, modifiers_prev::Vector{ModifierChangeDescription},
        modifiers_next::Vector{ModifierChangeDescription}, direction::Direction, coord_id::Int64)
    """ Compute growth function for 0-dimensional region based on growth function for adjacent 0-dimensional region"""
    output_polynomial = copy(reference_polynomial)
    for i in 1:length(modifiers_prev)
        if !modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end

    for i in 1:length(modifiers_next)
        if modifiers_next[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_next[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, modifiers_next[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end
    return output_polynomial
end

function calculate_function_for_edge(reference_polynomial::Vector{Int64},
    modifiers_prev::Vector{ModifierChangeDescription}, direction::Direction, coord_id::Int64)
    """Compute growth functions for 1d/2d region based on adjacent vertice"""
    output_polynomial = copy(reference_polynomial)
    for i in 1:length(modifiers_prev)
        if !modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end

    return output_polynomial
end

function calculate_function_for_rectangle(reference_polynomial::Vector{Int64},
    modifiers_prev::Vector{ModifierChangeDescription}, direction::Direction, coord_id::Int64)
    """Compute growth functions for 1d/2d region based on adjacent vertice"""
    output_polynomial = copy(reference_polynomial)
    for i in 1:length(modifiers_prev)
        if !modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_in_domain(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_in_domain(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end

    return output_polynomial
end

function find_boundary_lines(vertices::Vector{Tuple{A, A}}, direction::Direction) where A
    """Finds lines perpendicular to chosen axis defining domains"""
    index = direction
    unique_points_coord = Dict{A, Int64}()
    for p in vertices
        unique_points_coord[p[index]] = get(unique_points_coord, p[index], 0) + 1
    end
    endpoints = collect(keys(unique_points_coord))
    sort!(endpoints)
    vertices_at_lines = [unique_points_coord[c] for c in endpoints]
    return endpoints, vertices_at_lines
end

function describe_change_points_and_cells(cells, lengths::Tuple{Int64, Int64}, dicts::Tuple{Dict{A, Int64}, Dict{A, Int64}}) where A
    """Describe changes of modifiers and cells"""

    modifiers = Tuple(Vector{Vector{ModifierChangeDescription}}(undef, len) for len in lengths)
    fig_descriptions = Tuple(Vector{CellDescription}(undef, length(cells)) for _ in lengths)

    for i in 1:length(cells)
        for direction in instances(Direction)
            # Create description of cell
            min_c, max_c = get_min_max(cells[i], direction)
            lengths_in_directions = ceil(Int64, max_c - min_c)
            min_coordinate_id = dicts[direction][min_c]
            max_coordinate_id = dicts[direction][max_c]
            fig_descriptions[direction][i] = CellDescription(lengths_in_directions, min_coordinate_id, max_coordinate_id)

            # Create modifiers
            if min_coordinate_id != max_coordinate_id
                # if {min_x} < x_0 , then modifier += 1.
                modifiers[direction][min_coordinate_id].push(ModifierChangeDescription(1, false, false, fig_descriptions[direction]))
                # if {max_x} <= x_0 , then modifier += 1.
                modifiers[direction][max_coordinate_id].push(ModifierChangeDescription(-1, false, true, fig_descriptions[direction]))
            else
                # Min and max coordinates has the same mantissa, so we have only one -1 pick of modifier
                modifiers[direction][min_coordinate_id].push(ModifierChangeDescription(-1, true, true, fig_descriptions[direction]))
            end
        end
    end

    return modifiers, fig_descriptions
end
