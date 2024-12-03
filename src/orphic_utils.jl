import Base:-, +, getindex

struct CellDescription
    """Structure describing cell"""
    length_in_directions::Int64
    min_coordinate_id::Int64
    max_coordinate_id::Int64
end

struct ModifierChange
    """Structure describing change of modifier"""
    modifier_change::Int64
    is_change_temporary::Bool # True if new value is valid only on single line
    is_change_immediate::Bool # True if new value is valid also at this line - closed interval
    cell_description::CellDescription
end

@enum Direction X=1 Y=2
Base.to_index(d::Direction) = Int(d)
-(v::Int, d::Direction) = v - Int(d)
+(v::Int, d::Direction) = v + Int(d)
Base.getindex(t::Tuple, d::Direction) = t[Int(d)]

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
function get_coefs_poly_point(descriptions_x::Vector{CellDescription}, descriptions_y::Vector{CellDescription},
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

function calculate_polynomial_in_next_point(reference_polynomial::Vector{Int64}, modifiers_prev::Vector{ModifierChange},
        modifiers_next::Vector{ModifierChange}, direction::Direction, coord_id::Int64)
    """ Compute growth function for 0-dimensional region based on growth function for adjacent 0-dimensional region"""
    output_polynomial = copy(reference_polynomial)
    if direction == X && coord_id ==2
        println("Start: ", output_polynomial)
    end

    for i in 1:length(modifiers_prev)
        if modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            if direction == X && coord_id ==2
                println("Prev - :",perpendicular_direction_modifier, modifiers_prev[i].modifier_change)
                println(modifiers_prev[i].cell_description)
            end
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            if direction == X && coord_id ==2
                println("Prev + :",perpendicular_direction_modifier, modifiers_prev[i].modifier_change)
            end
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end

    for i in 1:length(modifiers_next)
        if modifiers_next[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_next[i].cell_description, coord_id)
            if direction == X && coord_id ==2
                println("Next + :",perpendicular_direction_modifier, modifiers_prev[i].modifier_change)
            end
            update_polynomial!(output_polynomial, modifiers_next[i].modifier_change, perpendicular_direction_modifier, direction)

        end
    end
    if direction == X && coord_id ==2
        println("Stop: ", output_polynomial)
    end
    return output_polynomial
end

function calculate_function_for_edge(reference_polynomial::Vector{Int64},
    modifiers_prev::Vector{ModifierChange}, direction::Direction, coord_id::Int64)
    """Compute growth functions for 1d/2d region based on adjacent vertice"""
    output_polynomial = copy(reference_polynomial)
    #println("Case start, ", output_polynomial)
    for i in 1:length(modifiers_prev)
        if modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            #println("Prev - :",perpendicular_direction_modifier, modifiers_prev[i].modifier_change)
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_at_point(modifiers_prev[i].cell_description, coord_id)
            #println("Prev + :",perpendicular_direction_modifier, modifiers_prev[i].modifier_change)
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, direction)
        end
    end
    println("Case end, ", output_polynomial)
    return output_polynomial
end

function calculate_function_for_rectangle(reference_polynomial::Vector{Int64},
    modifiers_prev::Vector{ModifierChange}, coord_id::Int64)
    """Compute growth functions for 1d/2d region based on adjacent vertice"""
    output_polynomial = copy(reference_polynomial)
    for i in 1:length(modifiers_prev)
        if modifiers_prev[i].is_change_temporary
            perpendicular_direction_modifier = get_modifier_in_domain(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, -modifiers_prev[i].modifier_change, perpendicular_direction_modifier, Y)
        elseif !modifiers_prev[i].is_change_immediate
            perpendicular_direction_modifier = get_modifier_in_domain(modifiers_prev[i].cell_description, coord_id)
            update_polynomial!(output_polynomial, modifiers_prev[i].modifier_change, perpendicular_direction_modifier, Y)
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


function describe_cell(cell, dicts::Tuple{Dict{A, Int64}, Dict{A, Int64}}, direction::Direction) where A
    """Describe cell in given direction"""
    min_c, max_c = cryst_private.get_min_max(cell, direction)
    lengths_in_directions = floor(Int64, max_c) - floor(Int64, min_c)
    min_coordinate_id = dicts[direction][cryst_private.mantisse(min_c)]
    max_coordinate_id = dicts[direction][cryst_private.mantisse(max_c)]
    return CellDescription(lengths_in_directions, min_coordinate_id, max_coordinate_id), min_coordinate_id, max_coordinate_id
end

function create_modifier_change!(buffer::Vector{Vector{ModifierChange}}, min_coordinate_id::Int64, max_coordinate_id::Int64,
    perpendicular_cell_description::CellDescription)
    """Create modifier change object"""
    if min_coordinate_id != max_coordinate_id
        # if {min_x} < x_0 , then modifier += 1.
        push!(buffer[min_coordinate_id], ModifierChange(1, false, false, perpendicular_cell_description))
        # if {max_x} <= x_0 , then modifier += 1.
        push!(buffer[max_coordinate_id], ModifierChange(-1, false, true, perpendicular_cell_description))
    else
            # Min and max coordinates has the same mantissa, so we have only one -1 pick of modifier
        push!(buffer[min_coordinate_id], ModifierChange(-1, true, true, perpendicular_cell_description))
    end
end

"""
Based on 1>=0 dim. cells describe changes in modifiers and create compressed representation of cells.


# Arguments
- `cells`: Edges and faces of the tessellation.
- `lengths::Tuple{Int64, Int64}`: Numbers of unique coordinates mantissas in each direction.
- `dicts::Tuple{Dict{A, Int64}, Dict{A, Int64}}`: Dictionaries mapping mantissas to unique coordinate ids in each direction.

# Returns
- `modifiers`: A tuple of vectors of vectors of `ModifierChange` objects, representing the changes in modifiers for each direction.
- `fig_descriptions`: A tuple of vectors of `CellDescription` objects, representing the descriptions of the cells for each direction.

# Description
This function initializes and populates arrays for modifiers and cell descriptions. It iterates over each cell and each direction, calculates the minimum and maximum coordinates, and creates descriptions and modifiers based on these coordinates. The function returns the populated arrays of modifiers and cell descriptions.
"""
function describe_change_points_and_cells(cells, lengths::Tuple{Int64, Int64},
    dicts::Tuple{Dict{A, Int64}, Dict{A, Int64}}) where A
    """Describe changes of modifiers and cells"""

    # Initialize the modifiers and fig_descriptions arrays
    modifiers = (Vector{Vector{ModifierChange}}(undef, lengths[1]), Vector{Vector{ModifierChange}}(undef, lengths[2]))
    for i in 1:lengths[1]
        modifiers[1][i] = Vector{ModifierChange}()
    end
    for i in 1:lengths[2]
        modifiers[2][i] = Vector{ModifierChange}()
    end

    fig_descriptions = (Vector{CellDescription}(undef, length(cells)), Vector{CellDescription}(undef, length(cells)))

    for i in 1:length(cells)
        # Create description of cell
        fig_descriptions[X][i], x_min_coordinate_id, x_max_coordinate_id = describe_cell(cells[i], dicts, X)
        fig_descriptions[Y][i], y_min_coordinate_id, y_max_coordinate_id = describe_cell(cells[i], dicts, Y)

        # Create modifier changes for cell (with description of cell in perpendicular direction)
        create_modifier_change!(modifiers[X], x_min_coordinate_id, x_max_coordinate_id, fig_descriptions[Y][i])
        create_modifier_change!(modifiers[Y], y_min_coordinate_id, y_max_coordinate_id, fig_descriptions[X][i])
    end

    return modifiers, fig_descriptions
end
