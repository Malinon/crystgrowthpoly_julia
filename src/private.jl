module cryst_private
using Nemo
import IterTools

function translate_vector(tup, vec)
    return Tuple(tup[i] + vec[i] for i in 1:length(tup))
end

## @param face Cell to translate
## @param vec Translation vector
function translate_face(face, vec)
    return Tuple(translate_vector(face[i], vec) for i in 1:length(face))
end

## @param cells_in_tessellation Cells, which already are in tessellation
## @param cells All repeating figure's cells of some dimmension
## @param vec Translation vector
## @param move_operator Function translating given cell by vec
function add_cells(cells_in_tessellation, cells, vec, move_operator)
    """Translate and add cells to tessellation model"""
    for c in cells
        push!(cells_in_tessellation, move_operator(c, vec))
    end
end

function multiply_vector(vec, scalar)
    """Returns vector multiplied by scalar"""
    return Tuple(a * scalar for a in vec)
end

## @param vectors List of vectors
## @param scalars List of scalars
function multiply_by_scalar_and_add(vectors, scalars)
    """Multiplies vectors by scalars and add them."""
    dim = length(vectors[1])
    return Tuple(sum(vectors[j][i] * scalars[j] for j in 1:length(vectors)) for i in 1:dim)
end

function subtract_vectors(minuend, subtrahend)
    """Return difference of two vectors"""
    return Tuple(z[1] - z[2] for z in zip(minuend, subtrahend))
end

function generate_matrix(args, symmetric_growth)
    """Creates matrix used for growth polynomial interpolation"""
    if symmetric_growth
        n = length(args) - 1
        return matrix(QQ, length(args), length(args), [arg[1]^(n - i) for arg in args, i in 0:n]) # Cache matrix 
    else
        dim = length(args[1])
        indexes  = 1:dim
        row_scheme = Tuple(itertools.chain(*(Tuple(IterTools.subsets(indexes, dim - n)) for n in 1:(dim+1)))) # TODO: Optimize
        function generate_row(arg)
            return Tuple(functools.reduce((acc, index) ->  acc * arg[index], chosen_indexes, 1) for chosen_indexes in row_scheme)
        end
        return matrix(QQ, length(args), length(args), [generate_row(arg) for arg in args])
    end
end

## @param gen_val Function generating value of polynomial in given points
## @param args Data points
## @param symmetric_growth True, if growth is the same in each direction, False otherwise.
##
## @return Coefficients of polynomial
function find_poly(gen_val, args, symmetric_growth::Bool)
    value_vector = matrix(QQ, length(args), 1, [gen_val(arg) for arg in args])
    mat = generate_matrix(args, symmetric_growth)
    return solve(mat, value_vector)
end


## @param arguments Numbers of steps
## @param cells List of cells with given dimmension in repeating unit
## @param translation_vectors List of translation vectors
## @return Number of cells with some dimmension after arguments[i] steps in direction translation_vectors[i]
function get_k_cells_num(arguments, cells, translation_vectors, k)
    """Function calculating number of cells with some dimmension after arguments[i] steps in direction translation_vectors[i]"""
    cells_in_tessellation = Set()
    if k == 0
        # Count 0 cells
        move_operator = translate_vector
    else
        # Count higher dimmension cells
        move_operator = translate_face
    end
    for translations_numbers in IterTools.product((0:(arg-1) for arg in arguments)...)
        add_cells(cells_in_tessellation, cells, multiply_by_scalar_and_add(translation_vectors, translations_numbers), move_operator)
    end
    return length(cells_in_tessellation)
end

## @param coefficient_lists List of lists containing polynomials' coefficients
## @return List of alternative sums
function get_alternative_sum(coefficient_lists, dim)
    """This function calculates alternative sums of corresponding coefficients."""
    len = length(coefficient_lists[1])
    return_list = Vector{fmpq}(undef, len)
    for i in 1:len
        acc = 0
        for j  in 1:(dim - 1)
            acc = acc - coefficient_lists[j][i] * ((-1) ^ (j % 2))
        end
        return_list[i] = acc
    end
    return return_list
end

## @param coefficient_lists List containing lists of growth polynomial coefficients corresponding to a cells with dimensions 0 through dim - 2 
## @param max_dimmension_coefficients Coefficients of growth polynomial corresponding to cells of dimmension dim
## @param dim Dimmension of tessellation
## @return Coefficients of growth polynomial corresponding to cells with dimmension dim - 1
function calculate_coefficients_based_ne_euler_charcteristics(coefficient_lists, max_dimmension_coefficients, dim)
    """Calculates coefficients of polynomial corresponding to a cell with dimension dim - 1"""
    return_list = get_alternative_sum(coefficient_lists, dim)
    return_list[length(return_list)] = return_list[length(return_list)] - 1
    if dim % 2 == 0
        for i in 1:length(return_list)
            return_list[i] = return_list[i] + max_dimmension_coefficients[i]
        end
    else
        for i in 1:length(return_list)
            return_list[i] = -(return_list[i] - max_dimmension_coefficients[i])
        end
    end
    return matrix(QQ, length(return_list), 1, return_list)
end


function get_best_realisation(v, x0)
    centralized_x  = v[1] - x0[1]
    centralized_y  = v[2] - x0[2]
    return (centralized_x - floor(centralized_x), centralized_y - floor(centralized_y))
end

function get_vertices_polynomial(vertices, x0, symmetric_growth::Bool=false)
    counter_x0 = 0
    counter_k = 0
    counter_l = 0
    for v in vertices
        new_p = get_best_realisation(v, x0)
        if new_p[1] == 0
            if new_p[2] == 0
                counter_x0 = counter_x0 + 1
            else
                counter_k = counter_k + 1
            end
        elseif new_p[2] == 0
            counter_l = counter_l + 1
        end
    end
    return (length(vertices), counter_l + counter_x0, counter_k + counter_x0, counter_x0)
end


function get_min_max(lis, index)
    max_val = lis[1][index]
    min_val = lis[1][index]
    for val in lis
        if max_val < val[index]
            max_val = val[index]
        elseif min_val > val[index]
            min_val = val[index]
        end
    end
    return (min_val, max_val)
end

function mantisse(x)
    return x - floor(x)
end

# 
# Function: get_modifier_for_cell
# 
# Description:
# This function calculates a modifier for a given cell based on its coordinates and a specified direction. 
# Impact of cell is equal I_(cell) = PRODUCT (n - modifier)
# 
# Parameters:
# - cell: The cell for which the modifier is being calculated.
# - x0_mantisse: A mantissa of selected coordinate of the anchor point
# - direction: The direction for modifier is calculated
# 
# Returns:
# - modifier: An integer value representing the calculated modifier for the cell.
function get_modifier_for_cell(cell, x0_mantisse, direction)
    min_coordinate, max_coordinate = get_min_max(cell, direction)
    floored_max = floor(max_coordinate)
    floored_min = floor(min_coordinate)

    modifier = floored_max - floored_min
    if ((max_coordinate - floored_max) <= x0_mantisse)
        modifier -= 1
    end
    if ((min_coordinate - floored_min) < x0_mantisse)
        modifier += 1
    end

    return modifier
end

"""
    get_high_dim_cell_polynomial(cells, x0)

Calculate the coefficients of a polynomial based on the provided cells and initial value `x0`.

# Arguments
- `cells::Vector`: A vector of cells reprersented by their vertices in cryst. coordinates
- `x0::Tuple`: Anchor point of frame

# Returns
- `coefficients::Vector{Any}`: A vector containing the coefficients of the polynomial in the form `[a_1, a_2, a_3, a_4]`, where:
  - f(k,l) = a_1 * kl + a_2 * k + a_3 * l + a_4
"""
function get_high_dim_cell_polynomial(cells, x0)
    coefficients = Any[0,0,0,0] # Polynomial a_1 * kl + a_2 * k + a_3 * l + a_4 
    x0_mantisse = mantisse(x0)
    for cell in cells
        modifier_k = get_modifier_for_cell(cell, x0_mantisse, 1)
        modifier_l = get_modifier_for_cell(cell, x0_mantisse, 2)
        coefficients[4] += modifier_k * modifier_l
        coefficients[3] -= modifier_k
        coefficients[2] -= modifier_l
    end
    coefficients[1] = len(edges)
    return coefficients
end

end # module cryst_private
