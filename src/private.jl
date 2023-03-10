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


## @param cells_in_tessellation Cells, which already are in tessellation
## @param cells All repeating figure's cells of some dimmension
## @param vec Translation vector
## @param move_operator Function translating given cell by vec
## @param censor  Bool function returning True if cell should be counted, False otherwise
function add_cells_censored(cells_in_tessellation, cells, vec, move_operator, censor)
    """Add cells to tessellation model and count number of new cells which satisfy given condition"""
    for c in cells
        cell = move_operator(c, vec)
        if censor(cell)
            push!(cells_in_tessellation, cell)
        end
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


## @param cell Tuple representing cell
## @returns Sorted tuple representing input cell !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sort_points_in_cell(cell)
    """Sorts points representing cell"""
    function is_grater_then(tup1, tup2)
        for i in 1:len(tup1)
            if tup1[i] > tup2[i]
                return true
            end
        end
        return false
    end
    cell_to_list = list(cell)
    cell_to_list.sort()
    return Tuple(cell_to_list)
end


## @param n Length of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
## @param m Lenght of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
## @param scale_v1 Scaling value
## @param scale_v2 Scaling value
## @param x0 Anchor point of parallelogram frame
function gen_is_point_in_parallelogram(n::Int, m::Int, scale_v1, scale_v2, x0)
    """Returns function checking if point is inside parallelogram"""
    limiter_n = n * scale_v1
    limiter_m = m * scale_v2
    function is_point_in_parallelogram(point)
        vec  = subtract_vectors(point, x0)
        return  (0 <= vec[1] && 0 <= vec[2] && vec[1] <= limiter_n && vec[2] <= limiter_m)
    end
    return  is_point_in_parallelogram
end


## @param n Length of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
## @param m Lenght of rectangle's side before scaling (In crystallographical coordinates parallelogram frames used by us are rectangles)
## @param scale_v1 Scaling value
## @param scale_v2 Scaling value
## @param x0 Anchor point of parallelogram frame !!!!!!!!!Capture
function gen_is_face_in_parallelogram(n::Int, m::Int, scale_v1, scale_v2, x0)
    """Returns function checking if cell (with dimmension higher then 1) is inside parallelogram"""
    is_point_in_parallelogram = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
    return (face -> all(is_point_in_parallelogram(point) for point in face))
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

## @param n Number of translations in v1 direction
## @param m Number of translations in v2 direction
## @param cells List of k-cells in repeating unit
## @param v1  Translation vector
## @param v2  Translation vector
## @param k Dimmension of cell
## @param x0 Anchor point of parallelogram frame
## @param scale_v1 Scaling value for vector v1
## @param scale_v2 Scaling value for vector v2
## @param additional_limits numbers of additional translations that need to be performed
## @return Number of k-cells inside parallelogram
function get_k_cells_num_parallelogram(n::Int, m::Int, cells, v1, v2, k, x0, scale_v1, scale_v2, additional_limits)
    """Function calculating number of cells with given dimmension inside parallelogram"""
    cells_in_tessellation = Set()
    if k == 0
        # Count 0 cells
        censor = gen_is_point_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_vector
    else
        # Count higher dimmension cells
        censor = gen_is_face_in_parallelogram(n, m, scale_v1, scale_v2, x0)
        move_operator = translate_face
    end
    for i in -additional_limits[1][1]:((1 + n + additional_limits[1][2]) * ceil(scale_v1))
        for j in -additional_limits[2][1] : ((1 + m + additional_limits[2][2]) * ceil(scale_v2))
            add_cells_censored(cells_in_tessellation, cells, multiply_by_scalar_and_add([v1, v2], [i, j]), move_operator, censor)
        end
    end
    return length(cells_in_tessellation)
 end

function my_numerator(x::ca)
    return numerator(QQ(x))
end

function my_numerator(x)
    return numerator(x)
end
## @param x0 Anchor point of the frame
## @param points Vertices in repeating motif
function get_limits_extenders(x0, points)
    """This function calculate how many additional translations are needed to be performed."""
    max_v1 = 0
    min_v1 = 0
    max_v2 = 0
    min_v2 = 0
    # Find extreme coordinates of repeating unit
    for p in points
        if p[1] - x0[1] > max_v1
            max_v1 = p[1] - x0[1]
        elseif p[1] - x0[1] < min_v1
            min_v1 = p[1] - x0[1]
            if floor(min_v1) == min_v1
                min_v1 = min_v1 -1
            end
        end
        if p[2] - x0[2] > max_v2
            max_v2 = p[2] - x0[2]
        elseif p[2] - x0[2] < min_v2
            min_v2 = p[2] - x0[2]
            if floor(min_v2) == min_v2
                min_v2 = min_v2 - 1
            end
        end
    end

    #if fits(Int64, numerator(floor(max_v1))) && fits(Int64, numerator(-floor(min_v1))) && fits(Int64, numerator(floor(max_v2))) && fits(Int64, numerator(-floor(min_v2)))
    return ((Int64(my_numerator(floor(max_v1))), Int64(my_numerator(-floor(min_v1)))), (Int64(my_numerator(floor(max_v2))), Int64(my_numerator(-floor(min_v2)))))
    # TODO: ADD Exception
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


function get_faces_polynomial(edges, x0, symmetric_growth::Bool=false)
    coefficients = Any[0,0,0,0]
    edgy_k_counter = 0::Int
    edgy_l_counter = 0::Int
    edgy_both_counter = 0::Int
    for e in edges
        face_in_new_system = translate_face(e, (-x0[1], -x0[2]))
        minmax_k = get_min_max(face_in_new_system, 1)
        minmax_l = get_min_max(face_in_new_system, 2)

        floored_max_k = (floor(minmax_k[2]))
        floored_max_l = (floor(minmax_l[2]))
        if minmax_k[2] == floored_max_k
            if minmax_l[2] == floored_max_l
                edgy_both_counter = edgy_both_counter + 1
            else
                edgy_k_counter = edgy_k_counter + 1
            end
        elseif minmax_l[2] == floored_max_l
            edgy_l_counter = edgy_l_counter + 1
        end
        mod_k = floored_max_k - (floor(minmax_k[1]))
        mod_l = floored_max_l - (floor(minmax_l[1]))
        coefficients[1] = coefficients[1] + 1
        coefficients[2]  = coefficients[2] + 1 - mod_l
        coefficients[3]  = coefficients[3] + 1 - mod_k
        coefficients[4] =  coefficients[3] + 1 - mod_l - mod_k - mod_k * mod_l
    end
    coefficients[2] = coefficients[2] + edgy_l_counter + edgy_both_counter
    coefficients[3] = coefficients[3] + edgy_k_counter + edgy_both_counter
    coefficients[4] = coefficients[4] - edgy_both_counter
    return coefficients
end
end


