"""This class represents 2-dimmensional domains"""
struct OpenRectangle{A}
    higher_right_corner::Tuple{A,A}
    lower_left_corner::Tuple{A,A}
    growth_f::Tuple{ Vector{Int64}, Vector{Int64}, Vector{Int64}}
end

mutable struct OpenLine{A}
    start::Tuple{A,A}
    end_point::Tuple{A,A}
    growth_f::Tuple{ Vector{Int64}, Vector{Int64}, Vector{Int64}}
end

mutable struct SpecialPoint{A}
    """This class represent 0-dimmensional domains"""
    point::Tuple{A,A}
    growth_f::Tuple{ Vector{Int64}, Vector{Int64}, Vector{Int64}}
end
