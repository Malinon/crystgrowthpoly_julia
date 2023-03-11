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

"""This class represents 2-dimmensional domains"""
struct OpenRectanglePrimitive{A}
    higher_right_corner::Tuple{A,A}
    lower_left_corner::Tuple{A,A}
    growth_f::Vector{Int64}
end

mutable struct OpenLinePrimitive{A}
    start::Tuple{A,A}
    end_point::Tuple{A,A}
    growth_f::Vector{Int64}
end

mutable struct SpecialPointPrimitive{A}
    """This class represent 0-dimmensional domains"""
    point::Tuple{A,A}
    growth_f::Vector{Int64}
end

const SEPERATOR = ";"

function to_string(rect::OpenRectanglePrimitive{A}) where A
    return string(rect.higher_right_corner) * SEPERATOR * string(rect.lower_left_corner) * SEPERATOR  * string(rect.growth_f)
end

function to_string(line::OpenLinePrimitive{A}) where A
    return string(line.start) * SEPERATOR * string(line.end_point) * SEPERATOR  * string(line.growth_f)
end

function to_string(point::SpecialPointPrimitive{A}) where A
    return string(point.point) * SEPERATOR * string(point.growth_f)
end
