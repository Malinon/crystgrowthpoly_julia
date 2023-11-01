# Goal of functions in this file is to make our data less tendentious
# We do not want rotation centers in the same places every time
#

function mantisse(num)
    return num - floor(num)
end

function translatePointModulo(point, trans_vec)
    return (mantisse(point[1] - trans_vec[1]), mantisse(point[2] - trans_vec[2]))
end

function translateVerticeModulo(vertice, trans_vec)
    println(translatePointModulo(vertice.point, trans_vec))
    return SpecialPointPrimitive(translatePointModulo(vertice.point, trans_vec), vertice.growth_f)
end

function translateEdgeModulo(edge, trans_vec)
    newStart = translatePointModulo(edge.start, trans_vec)
    return OpenLinePrimitive(newStart, (newStart[1] + edge.end_point[1] -edge.start[1],
    newStart[2] + edge.end_point[2] - edge.start[2]), edge.growth_f)
end

function translateRectangleModulo(rectangle, trans_vec)
    newStart = translatePointModulo(rectangle.lower_left_corner, trans_vec)
    return OpenRectanglePrimitive((newStart[1] + rectangle.higher_right_corner[1] -rectangle.lower_left_corner[1],
    newStart[2] + rectangle.higher_right_corner[2] -rectangle.lower_left_corner[2]), newStart, rectangle.growth_f)
end

function translateRegionsByRandomCorner(rectangles, edges, points)
    # Draw one corner of random rectangle
    new_origin = rectangles[mod(rand(Int), length(rectangles)) + 1].lower_left_corner
    # This new corner defines new origin
    # Translate all regions modulo by - origin
    return [[translateRectangleModulo(r, new_origin) for r in rectangles],
        [translateEdgeModulo(e, new_origin) for e in edges],
        [translateVerticeModulo(p, new_origin) for p in points]]
end