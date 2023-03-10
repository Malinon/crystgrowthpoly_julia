include("private.jl")

struct Polygon
    cells
    function Polygon(vertices, edges, faces)
        new_e = [Tuple(sort(collect(elem))) for elem in edges]
        new_f = [Tuple(sort(collect(elem))) for elem in faces]
        new((vertices, new_e, new_f))
    end
end

# TODO: Choose better type for A and B
struct Tessellation
    polygon::Polygon
end

function reduce_polygon(ver::Vector{Tuple{A,A}}, edge::Vector{Tuple{Tuple{A,A},Tuple{A,A}}}, faces) where A
    unique_vertices = Set{Tuple{A,A}}()
    for v in ver
        reduced_p = (v[1] - floor(v[1]), v[2] - floor(v[2]))
        push!(unique_vertices, reduced_p)
    end
    unique_edges = Set{Tuple{Tuple{A,A},Tuple{A,A}}}()
    for e in edge
        trans =  (-floor(e[1][1]), -floor(e[1][2]))
        push!(unique_edges, cryst_private.translate_face(e, trans))
    end
    return Polygon(collect(unique_vertices), collect(unique_edges), faces)
end