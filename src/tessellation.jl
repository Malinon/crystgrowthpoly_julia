"""
    struct Tessellation

A structure representing a tessellation, which consists of cells defined by vertices, edges, and faces.

# Fields
- `cells`: A tuple containing the vertices, edges, and faces of the tessellation.

# Inner Constructor
- `Polygon(vertices, edges, faces)`: Constructs a `Tessellation` object. The edges and faces are sorted and converted to tuples before being stored in the `cells` field.

    ## Arguments
    - `vertices`: A collection of vertices.
    - `edges`: A collection of edges, where each edge is represented as a collection of vertices.
    - `faces`: A collection of faces, where each face is represented as a collection of edges.
"""
struct Tessellation
    cells
    function Polygon(vertices, edges, faces)
        new_e = [Tuple(sort(collect(elem))) for elem in edges]
        new_f = [Tuple(sort(collect(elem))) for elem in faces]
        new((vertices, new_e, new_f))
    end
end


"""
    reduce_polygon(ver::Vector{Tuple{A,A}}, edge::Vector{Tuple{Tuple{A,A},Tuple{A,A}}}, faces) where A

Reduces the vertices and edges of a polygon by translating them (or their minimum coordinates) to a unit square.

# Arguments
- `ver::Vector{Tuple{A,A}}`: A vector of vertex coordinates, where each vertex is represented as a tuple of two coordinates.
- `edge::Vector{Tuple{Tuple{A,A},Tuple{A,A}}}`: A vector of edges, where each edge is represented as a tuple of two vertex tuples.
- `faces`: Additional face information (type unspecified).

# Returns
- `Tessellation`: A `Tessellation` object containing the reduced vertices and edges, and the original faces.

# Details
- The function translates each vertex to a unit square by subtracting the floor of its coordinates.
- Each edge is translated by the same amount as its starting vertex.
- The resulting vertices and edges are collected into sets to ensure uniqueness before being converted back to vectors.
"""
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
    return Tessellation(collect(unique_vertices), collect(unique_edges), faces)
end
