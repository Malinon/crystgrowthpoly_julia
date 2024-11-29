using Nemo

a = QQ(2)

module cryst_preprocessor


function string_coord_to_coord(s)
    if occursin("/", s)
        # Rational
        nom_denom = split(s, "/", keepempty=false)
        return QQ(parse(Int64, nom_denom[1])) // parse(Int64, nom_denom[2])
    else
        # Integer
        return QQ(parse(Int64, s))
    end
end

function line_to_point(s)
    s_coords = split(s, " ", keepempty=false)
    return (string_coord_to_coord(s_coords[1]), string_coord_to_coord(s_coords[2]))
end

function ids_to_faces(verts, ids)
    return Tuple(verts[parse(Int64, id)] for id in ids)
end

end # module cryst_preprocessor