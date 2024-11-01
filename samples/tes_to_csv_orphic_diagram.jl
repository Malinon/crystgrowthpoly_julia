import Pkg
import crystgrowthpoly

using DataFrames
using DelimitedFiles
using Nemo

function create_data_frame(rectangles, edges, vertices)
    df = DataFrame(Dimmension = Int64[],
    Start_x = Float64[], Start_y = Float64[], End_x = Float64[], End_y = Float64[],
    Growth_function_vertices_0 = Int64[], Growth_function_vertices_1= Int64[], Growth_function_vertices_2=Int64[], Growth_function_vertices_3=Int64[],
    Growth_function_edges_0=Int64[], Growth_function_edges_1=Int64[], Growth_function_edges_2=Int64[], Growth_function_edges_3=Int64[],
    Growth_function_faces_0=Int64[], Growth_function_faces_1=Int64[], Growth_function_faces_2=Int64[], Growth_function_faces_3=Int64[],
    )
    add_points_to_df!(df, vertices)
    add_edges_to_df!(df, edges)
    add_rectangles_to_df!(df, rectangles)
    return df
end

function save_df_to_files(out_name, vertices, edges, rectangles)
    df = create_data_frame(vertices, edges, rectangles)
    open(out_name, "w") do f
        # Step 3: Write the header (column names)
        writedlm(f, [names(df)], ',')
    
        # Step 4: Write each row of the DataFrame
        for row in eachrow(df)
            writedlm(f, [collect(row)], ',')
        end
    end
end

function add_points_to_df!(df, vertices)
    for v in vertices
        push!(df, (0,
        Float64(v.point[1]), Float64(v.point[2]), Float64(v.point[1]), Float64(v.point[2]),
        v.growth_f[1]..., v.growth_f[2]..., v.growth_f[3]...))
    end
end

function add_edges_to_df!(df, edges)
    for e in edges
        push!(df, (1,
        Float64(e.start[1]), Float64(e.start[2]), Float64(e.end_point[1]), Float64(e.end_point[2]),
        e.growth_f[1]..., e.growth_f[2]..., e.growth_f[3]...))
    end
end

function add_rectangles_to_df!(df, rectangles)
    for r in rectangles
        push!(df, (2,
        Float64(r.lower_left_corner[1]), Float64(r.lower_left_corner[2]), Float64(r.higher_right_corner[1]), Float64(r.higher_right_corner[2]),
        r.growth_f[1]..., r.growth_f[2]..., r.growth_f[3]...))
    end
end


const INPUT_PATH = ARGS[1]
const OUTPUT_PATH = ARGS[2]

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

tessellation = crystgrowthpoly.read_tessellation_from_file(INPUT_PATH)
reduced_tessellation = crystgrowthpoly.Tessellation(crystgrowthpoly.reduce_polygon(tessellation.polygon.cells[1], tessellation.polygon.cells[2], tessellation.polygon.cells[3]))
orphic = crystgrowthpoly.find_regions(reduced_tessellation.polygon.cells[1], reduced_tessellation.polygon.cells[2], reduced_tessellation.polygon.cells[3])
save_df_to_files(OUTPUT_PATH, orphic[1], orphic[2], orphic[3])

