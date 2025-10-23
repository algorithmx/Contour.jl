module TestHelpers

using Contour, Test, LinearAlgebra

export all_cell_configurations, create_test_matrix_from_config, generate_test_grid,
       simple_curvilinear_grid, verify_interpolation_accuracy, count_contour_vertices,
       is_closed_contour

"""
    generate_test_grid(nx=10, ny=10; func=(x,y) -> x^2 + y^2)

Generate a test grid with specified dimensions and function.
Returns (x, y, z) where x and y are coordinate arrays and z is the scalar field.
"""
function generate_test_grid(nx=10, ny=10; func=(x,y) -> x^2 + y^2)
    x = range(-1, 1, length=nx)
    y = range(-1, 1, length=ny)
    z = [func(xi, yi) for xi in x, yi in y]
    return x, y, z
end

"""
    all_cell_configurations()

Generate test cases for all 16 cell types.
Returns a vector of (case_number, z_values) tuples.
"""
function all_cell_configurations()
    configs = []
    for case in 0:15
        z_values = [case & 0x01 > 0, case & 0x02 > 0,
                    case & 0x04 > 0, case & 0x08 > 0]
        # Create values around level 0.5: below=0.0, above=1.0
        values = [z ? 1.0 : 0.0 for z in z_values]
        push!(configs, (case, values))
    end
    return configs
end

"""
    create_test_matrix_from_config(z_values, level=0.5)

Create a 2x2 test matrix from 4 vertex values relative to a contour level.
"""
function create_test_matrix_from_config(z_values, level=0.5)
    # z_values should be [z1, z2, z3, z4] corresponding to vertices 1-4
    # z_values are booleans indicating above/below level
    scaled_values = [z ? level + 0.5 : level - 0.5 for z in z_values]
    return reshape(scaled_values, 2, 2)
end

"""
    simple_curvilinear_grid(nr=10, nθ=20)

Create a simple polar coordinate grid for testing curvilinear interpolation.
Returns (x, y, z) where x and y are coordinate matrices and z is the radius.
"""
function simple_curvilinear_grid(nr=10, nθ=20)
    r = range(1.0, 2.0, length=nr)
    θ = range(0.0, 2π, length=nθ)

    x = [r_i * cos(θ_j) for r_i in r, θ_j in θ]
    y = [r_i * sin(θ_j) for r_i in r, θ_j in θ]
    z = sqrt.(x.^2 + y.^2)  # radius values

    return x, y, z
end

"""
    verify_interpolation_accuracy(result, expected, tolerance=1e-6)

Verify that interpolation result matches expected values within tolerance.
"""
function verify_interpolation_accuracy(result, expected, tolerance=1e-6)
    if isa(result, Tuple)
        return all(isapprox(r, e, atol=tolerance) for (r, e) in zip(result, expected))
    else
        return isapprox(result, expected, atol=tolerance)
    end
end

"""
    count_contour_vertices(contour_level)

Count total number of vertices across all lines in a contour level.
"""
function count_contour_vertices(contour_level)
    return sum(length(line.vertices) for line in contour_level.lines)
end

"""
    is_closed_contour(line)

Check if a contour line forms a closed loop.
"""
function is_closed_contour(line)
    if length(line.vertices) < 2
        return false
    end
    first_vertex = first(line.vertices)
    last_vertex = last(line.vertices)
    return isapprox(first_vertex[1], last_vertex[1], atol=1e-10) &&
           isapprox(first_vertex[2], last_vertex[2], atol=1e-10)
end

end # TestHelpers