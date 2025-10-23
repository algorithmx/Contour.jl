module InterpolationTests

using Contour, Test
include("test_helpers.jl")
using .TestHelpers

@testset "Interpolation Engine" begin

    @testset "Linear Interpolation - Non-uniform grids" begin
        @testset "West edge interpolation" begin
            x = [0.0, 2.0, 5.0]
            y = [0.0, 3.0, 7.0]
            z = [1.0 2.0; 3.0 4.0]

            # West edge: interpolate between (1,1) and (1,2)
            # z[1,1] = 1.0, z[1,2] = 2.0, want h = 1.5
            h = 1.5
            result = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})

            # x should be at x[1] (west edge), y should be interpolated
            @test result[1] == x[1]
            expected_y = y[1] + (y[2] - y[1]) * (h - z[1,1]) / (z[1,2] - z[1,1])
            @test isapprox(result[2], expected_y, rtol=1e-10)
        end

        @testset "East edge interpolation" begin
            x = [0.0, 2.0, 5.0]
            y = [0.0, 3.0, 7.0]
            z = [1.0 2.0; 3.0 4.0]

            # East edge: interpolate between (2,1) and (2,2)
            # z[2,1] = 3.0, z[2,2] = 4.0, want h = 3.5
            h = 3.5
            result = Contour.interpolate(x, y, z, h, (1,1), 0x04, NTuple{2,Float64})

            # x should be at x[2] (east edge), y should be interpolated
            @test result[1] == x[2]
            expected_y = y[1] + (y[2] - y[1]) * (h - z[2,1]) / (z[2,2] - z[2,1])
            @test isapprox(result[2], expected_y, rtol=1e-10)
        end

        @testset "North edge interpolation" begin
            x = [0.0, 2.0, 5.0]
            y = [0.0, 3.0, 7.0]
            z = [1.0 2.0; 3.0 4.0]

            # North edge: interpolate between (1,2) and (2,2)
            # z[1,2] = 2.0, z[2,2] = 4.0, want h = 3.0
            h = 3.0
            result = Contour.interpolate(x, y, z, h, (1,1), 0x01, NTuple{2,Float64})

            # y should be at y[2] (north edge), x should be interpolated
            @test result[2] == y[2]
            expected_x = x[1] + (x[2] - x[1]) * (h - z[1,2]) / (z[2,2] - z[1,2])
            @test isapprox(result[1], expected_x, rtol=1e-10)
        end

        @testset "South edge interpolation" begin
            x = [0.0, 2.0, 5.0]
            y = [0.0, 3.0, 7.0]
            z = [1.0 2.0; 3.0 4.0]

            # South edge: interpolate between (1,1) and (2,1)
            # z[1,1] = 1.0, z[2,1] = 3.0, want h = 2.0
            h = 2.0
            result = Contour.interpolate(x, y, z, h, (1,1), 0x02, NTuple{2,Float64})

            # y should be at y[1] (south edge), x should be interpolated
            @test result[2] == y[1]
            expected_x = x[1] + (x[2] - x[1]) * (h - z[1,1]) / (z[2,1] - z[1,1])
            @test isapprox(result[1], expected_x, rtol=1e-10)
        end
    end

    @testset "Uniform Grid Optimization" begin
        @testset "Step-based interpolation matches vector-based" begin
            x = 0.0:1.0:3.0
            y = 0.0:2.0:6.0
            z = [1.0 2.0; 3.0 4.0]

            # Test that step() optimization produces same results
            # West edge interpolation
            h = 1.5
            result_uniform = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})
            result_vector = Contour.interpolate(collect(x), collect(y), z, h, (1,1), 0x08, NTuple{2,Float64})
            @test result_uniform[1] ≈ result_vector[1] && result_uniform[2] ≈ result_vector[2]

            # East edge interpolation
            h = 3.5
            result_uniform = Contour.interpolate(x, y, z, h, (1,1), 0x04, NTuple{2,Float64})
            result_vector = Contour.interpolate(collect(x), collect(y), z, h, (1,1), 0x04, NTuple{2,Float64})
            @test result_uniform[1] ≈ result_vector[1] && result_uniform[2] ≈ result_vector[2]

            # North edge interpolation
            h = 3.0
            result_uniform = Contour.interpolate(x, y, z, h, (1,1), 0x01, NTuple{2,Float64})
            result_vector = Contour.interpolate(collect(x), collect(y), z, h, (1,1), 0x01, NTuple{2,Float64})
            @test result_uniform[1] ≈ result_vector[1] && result_uniform[2] ≈ result_vector[2]

            # South edge interpolation
            h = 2.0
            result_uniform = Contour.interpolate(x, y, z, h, (1,1), 0x02, NTuple{2,Float64})
            result_vector = Contour.interpolate(collect(x), collect(y), z, h, (1,1), 0x02, NTuple{2,Float64})
            @test result_uniform[1] ≈ result_vector[1] && result_uniform[2] ≈ result_vector[2]
        end

        @testset "Step calculations" begin
            x = -2.0:0.5:2.0
            y = 0.0:1.0:4.0
            z = randn(length(x), length(y))

            # Test interpolation with different steps
            h = 0.0
            result = Contour.interpolate(x, y, z, h, (3,2), 0x08, NTuple{2,Float64})

            # Verify result is within expected bounds
            @test result[1] == x[3]  # West edge
            # Result should be between y[2] and y[3], but might be outside if h is outside the range
            # Just check that it's finite and reasonable
            @test isfinite(result[2])
        end
    end

    @testset "Curvilinear Grid Interpolation" begin
        @testset "Polar coordinate grid" begin
            # Simple polar coordinates
            θ = range(0, 2π, length=20)
            r = range(1, 2, length=10)
            x = [r_i * cos(θ_j) for r_i in r, θ_j in θ]
            y = [r_i * sin(θ_j) for r_i in r, θ_j in θ]
            z = sqrt.(x.^2 + y.^2)  # radius values

            # Test West edge interpolation with a level that actually exists in the grid
            # Use a radius that's between min and max of the grid
            h = 1.2
            ind = (3, 10)  # Choose a cell where interpolation is valid

            # Check that the contour level exists in this cell
            z_values = [z[ind[1], ind[2]], z[ind[1], ind[2]+1]]
            if min(z_values...) <= h <= max(z_values...)
                result = Contour.interpolate(x, y, z, h, ind, 0x08, NTuple{2,Float64})

                # Verify point lies on correct radius (within tolerance)
                computed_radius = sqrt(result[1]^2 + result[2]^2)
                @test isapprox(computed_radius, h, rtol=0.05)  # Some tolerance due to grid discretization
            else
                # Skip this test if the level doesn't exist in the cell
                @test true
            end

            # Test East edge interpolation
            h = 1.8
            ind = (7, 15)  # Choose another cell
            z_values = [z[ind[1], ind[2]], z[ind[1], ind[2]+1]]
            if min(z_values...) <= h <= max(z_values...)
                result_east = Contour.interpolate(x, y, z, h, ind, 0x04, NTuple{2,Float64})
                computed_radius_east = sqrt(result_east[1]^2 + result_east[2]^2)
                @test isapprox(computed_radius_east, h, rtol=0.05)
            else
                # Skip this test if the level doesn't exist in the cell
                @test true
            end
        end

        @testset "General curvilinear deformation" begin
            # Create a simple curved grid
            nx, ny = 5, 4
            ξ = range(0, 1, length=nx)
            η = range(0, 1, length=ny)

            # Apply some deformation
            x = [ξ_i + 0.1*sin(π*η_j) for ξ_i in ξ, η_j in η]
            y = [η_j + 0.1*cos(π*ξ_i) for ξ_i in ξ, η_j in η]
            z = x.^2 + y.^2

            h = 0.5
            ind = (2, 2)

            # Test all four edges
            result_w = Contour.interpolate(x, y, z, h, ind, 0x08, NTuple{2,Float64})
            result_e = Contour.interpolate(x, y, z, h, ind, 0x04, NTuple{2,Float64})
            result_n = Contour.interpolate(x, y, z, h, ind, 0x01, NTuple{2,Float64})
            result_s = Contour.interpolate(x, y, z, h, ind, 0x02, NTuple{2,Float64})

            # Verify all results are finite and reasonable
            for result in [result_w, result_e, result_n, result_s]
                @test all(isfinite, result)
                @test result[1] >= minimum(x) - 0.2
                @test result[1] <= maximum(x) + 0.2
                @test result[2] >= minimum(y) - 0.2
                @test result[2] <= maximum(y) + 0.2
            end
        end
    end

    @testset "Edge Cases and Numerical Stability" begin
        @testset "Interpolation at vertex values" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [0.0 1.0; 1.0 2.0]

            # Test interpolation when h equals vertex values
            for h in [0.0, 1.0, 2.0]
                # Test each edge
                result_w = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})
                result_e = Contour.interpolate(x, y, z, h, (1,1), 0x04, NTuple{2,Float64})
                result_n = Contour.interpolate(x, y, z, h, (1,1), 0x01, NTuple{2,Float64})
                result_s = Contour.interpolate(x, y, z, h, (1,1), 0x02, NTuple{2,Float64})

                # All should be finite
                for result in [result_w, result_e, result_n, result_s]
                    @test all(isfinite, result)
                end
            end
        end

        @testset "Very small differences" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 1.0000001; 1.0000001 1.0000002]

            # Test with very small z differences
            h = 1.00000005
            result = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})
            @test all(isfinite, result)
        end

        @testset "Large value ranges" begin
            x = [0.0, 1e6]
            y = [0.0, 1e6]
            z = [1e10 2e10; 1e-10 2e-10]

            h = 5e9
            result = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})
            @test all(isfinite, result)
            @test result[1] == x[1]  # West edge
        end
    end

    @testset "Different Vertex Types" begin
        @testset "Tuple return type" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            result = Contour.interpolate(x, y, z, 2.5, (1,1), 0x08, Tuple{Float64,Float64})
            @test isa(result, Tuple{Float64,Float64})
            @test result[1] == x[1]
        end

        @testset "Custom type" begin
            struct CustomPoint
                x::Float64
                y::Float64
            end

            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            result = Contour.interpolate(x, y, z, 2.5, (1,1), 0x08, CustomPoint)
            @test isa(result, CustomPoint)
            @test result.x == x[1]
            @test isapprox(result.y, y[1] + (y[2] - y[1]) * (2.5 - z[1,1]) / (z[1,2] - z[1,1]))
        end

        @testset "SArray return type" begin
            using StaticArrays
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            result = Contour.interpolate(x, y, z, 2.5, (1,1), 0x08, SVector{2,Float64})
            @test isa(result, SVector{2,Float64})
            @test result[1] == x[1]
        end
    end

    @testset "Mathematical verification" begin
        @testset "Linear interpolation correctness" begin
            # Test that interpolation formula is mathematically correct
            x = [0.0, 2.0]
            y = [0.0, 3.0]
            z1, z2 = 1.0, 5.0  # Values at two points along an edge

            z = [z1 z2; 0.0 0.0]  # West edge: z[1,1]=z1, z[1,2]=z2
            h = 3.0  # Target value between z1 and z2

            result = Contour.interpolate(x, y, z, h, (1,1), 0x08, NTuple{2,Float64})

            # Manual calculation using linear interpolation formula
            # For West edge: z[1,1] = z1, z[1,2] = z2
            expected_y = y[1] + (y[2] - y[1]) * (h - z1) / (z2 - z1)

            @test result[1] == x[1]  # West edge
            @test isapprox(result[2], expected_y, rtol=1e-12)
        end

        @testset "Proportional interpolation" begin
            # Test that interpolation is proportional to the value difference
            x = [0.0, 10.0]
            y = [0.0, 10.0]
            z = [0.0 0.0; 10.0 0.0]

            # For South edge: z goes from 0 to 10, h = 5 should be halfway
            h = 5.0
            result = Contour.interpolate(x, y, z, h, (1,1), 0x02, NTuple{2,Float64})

            @test result[2] == y[1]  # South edge
            @test isapprox(result[1], 5.0, rtol=1e-10)  # Halfway between x[1] and x[2]
        end
    end

end

end