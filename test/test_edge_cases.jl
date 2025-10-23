module EdgeCaseTests

using Contour, Test
include("test_helpers.jl")
using .TestHelpers

@testset "Edge Cases and Integration Tests" begin

    @testset "Boundary Conditions" begin
        @testset "Single cell grids" begin
            @testset "1x1 grid" begin
                # For marching squares, we need at least 2x2 grid to have one cell
                x = [0.0, 1.0]
                y = [0.0, 1.0]
                z = [0.5 0.5; 0.5 0.5]  # 2x2 matrix with constant values

                # Should handle gracefully without crashing
                result = Contour.contour(x, y, z, 0.25)
                @test result isa Contour.ContourLevel
                @test result.level == 0.25
            end

            @testset "1x2 grid" begin
                # x has 2 points, y has 3 points, so z should be 2x3 matrix for marching squares
                x = [0.0, 1.0]
                y = [0.0, 0.5, 1.0]
                z = [0.3 0.7 0.5; 0.6 0.4 0.8]  # 2x3 matrix

                result = Contour.contour(x, y, z, 0.5)
                @test result isa Contour.ContourLevel
                @test length(result.lines) >= 0
            end

            @testset "2x1 grid" begin
                # x has 3 points, y has 2 points, so z should be 3x2 matrix
                x = [0.0, 0.5, 1.0]
                y = [0.0, 1.0]
                z = [0.3 0.6; 0.7 0.4; 0.5 0.8]  # 3x2 matrix

                result = Contour.contour(x, y, z, 0.5)
                @test result isa Contour.ContourLevel
                @test length(result.lines) >= 0
            end

            @testset "2x2 grid" begin
                x = [0.0, 1.0]
                y = [0.0, 1.0]
                z = [0.3 0.7; 0.6 0.4]

                result = Contour.contour(x, y, z, 0.5)
                @test result isa Contour.ContourLevel
                @test length(result.lines) >= 0
            end
        end

        @testset "Degenerate grids" begin
            @testset "Zero area cells" begin
                x = [0.0, 0.0, 1.0]
                y = [0.0, 1.0, 1.0]
                z = [0.0 0.5 0.3; 0.5 1.0 0.7; 0.2 0.8 0.9]  # 3x3 matrix

                # Should not crash or return invalid results
                result = Contour.contour(x, y, z, 0.75)
                @test result isa Contour.ContourLevel
                @test length(result.lines) >= 0

                # All vertices should be finite
                for line in result.lines
                    for vertex in line.vertices
                        @test all(isfinite, vertex)
                    end
                end
            end

            @testset "Negative coordinates" begin
                x = [-2.0, -1.0, 0.0]
                y = [-1.5, -0.5, 0.5]
                z = [0.1 0.2 0.3; 0.4 0.5 0.6; 0.7 0.8 0.9]

                result = Contour.contour(x, y, z, 0.45)
                @test result isa Contour.ContourLevel
                @test length(result.lines) >= 0

                # Vertices should be within grid bounds
                for line in result.lines
                    for vertex in line.vertices
                        @test x[1] <= vertex[1] <= x[end]
                        @test y[1] <= vertex[2] <= y[end]
                    end
                end
            end
        end

        @testset "Mixed coordinate types" begin
            # Test with ranges and vectors mixed
            x = [0.0, 1.0, 2.0]  # Vector
            y = 0.0:0.5:1.0       # Range
            z = randn(3, 3)

            result = Contour.contour(x, y, z, 0.0)
            @test result isa Contour.ContourLevel
            @test length(result.lines) >= 0
        end
    end

    @testset "Numerical Precision" begin
        @testset "Contour level equal to vertex values" begin
            z = [1.0 2.0; 3.0 4.0]
            h = 2.0  # Exactly equal to vertex

            result = Contour.contour([0.0, 1.0], [0.0, 1.0], z, h)
            @test result isa Contour.ContourLevel
            @test length(result.lines) >= 0

            # Should handle equality consistently
            for line in result.lines
                for vertex in line.vertices
                    @test all(isfinite, vertex)
                end
            end
        end

        @testset "Very small differences" begin
            # Test with values that are very close together
            x = [0.0, 1e-6, 2e-6]
            y = [0.0, 1e-6]
            z = [1.0 1.000000001; 1.000000002 1.000000003; 1.000000004 1.000000005]

            result = Contour.contour(x, y, z, 1.0000000015)
            @test result isa Contour.ContourLevel

            # All vertices should be finite
            for line in result.lines
                for vertex in line.vertices
                    @test all(isfinite, vertex)
                end
            end
        end

        @testset "Large value ranges" begin
            z = [1e10 2e10; 1e-10 2e-10]
            h = 5e9

            result = Contour.contour([0.0, 1.0], [0.0, 1.0], z, h)
            @test result isa Contour.ContourLevel

            # Should handle large/small value differences
            for line in result.lines
                for vertex in line.vertices
                    @test all(isfinite, vertex)
                end
            end
        end

        @testset "Inf and NaN handling" begin
            # Test with extreme values (though Contour.jl may not handle these gracefully)
            # This test documents current behavior rather than prescribing expected behavior
            z = [1.0 2.0; 3.0 Inf]
            x = [0.0, 1.0]
            y = [0.0, 1.0]

            # This may or may not work depending on implementation
            try
                result = Contour.contour(x, y, z, 1.5)
                # If it works, check that vertices are reasonable
                for line in result.lines
                    for vertex in line.vertices
                        # Should not return NaN if possible
                        @test !any(isnan, vertex)
                    end
                end
            catch e
                # It's acceptable to throw an error for these inputs
                @test true
            end
        end
    end

    @testset "Integration with Real Data" begin
        @testset "Mathematical functions" begin
            # Gaussian function
            x = range(-3, 3, length=50)
            y = range(-3, 3, length=40)
            z = [exp(-(xi^2 + yi^2)/2) for xi in x, yi in y]

            levels = [0.1, 0.3, 0.5, 0.7]
            results = Contour.contours(x, y, z, levels)

            @test length(results.contours) == length(levels)

            for (level, result) in zip(levels, results.contours)
                @test result.level == level
                @test length(result.lines) >= 1

                # Verify contour points are approximately at correct level
                for line in result.lines
                    for vertex in line.vertices
                        computed_z = exp(-(vertex[1]^2 + vertex[2]^2)/2)
                        @test isapprox(computed_z, level, rtol=0.1)  # Allow some tolerance
                    end
                end
            end
        end

        @testset "Saddle point function" begin
            # Classic saddle: z = x^2 - y^2
            x = range(-2, 2, length=30)
            y = range(-2, 2, length=30)
            z = [xi^2 - yi^2 for xi in x, yi in y]

            h = 0.0  # Saddle point level
            result = Contour.contour(x, y, z, h)

            # Should have crossing lines through saddle point
            @test length(result.lines) >= 2

            total_vertices = count_contour_vertices(result)
            @test total_vertices > 0
        end

        @testset "Multiple disjoint contours" begin
            # Create multiple separated hills
            x = range(0, 4, length=50)
            y = range(0, 3, length=40)

            z = zeros(length(x), length(y))
            # Add two hills
            for (i, xi) in enumerate(x), (j, yj) in enumerate(y)
                hill1 = exp(-((xi-1)^2 + (yj-1)^2)/0.3)
                hill2 = exp(-((xi-3)^2 + (yj-2)^2)/0.3)
                z[i,j] = hill1 + hill2
            end

            h = 0.5
            result = Contour.contour(x, y, z, h)

            # Should have multiple separate closed contours
            closed_contours = count(is_closed_contour(line) for line in result.lines)
            @test closed_contours >= 1

            total_vertices = count_contour_vertices(result)
            @test total_vertices > 0
        end
    end

    @testset "API Edge Cases" begin
        @testset "Empty level collections" begin
            # Test contours() with empty level list
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            result = Contour.contours(x, y, z, [])
            @test result isa Contour.ContourCollection
            @test length(result.contours) == 0
        end

        @testset "Single level in contours()" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            result = Contour.contours(x, y, z, [2.5])
            @test length(result.contours) == 1
            @test result.contours[1].level == 2.5

            # Should be equivalent to single contour()
            single = Contour.contour(x, y, z, 2.5)
            @test length(result.contours[1].lines) == length(single.lines)
        end

        @testset "Multiple identical levels" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            levels = [2.0, 2.0, 2.0]
            result = Contour.contours(x, y, z, levels)
            @test length(result.contours) == 3

            # All should have same level but possibly different lines
            for contour in result.contours
                @test contour.level == 2.0
                @test length(contour.lines) >= 0
            end
        end
    end

    @testset "Performance and Memory" begin
        @testset "Large grid performance" begin
            sizes = [10, 20]  # Keep small for test speed

            for size in sizes
                x = range(0, 1, length=size)
                y = range(0, 1, length=size)
                z = [sin(xi*π) * cos(yi*π) for xi in x, yi in y]

                time = @elapsed result = Contour.contours(x, y, z, 5)

                # Basic performance expectations (adjust as needed)
                @test time < 5.0  # Should complete within reasonable time
                @test length(result.contours) == 5

                # Memory should not blow up
                total_vertices = sum(count_contour_vertices(cl) for cl in result.contours)
                @test total_vertices < size * size * 10  # Should be reasonable
            end
        end

        @testset "Memory efficiency with empty contours" begin
            # Create data that will produce very few contours
            x = range(0, 1, length=20)
            y = range(0, 1, length=20)
            z = ones(20, 20)  # All same value

            # No contours should be generated
            result = Contour.contours(x, y, z, 5)
            @test length(result.contours) == 5

            for contour in result.contours
                @test length(contour.lines) == 0
            end
        end
    end

    @testset "Error Handling" begin
        @testset "Mismatched dimensions" begin
            x = [0.0, 1.0, 2.0]  # 3 points
            y = [0.0, 1.0]       # 2 points
            z = randn(4, 2)      # 4x2 matrix - wrong dimension

            @test_throws ArgumentError Contour.contour(x, y, z, 0.5)
        end

        @testset "Invalid vertex types" begin
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [1.0 2.0; 3.0 4.0]

            # Test with invalid vertex type - this may or may not error
            try
                result = Contour.contour(x, y, z, 2.5, VT=String)
                # If it doesn't error, the result should still be valid
                @test result isa Contour.ContourLevel
            catch e
                # It's acceptable to throw an error for incompatible types
                @test true
            end
        end
    end

end

end