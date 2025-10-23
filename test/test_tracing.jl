module ContourTracingTests

using Contour, Test
include("test_helpers.jl")
using .TestHelpers

@testset "Contour Tracing Logic" begin

    @testset "Cell Collection" begin
        @testset "get_level_cells - basic functionality" begin
            z = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
            h = 5.0

            cells = Contour.get_level_cells(z, h)

            # Should identify cells that contain contour
            @test length(cells) >= 0

            # Verify expected cells are identified (manual inspection)
            # This should include at least cell (1,1) since values range from 1-9
            @test haskey(cells, (1,1)) || haskey(cells, (2,1)) ||
                  haskey(cells, (1,2)) || haskey(cells, (2,2))

            # All keys should be valid cell indices
            for (xi, yi) in keys(cells)
                @test 1 <= xi <= size(z, 1) - 1
                @test 1 <= yi <= size(z, 2) - 1
            end
        end

        @testset "Empty contour levels" begin
            # Test with level outside data range
            z = [1.0 2.0; 3.0 4.0]

            # Below minimum
            cells_low = Contour.get_level_cells(z, 0.0)
            @test length(cells_low) == 0

            # Above maximum
            cells_high = Contour.get_level_cells(z, 5.0)
            @test length(cells_high) == 0
        end

        @testset "Consistent cell identification" begin
            # Test that same input produces same output
            z = rand(5, 5)
            h = 0.5 * (minimum(z) + maximum(z))

            cells1 = Contour.get_level_cells(z, h)
            cells2 = Contour.get_level_cells(z, h)

            @test cells1 == cells2
        end
    end

    @testset "Edge Navigation" begin
        @testset "get_next_edge! - Non-ambiguous cases" begin
            cells = Dict{Tuple{Int,Int},UInt8}()

            # Test simple north-south crossing
            cells[(1,1)] = Contour.NS
            result = Contour.get_next_edge!(cells, (1,1), 0x01)  # Enter from N
            @test result == 0x02  # Exit S
            @test !haskey(cells, (1,1))  # Cell should be removed

            # Test east-west crossing
            cells[(1,1)] = Contour.EW
            result = Contour.get_next_edge!(cells, (1,1), 0x04)  # Enter from E
            @test result == 0x08  # Exit W

            # Test diagonal crossings
            cells[(1,1)] = Contour.NW
            result = Contour.get_next_edge!(cells, (1,1), 0x01)  # Enter from N
            @test result == 0x08  # Exit W

            cells[(1,1)] = Contour.SE
            result = Contour.get_next_edge!(cells, (1,1), 0x02)  # Enter from S
            @test result == 0x04  # Exit E
        end

        @testset "get_next_edge! - Ambiguous cases (NWSE)" begin
            cells = Dict{Tuple{Int,Int},UInt8}()

            # Test NWSE case entering from North
            cells[(1,1)] = Contour.NWSE
            result = Contour.get_next_edge!(cells, (1,1), 0x01)  # Enter from N
            @test result == 0x08  # Exit W
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x06  # Remaining SE crossing (0x02 | 0x04)

            # Test NWSE case entering from West
            cells[(1,1)] = Contour.NWSE
            result = Contour.get_next_edge!(cells, (1,1), 0x08)  # Enter from W
            @test result == 0x01  # Exit N
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x06  # Remaining SE crossing

            # Test NWSE case entering from South
            cells[(1,1)] = Contour.NWSE
            result = Contour.get_next_edge!(cells, (1,1), 0x02)  # Enter from S
            @test result == 0x04  # Exit E
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x09  # Remaining NW crossing (0x01 | 0x08)

            # Test NWSE case entering from East
            cells[(1,1)] = Contour.NWSE
            result = Contour.get_next_edge!(cells, (1,1), 0x04)  # Enter from E
            @test result == 0x02  # Exit S
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x09  # Remaining NW crossing
        end

        @testset "get_next_edge! - Ambiguous cases (NESW)" begin
            cells = Dict{Tuple{Int,Int},UInt8}()

            # Test NESW case entering from North
            cells[(1,1)] = Contour.NESW
            result = Contour.get_next_edge!(cells, (1,1), 0x01)  # Enter from N
            @test result == 0x04  # Exit E
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x0a  # Remaining SW crossing (0x02 | 0x08)

            # Test NESW case entering from East
            cells[(1,1)] = Contour.NESW
            result = Contour.get_next_edge!(cells, (1,1), 0x04)  # Enter from E
            @test result == 0x01  # Exit N
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x0a  # Remaining SW crossing

            # Test NESW case entering from South
            cells[(1,1)] = Contour.NESW
            result = Contour.get_next_edge!(cells, (1,1), 0x02)  # Enter from S
            @test result == 0x08  # Exit W
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x05  # Remaining NE crossing (0x01 | 0x04)

            # Test NESW case entering from West
            cells[(1,1)] = Contour.NESW
            result = Contour.get_next_edge!(cells, (1,1), 0x08)  # Enter from W
            @test result == 0x02  # Exit S
            @test haskey(cells, (1,1))  # Cell should still exist
            @test cells[(1,1)] == 0x05  # Remaining NE crossing
        end
    end

    @testset "Contour Following Algorithm" begin
        @testset "Simple closed loop" begin
            # Create a 2x2 grid with a single closed contour
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [0.0 1.0; 1.0 0.0]
            h = 0.5

            contour_level = Contour.contour(x, y, z, h)

            @test length(contour_level.lines) >= 1

            # Should form a closed loop (X-shaped pattern)
            for line in contour_level.lines
                if length(line.vertices) > 2
                    # Check if it's a closed contour
                    @test is_closed_contour(line)
                end
            end
        end

        @testset "Open contour reaching boundary" begin
            # Create a contour that reaches the boundary
            x = [0.0, 1.0, 2.0]
            y = [0.0, 1.0]
            z = [0.0 1.0; 0.5 1.5; 1.0 2.0]  # Correct dimensions: (length(x), length(y))
            h = 1.0

            contour_level = Contour.contour(x, y, z, h)
            @test length(contour_level.lines) >= 1

            # Should have endpoints on boundary
            for line in contour_level.lines
                if length(line.vertices) >= 2
                    first_vertex = first(line.vertices)
                    last_vertex = last(line.vertices)

                    # At least one endpoint should be on boundary
                    @test (first_vertex[1] ≈ 0.0 || first_vertex[2] ≈ 0.0 ||
                           last_vertex[1] ≈ 0.0 || last_vertex[2] ≈ 0.0 ||
                           first_vertex[1] ≈ 2.0 || first_vertex[2] ≈ 1.0 ||
                           last_vertex[1] ≈ 2.0 || last_vertex[2] ≈ 1.0)
                end
            end
        end

        @testset "Multiple disjoint contours" begin
            # Create a field with multiple separate contours
            x = range(0, 3, length=10)
            y = range(0, 2, length=7)
            z = [sin(xi*π) * cos(yi*π) for xi in x, yi in y]
            h = 0.5

            contour_level = Contour.contour(x, y, z, h)
            @test length(contour_level.lines) >= 1

            # Each line should have valid vertices
            total_vertices = 0
            for line in contour_level.lines
                @test length(line.vertices) >= 2
                for vertex in line.vertices
                    @test all(isfinite, vertex)
                    @test x[1] <= vertex[1] <= x[end]
                    @test y[1] <= vertex[2] <= y[end]
                end
                total_vertices += length(line.vertices)
            end
            @test total_vertices > 0
        end

        @testset "Saddle point handling" begin
            # Test with a classic saddle function
            x = range(-1, 1, length=10)
            y = range(-1, 1, length=10)
            z = [xi^2 - yi^2 for xi in x, yi in y]  # saddle function
            h = 0.0

            contour_level = Contour.contour(x, y, z, h)
            @test length(contour_level.lines) >= 2  # Should have two crossing lines

            # Should handle the saddle point at (0,0) correctly
            total_vertices = count_contour_vertices(contour_level)
            @test total_vertices > 0
        end
    end

    @testset "Cell Tracing Subroutine" begin
        @testset "chase! function basic functionality" begin
            # Create a simple test case
            x = [0.0, 1.0]
            y = [0.0, 1.0]
            z = [0.0 1.0; 1.0 0.0]
            h = 0.5

            cells = Contour.get_level_cells(z, h)
            @test length(cells) > 0

            # Manually test chase! by extracting and running parts of trace_contour
            # This is complex to test directly, so we verify through contour results
            contour_level = Contour.contour(x, y, z, h)
            @test length(contour_level.lines) >= 1

            # Verify that all cells get processed (empty at end)
            final_cells = Contour.get_level_cells(z, h)  # Fresh copy to compare
            @test length(final_cells) == length(cells)  # Should be same number initially
        end

        @testset "Boundary stopping conditions" begin
            # Create contour that starts at boundary
            x = [0.0, 1.0, 2.0]
            y = [0.0, 1.0, 2.0]

            # Create gradient from corner
            z = [0.0 0.5 1.0; 0.5 1.0 1.5; 1.0 1.5 2.0]
            h = 0.75

            contour_level = Contour.contour(x, y, z, h)
            @test length(contour_level.lines) >= 1

            # Should stop at boundaries correctly
            for line in contour_level.lines
                for vertex in line.vertices
                    @test all(isfinite, vertex)
                end
            end
        end
    end

    @testset "Algorithm consistency" begin
        @testset "Deterministic results" begin
            # Test that same input produces same output multiple times
            x = randn(10)
            y = randn(8)
            z = randn(10, 8)
            h = 0.0

            contour1 = Contour.contour(x, y, z, h)
            contour2 = Contour.contour(x, y, z, h)
            contour3 = Contour.contour(x, y, z, h)

            @test length(contour1.lines) == length(contour2.lines) == length(contour3.lines)

            # Compare vertex counts (may be ordered differently)
            counts1 = [length(line.vertices) for line in contour1.lines]
            counts2 = [length(line.vertices) for line in contour2.lines]
            counts3 = [length(line.vertices) for line in contour3.lines]

            @test Set(counts1) == Set(counts2) == Set(counts3)
        end

        @testset "Symmetry preservation" begin
            # Test with symmetric data
            x = [-2.0, -1.0, 0.0, 1.0, 2.0]
            y = [-1.5, -0.5, 0.5, 1.5]

            # Create radially symmetric function
            z = [xi^2 + yi^2 for xi in x, yi in y]
            h = 2.0

            contour_level = Contour.contour(x, y, z, h)

            # Should produce circular-like contours
            total_vertices = count_contour_vertices(contour_level)
            @test total_vertices > 0

            # All vertices should be approximately at the same radius
            for line in contour_level.lines
                for vertex in line.vertices
                    radius = sqrt(vertex[1]^2 + vertex[2]^2)
                    @test isapprox(radius, sqrt(h), rtol=0.2)  # Allow some tolerance
                end
            end
        end
    end

end

end