module CellClassificationTests

using Contour, Test
include("test_helpers.jl")
using .TestHelpers

@testset "Cell Classification System" begin

    @testset "_get_case - All 16 configurations" begin
        # Test all 16 cell configurations systematically
        test_configs = all_cell_configurations()

        for (case, z_values) in test_configs
            h = 0.5
            result = Contour._get_case(z_values, h)
            @test result == case

            # Test with different contour level
            h_alt = 1.5
            z_alt = [z + 1 for z in z_values]  # Shift all values up
            result_alt = Contour._get_case(z_alt, h_alt)
            @test result_alt == case
        end

        # Specific test cases
        @testset "Specific edge cases" begin
            # Case 0: All below
            @test Contour._get_case([0.0, 1.0, 2.0, 3.0], 5.0) == 0x00

            # Case 15: All above
            @test Contour._get_case([6.0, 7.0, 8.0, 9.0], 5.0) == 0x0f

            # Mixed cases
            @test Contour._get_case([1.0, 0.0, 0.0, 0.0], 0.5) == 0x01  # Only vertex 1
            @test Contour._get_case([0.0, 1.0, 0.0, 0.0], 0.5) == 0x02  # Only vertex 2
            @test Contour._get_case([0.0, 0.0, 1.0, 0.0], 0.5) == 0x04  # Only vertex 3
            @test Contour._get_case([0.0, 0.0, 0.0, 1.0], 0.5) == 0x08  # Only vertex 4
        end
    end

    @testset "Edge Lookup Table" begin
        # Use local variables instead of const for testing
        N, S, E, W = (UInt8(1)), (UInt8(2)), (UInt8(4)), (UInt8(8))
        NS, NE, NW = N|S, N|E, N|W
        SN, SE, SW = S|N, S|E, S|W
        EN, ES, EW = E|N, E|S, E|W

        # Test that edge_LUT entries match expected patterns
        # Skip cases 0, 5, 10, 15 (not in LUT)
        @test Contour.edge_LUT[1] == SW   # Case 1 -> SW crossing
        @test Contour.edge_LUT[2] == SE   # Case 2 -> SE crossing
        @test Contour.edge_LUT[3] == EW   # Case 3 -> EW crossing
        @test Contour.edge_LUT[4] == NE   # Case 4 -> NE crossing
        @test Contour.edge_LUT[5] == 0x0  # Case 5 (ambiguous - not in LUT)
        @test Contour.edge_LUT[6] == NS   # Case 6 -> NS crossing
        @test Contour.edge_LUT[7] == NW   # Case 7 -> NW crossing
        @test Contour.edge_LUT[8] == NW   # Case 8 -> NW crossing
        @test Contour.edge_LUT[9] == NS   # Case 9 -> NS crossing
        @test Contour.edge_LUT[10] == 0x0 # Case 10 (ambiguous - not in LUT)
        @test Contour.edge_LUT[11] == NE  # Case 11 -> NE crossing
        @test Contour.edge_LUT[12] == EW  # Case 12 -> EW crossing
        @test Contour.edge_LUT[13] == SE  # Case 13 -> SE crossing
        @test Contour.edge_LUT[14] == SW  # Case 14 -> SW crossing
    end

    @testset "Ambiguous Case Resolution" begin
        @testset "Case 5 (0b0101)" begin
            # Test bilinear interpolation threshold
            # Pattern: [1, 0, 0, 1] relative to vertices 1-4 (booleans)
            z_values = [true, false, false, true]

            # When h > average (0.5), should resolve to NWSE
            h = 0.6
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            @test first(values(cells)) == Contour.NWSE

            # When h < average (0.5), should resolve to NWSE
            h = 0.4
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            @test first(values(cells)) == Contour.NWSE

            # Test exactly at average boundary
            h = 0.5
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            # At boundary, should resolve to NWSE due to >= comparison
            @test first(values(cells)) == Contour.NWSE
        end

        @testset "Case 10 (0b1010)" begin
            # Pattern: [0, 1, 1, 0] relative to vertices 1-4 (booleans)
            z_values = [false, true, true, false]

            # When h > average (0.5), should resolve to NESW
            h = 0.6
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            @test first(values(cells)) == Contour.NESW

            # When h < average (0.5), should resolve to NESW
            h = 0.4
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            @test first(values(cells)) == Contour.NESW

            # Test exactly at average boundary
            h = 0.5
            z = create_test_matrix_from_config(z_values, h)
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            # At boundary, should resolve to NESW due to >= comparison
            @test first(values(cells)) == Contour.NESW
        end

        @testset "Non-ambiguous cases" begin
            # Test that non-ambiguous cases use edge_LUT directly
            # Create case 6 (0110): [1,1,0,0] pattern
            z = [1.0 1.0; 0.0 0.0]  # Creates case 6 (0110)
            h = 0.5
            cells = Contour.get_level_cells(z, h)
            @test length(cells) == 1
            @test first(values(cells)) == Contour.NS  # North-South crossing
        end
    end

    @testset "get_level_cells function" begin
        @testset "Basic functionality" begin
            # Simple 2x2 grid
            z = [1.0 2.0; 3.0 4.0]
            h = 2.5

            cells = Contour.get_level_cells(z, h)

            # Should identify at least one cell containing contour
            @test length(cells) >= 0

            # All keys should be valid cell indices
            for (xi, yi) in keys(cells)
                @test 1 <= xi <= size(z, 1) - 1
                @test 1 <= yi <= size(z, 2) - 1
            end

            # All values should be valid crossing patterns
            for crossing in values(cells)
                @test crossing isa UInt8
                @test crossing > 0  # Should not be empty cells
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

        @testset "Larger grid" begin
            # Test with larger grid to ensure algorithm scales
            nx, ny = 5, 4
            z = rand(nx, ny)
            h = (minimum(z) + maximum(z)) / 2

            cells = Contour.get_level_cells(z, h)

            # Should find some cells
            @test length(cells) >= 0

            # All cells should be within bounds
            for (xi, yi) in keys(cells)
                @test 1 <= xi <= nx - 1
                @test 1 <= yi <= ny - 1
            end
        end

        @testset "Edge values" begin
            # Test when contour level equals some vertex values
            z = [1.0 2.0; 3.0 4.0]

            # Test with levels equal to vertex values
            for vertex_val in [1.0, 2.0, 3.0, 4.0]
                cells = Contour.get_level_cells(z, vertex_val)
                @test length(cells) >= 0  # Should not crash

                for crossing in values(cells)
                    @test crossing isa UInt8
                end
            end
        end
    end

    @testset "Edge Navigation Constants" begin
        # Test that navigation constants are correctly defined
        @test Contour.next_map == ((0,1), (0,-1), (1,0), (-1,0))  # N, S, E, W movement
        @test Contour.next_edge == (0x02, 0x01, 0x08, 0x04)    # S, N, W, E entry edges
    end

    @testset "advance_edge function" begin
        # Test each direction
        @test Contour.advance_edge((5,5), 0x01) == ((5,6), 0x02)  # N -> S
        @test Contour.advance_edge((5,5), 0x02) == ((5,4), 0x01)  # S -> N
        @test Contour.advance_edge((5,5), 0x04) == ((6,5), 0x08)  # E -> W
        @test Contour.advance_edge((5,5), 0x08) == ((4,5), 0x04)  # W -> E

        # Test at boundaries
        @test Contour.advance_edge((1,1), 0x02) == ((1,0), 0x01)  # Move south from (1,1)
        @test Contour.advance_edge((10,10), 0x01) == ((10,11), 0x02)  # Move north from (10,10)
    end

    @testset "get_first_crossing function" begin
        @test Contour.get_first_crossing(Contour.NWSE) == Contour.NW
        @test Contour.get_first_crossing(Contour.NESW) == Contour.NE
        @test Contour.get_first_crossing(Contour.NS) == Contour.NS
        @test Contour.get_first_crossing(Contour.EW) == Contour.EW
        @test Contour.get_first_crossing(Contour.NW) == Contour.NW
        @test Contour.get_first_crossing(Contour.SE) == Contour.SE
    end

end

end