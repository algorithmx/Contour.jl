using Contour, Test

@test length(detect_ambiguities(Base, Contour)) == 0

# Include test helpers first
include("test_helpers.jl")

# Include existing tests
include("verify_vertices.jl")
include("interface.jl")

# Include new comprehensive unit tests
include("test_cell_classification.jl")
include("test_interpolation.jl")
include("test_tracing.jl")
include("test_edge_cases.jl")
include("test_types.jl")

#issue 59
@inferred collect(())

using Aqua
# Aqua tests
# Intervals brings a bunch of ambiquities unfortunately
Aqua.test_all(Contour)


@static if Base.VERSION >= v"1.7"

    @info "Running JET..."

    using JET
    display(JET.report_package(Contour))
end