"GSpace."
struct GSpace{G}
    "Group."
    group::G

    "Dimension of base space."
    dim::Int
end

basespace_action((; group)::GSpace) = Irrep(group, 1)

function build_kernel_basis(
    gspace::GSpace{CyclicGroup},
    ρin,
    ρout,
    σ,
    rings,
    maximum_frequency,
)
    (; group) = gspace
    (; n) = group
    sg_id = (nothing, n)
    radial_profile = GaussianRadialProfile(radii, sigma)
    basis = SteerableKernelBasis(
        CircularShellsBasis(maximum_frequency, radial_profile),
        in_repr,
        out_repr,
        RestrictedWignerEckartBasis,
        sg_id,
    )
end

GaussianRadialProfile(radii, σ) = r -> map((r0, σ) -> exp(-(r - r0)^2 / 2σ^2), radii, σ)
function CircularShellsBasis(maximum_frequency, radial_profile) end

"Field type."
struct FieldType{G,R}
    "GSpace."
    gspace::G

    "Representation."
    representation::R
end

"Fiber field."
struct FiberField{F,A}
    "Field type."
    fieldtype::F

    "Field content."
    x::A

    function FiberField(fieldtype, x)
        @assert size(x, 1) == size(x, 2) "Domain must be square"
        new{typeof(fieldtype),typeof(x)}(fieldtype, x)
    end
end

function Base.:*(g::Element, f::FiberField)
    (; fieldtype, x) = f
    (; gspace, representation) = fieldtype
    (; group, dim) = gspace
    n = size(x, 1)
    A = basespace_action(gspace)(inv(g))
    @assert dim == 2 "Only 2D for now."
    R = representation(g)

    # Transform fibers
    @tensor y[i, j, a] := R[a, b] * x[i, j, b]

    ispermutation =
        all(≈(1), sum(abs, A; dims = 1)) &&
        all(≈(1), sum(abs, A; dims = 2)) &&
        A ≈ round.(Int, A)

    if ispermutation
        A = round.(Int, A)

        # Transform base space
        Rx = findfirst(a -> abs(a) == 1, A[1, :])
        Ry = findfirst(a -> abs(a) == 1, A[2, :])
        sx = sign(A[1, Rx])
        sy = sign(A[2, Ry])
        i = sx == 1 ? (1:n) : (n:-1:1)
        i = Rx == 1 ? i : i'
        j = sy == 1 ? (1:n) : (n:-1:1)
        j = Ry == 1 ? j : j'
        I = CartesianIndex.(i, j)
        y = y[I, :]
    else
        # Get coordinates of the new grid (Cartesian)
        # Transform coordinates with inverse element
        # Interpolate old field to these coordinates
        grid = get_grid_coords(dim, n, 1)
        Agrid = A * grid
        bounds = extrema(grid)
        scalargrid = range(bounds..., n)
        i = @. (bounds[1] < Agrid[1, :] < bounds[2]) && (bounds[1] < Agrid[2, :] < bounds[2])
        y = zero(x)
        y = stack(1:size(x, 3)) do a
            field = x[:, :, a]
            p = linear_interpolation((scalargrid, scalargrid), field)
            ya = zeros(n^2)
            @. ya[i] = p(Agrid[1, i], Agrid[2, i])
            reshape(ya, n, n)
        end
    end

    FiberField(fieldtype, y)
end

# Tensors/einsums:
# - https://github.com/mcabbott/Tullio.jl
# - https://github.com/under-Peter/OMEinsum.jl
# - https://github.com/ITensor/ITensors.jl
# - https://github.com/Jutho/TensorKit.jl
# - https://github.com/Jutho/TensorOperations.jl
