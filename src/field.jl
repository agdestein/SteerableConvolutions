"GSpace."
struct GSpace{G}
    "Group."
    group::G

    "Dimension of base space."
    dim::Int
end

basespace_action(gspace::GSpace) = irrep(gspace.group, 1)

"Field type."
struct FieldType{G,R}
    "GSpace."
    gspace::G

    "List of representations."
    representations::R
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
    (; gspace, representations) = fieldtype
    (; group, dim) = gspace
    n = size(x, 1)
    A = basespace_action(gspace)(inv(g))
    @assert dim == 2 "Only 2D for now."
    ρ = directsum(representations...)
    R = ρ(g)

    # Transform fibers
    @tensor y[i, j, a] := R[a, b] * x[i, j, b]

    ispermutation =
        all(≈(1), sum(abs, A; dims = 1)) &&
        all(≈(1), sum(abs, A; dims = 2)) &&
        A ≈ round.(Int, A)

    @assert ispermutation "Only pure permutation works for now."
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

    FiberField(fieldtype, y)
end

# Tensors/einsums:
# - https://github.com/mcabbott/Tullio.jl
# - https://github.com/under-Peter/OMEinsum.jl
# - https://github.com/ITensor/ITensors.jl
# - https://github.com/Jutho/TensorKit.jl
# - https://github.com/Jutho/TensorOperations.jl
