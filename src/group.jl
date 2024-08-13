abstract type AbstractGroup end
abstract type AbstractFiniteGroup <: AbstractGroup end

struct Element{G}
    group::G
    n::Int
end

(group::AbstractGroup)(n) = Element(group, n)
Base.show(io::IO, e::Element) = print(io, e.group, "(", e.n, ")")

Base.one(e::Element) = one(e.group)
Base.:(^)(e::Element, n) = prod(fill(e, n); init = one(e))

abstract type AbstractRepresentation end
struct IrreducibleRepresentation{G,A} <: AbstractRepresentation
    group::G
    id::Int
    mat::A
    type::Char
end
struct Representation{A} <: AbstractRepresentation
    irreps::Vector{Int}
    basis::A
end

((; mat)::IrreducibleRepresentation)(g) = mat(g)
function ((; irreps, basis)::Representation)(g)
    i = directsum(map(i -> irrepmat(g, i), irreps)...)
    basis * i * inv(basis)
end

irreps(r::IrreducibleRepresentation) = [r.id]
basis(r::IrreducibleRepresentation) = one(r.mat(one(r.group)))

irreps(r::Representation) = r.irreps
basis(r::Representation) = r.basis

function directsum(x...)
    n = size.(x, 1)
    y = fill!(similar(x[1], sum(n), sum(n)), 0)
    j = 0
    for (i, x) in enumerate(x)
        k = j .+ (1:n[i])
        j += n[i]
        y[k, k] .= x
    end
    y
end
directsum(r::AbstractRepresentation...) =
    Representation(vcat(irreps.(r)...), directsum(basis.(r)...))
⊕ = directsum

sum_of_squares_constituents(type) =
    if type == 'R'
        1
    elseif type == 'C'
        2
    elseif type == 'Q'
        4
    else
        error("irrep: unknown type $(i.type)")
    end

function regular_representation(group)
    N = order(group)
    e = elements(group)
    representations = map(e) do g
        r = zeros(N, N)
        for (j, h) in enumerate(e)
            i = findfirst(==(h * g), e)
            r[i, j] = 1
        end
        r
    end
    characters = tr.(representations)

    i = irreps(group)
    multiplicities = map(i) do i
        (; type) = irrep(group, i)

        m = map(characters, e) do c, g
            c * tr(irrepmat(inv(g), i))
        end
        m = sum(m) / length(characters) / sum_of_squares_constituents(type)
        mint = round(Int, m)
        @assert mint ≈ m
        mint
    end
    ilist = vcat(fill.(i, multiplicities)...)
    irrs = map(i -> irrep(group, i), ilist)
    P = directsum(irrs...)

    v = zeros(N)
    p = 0
    for (irr, m) in zip(i, multiplicities)
        s = size(irrepmat(one(group), irr), 1)
        @assert s ≥ m
        M = Matrix(I * sqrt(s), s, m)
        v[p+1:p+m*s] .= M[:]
        p += m * s
    end

    basis = zeros(N, N)
    for (i, g) in enumerate(e)
        basis[i, :] .= P(g) * v
    end
    basis ./= sqrt(N)

    Representation(ilist, basis)
end

struct Rotation{T} <: AbstractGroup
    angle::T
    Rotation(angle) = new{typeof(angle)}(mod2pi(angle))
end

Base.one(r::Rotation) = Rotation(zero(r.angle))
Base.:*(f::Rotation, g::Rotation) = Rotation(f.angle + g.angle)
Base.inv(f::Rotation) = Rotation(-f.angle)

struct CyclicGroup <: AbstractFiniteGroup
    N::Int
end

(group::CyclicGroup)(n) = Element(group, mod(n, group.N))
order(group::CyclicGroup) = group.N
Base.one(group::CyclicGroup) = Element(group, 0)
Base.:*(g::Element{CyclicGroup}, h::Element{CyclicGroup}) = g.group(g.n + h.n)
Base.inv(g::Element{CyclicGroup}) = g.group(-g.n)
elements(group::CyclicGroup) = group.(0:group.N-1)
irrepmat(g::Element{CyclicGroup}, i) =
    if 2i > g.group.N
        error("irrepmat: 2i must be less N")
    elseif i == 0
        fill(1.0, 1, 1)
    elseif iseven(g.group.N) && 2i == g.group.N
        fill(cospi(g.n), 1, 1)
    else
        θ = 2 * i * g.n / g.group.N
        c, s = cospi(θ), sinpi(θ)
        [c -s; s c]
    end
irreptype(group::CyclicGroup, i) =
    if i == 0
        'R'
    elseif iseven(group.N) && 2i == group.N
        'R'
    else
        'C'
    end
irreps(group::CyclicGroup) = 0:div(group.N, 2)

irrep(group, i) =
    IrreducibleRepresentation(group, i, g -> irrepmat(g, i), irreptype(group, i))

struct FiberField{G,A}
    group::G
    x::A
end
