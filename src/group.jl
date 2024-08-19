"Abstract group."
abstract type AbstractGroup end

"Abstract finite group."
abstract type AbstractFiniteGroup <: AbstractGroup end

"""
Group element.
Group elements can be multiplied with `*`.
To build group elements, call `group(n)` instead of `Element(group, n)`.
"""
struct Element{G,N}
    "Group to which the element belongs."
    group::G

    "Group element content."
    n::N
end

"Get elements of a finite group."
function elements end

"Get irrep of frequency `i`."
function irrep end

(group::AbstractGroup)(n) = Element(group, n)
Base.show(io::IO, e::Element) = print(io, e.group, "(", e.n, ")")
Base.one(e::Element) = one(e.group)
Base.:(^)(e::Element, n) = prod(fill(e, n); init = one(e))

"""
Abstract group representation.

A linear representation

```math
    \\rho : G \\to \\operatorname{GL}(\\mathbb{R}^d)
```

maps group elements ``g`` to invertible matrices ``\\rho(g)``
of size ``d \\times d``.
The representation must be such that

```math
    \\rho(gh) = \\rho(g) \\rho(h) \\quad \\forall (g, h) \\in G^2.
```
"""
abstract type AbstractRepresentation end

"Irreducible group representation (\"irrep\")."
struct IrreducibleRepresentation{G,A} <: AbstractRepresentation
    "Group."
    group::G

    "Irrep ID."
    id::Int

    "Irrep matrix constructor. Can be called as `mat(g)`."
    mat::A

    "Irrep type, can be `'R'`, `'C'`, or `'Q'`."
    type::Char
end

"General group representation. It is stored as a direct sum of irreps with a basis change."
struct Representation{A} <: AbstractRepresentation
    "List of irrep frequencies by order of appearance in direct sum."
    frequencies::Vector{Int}

    "Change of basis matrix."
    basis::A
end

(ψ::IrreducibleRepresentation)(g) = ψ.mat(g)
function ((; frequencies, basis)::Representation)(g)
    i = directsum(map(i -> irrepmat(g, i), frequencies)...)
    basis * i * inv(basis)
end

"Get irrep frequencies used to block-diagonalize a representation."
function frequencies end

"Get change of basis matrix for block-diagonal representation decomposition."
function basis end

frequencies(r::IrreducibleRepresentation) = [r.id]
basis(r::IrreducibleRepresentation) = one(r.mat(one(r.group)))

frequencies(r::Representation) = r.frequencies
basis(r::Representation) = r.basis

"Direct sum of matrices or representations."
function directsum end
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
    Representation(vcat(frequencies.(r)...), directsum(basis.(r)...))

"Alias for [`directsum`](@ref)."
const ⊕ = directsum

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

"Get the regular representation of a group."
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

    i = frequencies(group)
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

"Rotation group in the plane."
struct RotationGroup <: AbstractGroup end

(G::RotationGroup)(angle) = Element(G, mod2pi(angle))
Base.one(G::RotationGroup) = G(0)
Base.:*(g::Element{RotationGroup}...) = g[1].group(sum(getfield.(g, :angle)))
Base.inv(g::Element{RotationGroup}) = g.group(-g.angle)

"Cyclic group ``C_N`` of order ``N``."
struct CyclicGroup <: AbstractFiniteGroup
    "Number of rotations in the group."
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
frequencies(group::CyclicGroup) = 0:div(group.N, 2)

irrep(group, i) =
    IrreducibleRepresentation(group, i, g -> irrepmat(g, i), irreptype(group, i))
