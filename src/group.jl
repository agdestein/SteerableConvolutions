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

# (group::AbstractGroup)(n) = Element(group, n)
function Base.show(io::IO, e::Element)
    print(io, e.group)
    if e.n isa Tuple
        print(io, e.n)
    else
        print(io, "(", e.n, ")")
    end
end
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
struct Irrep{G,F} <: AbstractRepresentation
    "Group."
    group::G

    "Irrep frequency."
    freq::F
end

"General group representation. It is stored as a direct sum of irreps with a basis change."
struct Representation{F,A} <: AbstractRepresentation
    "List of irrep frequencies by order of appearance in direct sum."
    frequencies::F

    "Change of basis matrix."
    basis::A
end

function ((; frequencies, basis)::Representation)(g)
    S = directsum(map(i -> Irrep(g.group, i)(g), frequencies)...)
    basis * S * inv(basis)
end

"""
    istrivial(::Irrep)

Return `true` if the irrep is trivial.
"""
function istrivial end

"Get irrep frequencies used to block-diagonalize a representation."
function frequencies end

"Get change of basis matrix for block-diagonal representation decomposition."
function basis end

frequencies(ψ::Irrep) = [ψ.freq]
basis(ψ::Irrep) = ψ(one(ψ.group))

frequencies(ρ::Representation) = ρ.frequencies
basis(ρ::Representation) = ρ.basis

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
        error("irrep: unknown type $(type)")
    end

"Get the regular representation of a group."
function regular_representation(group)
    N = order(group)
    e = elements(group) |> vec
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
        ψ = Irrep(group, i)
        type = irreptype(ψ)
        m = map(characters, e) do c, g
            c * tr(ψ(inv(g)))
        end
        m = sum(m) / length(characters) / sum_of_squares_constituents(type)
        mint = round(Int, m)
        @assert mint ≈ m
        mint
    end
    ilist = vcat(fill.(i, multiplicities)...)
    irrs = map(i -> Irrep(group, i), ilist)
    P = directsum(irrs...)

    v = zeros(N)
    p = 0
    for (irr, m) in zip(i, multiplicities)
        s = size(Irrep(group, irr)(one(group)), 1)
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

"""
Cyclic group ``C_N`` of order ``N``.

# Example

```jldoctest
julia> G = CyclicGroup(8)
CyclicGroup(8)

julia> one(G)
CyclicGroup(8)(0)

julia> g, h = G(1), G(2)
(CyclicGroup(8)(1), CyclicGroup(8)(2))

julia> inv(g)
CyclicGroup(8)(7)

julia> g * h
CyclicGroup(8)(3)
```
"""
struct CyclicGroup <: AbstractFiniteGroup
    "Number of rotations in the group."
    N::Int
end

(group::CyclicGroup)(n) = Element(group, mod(n, group.N))
order(group::CyclicGroup) = group.N
Base.one(group::CyclicGroup) = group(0)
Base.:*(g::Element{CyclicGroup}, h::Element{CyclicGroup}) = g.group(g.n + h.n)
Base.inv(g::Element{CyclicGroup}) = g.group(-g.n)
elements(group::CyclicGroup) = group.(0:group.N-1)
function (ψ::Irrep{CyclicGroup})(g)
    (; group, freq) = ψ
    @assert group == g.group
    if 2 * freq > group.N
        error("Irrep{CyclicGroup}: 2 * freq must be less N")
    elseif freq == 0
        fill(1.0, 1, 1)
    elseif iseven(group.N) && 2 * freq == group.N
        fill(cospi(g.n), 1, 1)
    else
        rotmat(2π * freq * g.n / group.N)
    end
end
irreptype(ψ::Irrep{CyclicGroup}) =
    if ψ.freq == 0
        'R'
    elseif iseven(ψ.group.N) && 2 * ψ.freq == ψ.group.N
        'R'
    else
        'C'
    end
frequencies(group::CyclicGroup) = 0:div(group.N, 2)
istrivial(ψ::Irrep{CyclicGroup}) = ψ.freq == 0

"""
Dihedral group ``D_N``.

# Example

```jldoctest
julia> G = DihedralGroup(8)
DihedralGroup(8)

julia> one(G)
DihedralGroup(8)(false, 0)

julia> g, h = G(false, 5), G(true, 3)
(DihedralGroup(8)(false, 5), DihedralGroup(8)(true, 3))

julia> inv(g)
DihedralGroup(8)(false, 3)

julia> g * h
DihedralGroup(8)(true, 0)
```
"""
struct DihedralGroup <: AbstractFiniteGroup
    "Number of rotations in the group."
    N::Int
end

(group::DihedralGroup)(t::Tuple) = group(t...)
(group::DihedralGroup)(flip::Bool, n::Int) = Element(group, (flip, mod(n, group.N)))
order(group::DihedralGroup) = 2 * group.N
Base.one(group::DihedralGroup) = group(false, 0)
Base.:*(g::Element{DihedralGroup}, h::Element{DihedralGroup}) =
    g.group(g.n[1] ⊻ h.n[1], g.n[2] + (g.n[1] ? -1 : 1) * h.n[2])
Base.inv(g::Element{DihedralGroup}) = g.group(g.n[1], -(g.n[1] ? -1 : 1) * g.n[2])
elements(group::DihedralGroup) = group.(Iterators.product((false, true), 0:group.N-1))
function (ψ::Irrep{DihedralGroup})(g)
    (; group, freq) = ψ
    (; N) = group
    flip, i = freq
    @assert g.group == group
    if i == 0
        (flip && g.n[1] ? -1 : 1) * fill(1.0, 1, 1)
    elseif iseven(N) && 2i == N
        (flip && g.n[1] ? -1 : 1) * fill(cospi(g.n[2]), 1, 1)
    elseif flip && i > 0 && 2i < N
        fliprotmat(g.n[1], 2π * i * g.n[2] / N)
    else
        error("This combination is not an irrep")
    end
end
irreptype(::Irrep{DihedralGroup}) = 'R'
frequencies(group::DihedralGroup) = vcat(
    (false, 0),
    map(i -> (true, i), 0:div(group.N, 2)),
    fill((false, div(group.N, 2)), iseven(group.N)),
)
istrivial(ψ::Irrep{DihedralGroup}) = ψ.freq == (false, 0)
