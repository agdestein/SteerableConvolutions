abstract type AbstractGroup end
abstract type AbstractFiniteGroup <: AbstractGroup end
abstract type AbstractRepresentation end
struct IrreducibleRepresentation{G,A} <: AbstractRepresentation
    id::Int
    mat::A
    type::Char
    IrreducibleRepresentation(::Type{G}, id, mat::A, type) where {G,A} =
        new{G,A}(id, mat, type)
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
basis(r::IrreducibleRepresentation{G}) where {G} = one(r.mat(one(G)))

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

function regular_representation(G)
    N = order(G)
    e = elements(G)
    representations = map(e) do g
        r = zeros(N, N)
        for (j, h) in enumerate(e)
            i = findfirst(==(h * g), e)
            r[i, j] = 1
        end
        r
    end
    characters = tr.(representations)

    i = irreps(G)
    multiplicities = map(i) do i
        (; type) = irrep(G, i)

        m = map(characters, e) do c, g
            c * tr(irrepmat(inv(g), i))
        end
        m = sum(m) / length(characters) / sum_of_squares_constituents(type)
        mint = round(Int, m)
        @assert mint ≈ m
        mint
    end
    ilist = vcat(fill.(i, multiplicities)...)
    irrs = map(i -> irrep(G, i), ilist)
    P = directsum(irrs...)

    v = zeros(N)
    p = 0
    for (irr, m) in zip(i, multiplicities)
        s = size(irrepmat(one(G), irr), 1)
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

struct CyclicGroup{N} <: AbstractFiniteGroup
    n::Int
    CyclicGroup{N}(n) where {N} = new{N}(mod(n, N))
end

order(::Type{CyclicGroup{N}}) where {N} = N
Base.one(G::Type{CyclicGroup{N}}) where {N} = G(0)
Base.:*(g::CyclicGroup{N}, h::CyclicGroup{N}) where {N} = CyclicGroup{N}(g.n + h.n)
Base.inv(g::CyclicGroup) = typeof(g)(-g.n)
elements(G::Type{CyclicGroup{N}}) where {N} = G.(0:N-1)
irrepmat(g::CyclicGroup{N}, i) where {N} =
    if 2i > N
        error("irrepmat: 2i must be less N")
    elseif i == 0
        fill(1.0, 1, 1)
    elseif iseven(N) && 2i == N
        fill(cospi(g.n), 1, 1)
    else
        θ = 2 * i * g.n / N
        c, s = cospi(θ), sinpi(θ)
        [c -s; s c]
    end
irreptype(::Type{CyclicGroup{N}}, i) where {N} =
    if i == 0
        'R'
    elseif iseven(N) && 2i == N
        'R'
    else
        'C'
    end
irreps(::Type{CyclicGroup{N}}) where {N} = 0:N÷2

irrep(G, i) = IrreducibleRepresentation(G, i, g -> irrepmat(g, i), irreptype(G, i))
