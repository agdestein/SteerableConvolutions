"GSpace."
struct GSpace{G}
    "Group."
    group::G

    "Dimension of base space."
    dim::Int
end

basespace_action((; group)::GSpace) = Irrep(group, 1)

function build_angular_basis(; gspace, angles, ρ_in, ρ_out, σ, maximum_frequency)
    (; group) = gspace

    # Representation stuff
    Q_in, Q_out = basis(ρ_in), basis(ρ_out)
    Q_in_inv = inv(Q_in)
    f_in, f_out = frequencies(ρ_in), frequencies(ρ_out)
    s_in = map(f -> size(Irrep(group, f), 1), f_in)
    s_out = map(f -> size(Irrep(group, f), 1), f_out)
    stop_in, stop_out = cumsum(s_in), cumsum(s_out)
    ind_in, ind_out = eachindex(f_in), eachindex(f_out)

    # Build full basis - with super agressive broadcasting
    b = map(Iterators.product(ind_in, ind_out)) do (i, j)
        m, n = f_out[j], f_in[i]
        A = Q_in_inv[stop_in[i]-s_in[i]+1:stop_in[i], :]
        B = Q_out[:, stop_out[j]-s_out[j]+1:stop_out[j]]
        ψm = Irrep(group, m)
        ψn = Irrep(group, n)

        # Resulting size: (npoint, β, α, dϕ),
        # where (α, β) is the size of the subkernel
        # and dϕ is the number of angular basis functions
        # for the irrep combination (ψm, ψn)
        κ = get_basis_block(ψm, ψn, angles; maximum_frequency)

        # Expand subkernel κ to full basis kernel k
        # Resulting size: (npoint, a, b, d)
        # where a is the size of ρ_in and b is the size of ρ_out
        s = size(κ)
        if s[4] == 0
            k = zeros(s[1], size(A, 2), size(B, 1), 0)
        else
            @tensor k[x, a, b, d] := B[b, β] * κ[x, β, α, d] * A[α, a]
        end

        k
    end
    b = reshape(b, :)
    b = cat(b...; dims = 4)
    reshape(b, size(b, 1), size(b, 2), size(b, 3), :)
end

GaussianRadialProfile(r0, σ) = r -> exp(-(r - r0)^2 / 2σ^2)

function get_basis_block(ψm, ψn, angles; maximum_frequency)
    @assert ψm.group == ψn.group
    get_basis_block(ψm.group, ψm, ψn, angles; maximum_frequency)
end

function get_basis_block(group::CyclicGroup, ψm, ψn, angles; maximum_frequency)
    N = order(group)
    M = maximum_frequency
    m, n = ψm.freq, ψn.freq
    ϕ = angles
    mscalar = m == 0 || 2 * m == N
    nscalar = n == 0 || 2 * n == N
    get_a(i, j, t) =
        if 2 * i == N ⊻ 2 * j == N
            (t + 1) / 2
        else
            t
        end
    if mscalar && nscalar
        nt = if 2 * m == N ⊻ 2 * n == N
            div(2 * M, N) - 1
        else
            div(M, N)
        end
        b1 = map(0:nt) do t
            a = get_a(m, n, t)
            @. cos(a * N * ϕ)
        end
        b2 = map(1:nt) do t
            a = get_a(m, n, t)
            @. sin(a * N * ϕ)
        end
        b = cat(b1..., b2...; dims = 2)
    elseif mscalar || nscalar
        i, j = mscalar ? (m, n) : (n, m)
        start, stop = if 2 * i == N
            -div(2 * (M + m), N) - 1, div(2 * (M - m), N) - 1
        else
            -div(M + m, N), div(M - m, N)
        end
        b1 = map(start:stop) do t
            a = get_a(i, j, t)
            s = @. sin((j + a * N) * ϕ)
            c = @. cos((j + a * N) * ϕ)
            cat(-s, c; dims = 2)
        end
        b2 = map(start:stop) do t
            a = get_a(i, j, t)
            s = @. sin((j + a * N) * ϕ)
            c = @. cos((j + a * N) * ϕ)
            cat(c, s; dims = 2)
        end
        b = cat(stack(b1), stack(b2); dims = 3)
    else
        b1 = map(-div(M + m - n, N):div(M - m + n, N)) do t
            s = @. sin((m - n + t * N) * ϕ)
            c = @. cos((m - n + t * N) * ϕ)
            cat(c, s, -s, c; dims = 2)
        end
        b2 = map(-div(M + m - n, N):div(M - m + n, N)) do t
            s = @. sin((m - n + t * N) * ϕ)
            c = @. cos((m - n + t * N) * ϕ)
            cat(-s, c, -c, -s; dims = 2)
        end
        b3 = map(-div(M + m + n, N):div(M - m - n, N)) do t
            s = @. sin((m + n + t * N) * ϕ)
            c = @. cos((m + n + t * N) * ϕ)
            cat(c, s, s, -c; dims = 2)
        end
        b4 = map(-div(M + m + n, N):div(M - m - n, N)) do t
            s = @. sin((m + n + t * N) * ϕ)
            c = @. cos((m + n + t * N) * ϕ)
            cat(-s, c, -c, -s; dims = 2)
        end
        b = cat(stack(b1), stack(b2), stack(b3), stack(b4); dims = 3)
    end
    sm, sn = size(ψm, 1), size(ψn, 1)
    reshape(b, length(ϕ), sm, sn, :)
end

"Fiber field."
struct FiberField{F,R,A}
    "Field type."
    gspace::F

    "Representation."
    representation::R

    "Field content."
    x::A

    function FiberField(gspace, ρ, x)
        @assert size(x, 1) == size(x, 2) "Domain must be square"
        new{typeof(gspace),typeof(ρ),typeof(x)}(gspace, ρ, x)
    end
end

function Base.:*(g::Element, f::FiberField)
    (; gspace, representation, x) = f
    (; group, dim) = gspace
    k, _, nfeat, s... = size(x)
    x = reshape(x, k, k, nfeat, :)

    A = basespace_action(gspace)(inv(g))
    @assert dim == 2 "Only 2D for now."
    R = representation(g)

    # Transform fibers
    @tensor y[i, j, a, sample] := R[a, b] * x[i, j, b, sample]

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
        i = sx == 1 ? (1:k) : (k:-1:1)
        i = Rx == 1 ? i : i'
        j = sy == 1 ? (1:k) : (k:-1:1)
        j = Ry == 1 ? j : j'
        I = CartesianIndex.(i, j)
        z = y[I, :, :]
    else
        # Get coordinates of the new grid (Cartesian)
        # Transform coordinates with inverse element
        # Interpolate old field to these coordinates
        grid = get_grid_coords(dim, k, 1)
        Agrid = A * grid
        bounds = extrema(grid)
        scalargrid = range(bounds..., k)
        i = @. (bounds[1] < Agrid[1, :] < bounds[2]) &&
           (bounds[1] < Agrid[2, :] < bounds[2])
        y = reshape(y, k, k, :)
        z = stack(1:size(y, 3)) do a
            field = y[:, :, a]
            p = linear_interpolation((scalargrid, scalargrid), field)
            za = zeros(k^2)
            @. za[i] = p(Agrid[1, i], Agrid[2, i])
            reshape(za, k, k)
        end
    end
    z = reshape(z, k, k, nfeat, s...)

    FiberField(gspace, representation, z)
end

# Tensors/einsums:
# - https://github.com/mcabbott/Tullio.jl
# - https://github.com/under-Peter/OMEinsum.jl
# - https://github.com/ITensor/ITensors.jl
# - https://github.com/Jutho/TensorKit.jl
# - https://github.com/Jutho/TensorOperations.jl
