"GSpace."
struct GSpace{G}
    "Group."
    group::G

    "Dimension of base space."
    dim::Int
end

basespace_action((; group)::GSpace) = Irrep(group, 1)

function build_kernel_basis(;
    gspace,
    kernel_size,
    ρ_in,
    ρ_out,
    σ,
    rings,
    maximum_frequency,
    dilation,
)
    (; group) = gspace

    # Representation stuff
    Q_in, Q_out = (ρ_in, ρ_out) |> basis
    Q_in_inv = inv(Q_in)
    f_in, f_out = (ρ_in, ρ_out) |> frequencies
    s_in = map(f -> size(Irrep(group, f), 1), f_in)
    s_out = map(f -> size(Irrep(group, f), 1), f_out)
    stop_in, stop_out = (s_in, s_out) |> cumsum
    ind_in, ind_out = (f_in, f_out) |> eachindex

    # Build coordinates
    grid = get_grid_coords(2, kernel_size, 1)
    x, y = grid[1, :], grid[2, :]
    r = @. sqrt(x^2 + y^2)
    ϕ = @. atan(y, x)

    # Build radial basis
    radial_basis = map(rings, σ) do ring, σ
        GaussianRadialProfile(ring, σ).(r)
    end
    radial_basis = reshape(radial_basis, kernel_size^2, 1, 1, 1, :)

    # Build full basis - with super agressive broadcasting
    b = map(Iterators.product(ind_in, ind_out)) do i, j
        m, n = fout[i], fin[j]
        A = Q_out[:, stop_out[i]-s_out[i]+1:stop_out[i]]
        B = Q_in_inv[stop_in[j]-s_in[j]+1:stop_in[j], :]
        ψm = Irrep(group, m)
        ψn = Irrep(group, n)

        # Resulting size: (kernel_size^2, α, β, dϕ),
        # where (α, β) is the size of the subkernel
        # and dϕ is the number of angular basis functions
        # for the irrep combination (ψm, ψn)
        block = get_basis_block(ψm, ψn, ϕ; maximum_frequency)

        # Resulting size: (kernel_size^2, α, β, dϕ, dr),
        # where dr is the number of radial basis functions (rings)
        κ = block .* radial_basis

        # Resulting size: (kernel_size^2, α, β, d),
        # where d is now the total number of basis kernels
        # for the irrep combination (ψm, ψn)
        κ = reshape(κ, kernel_size, kernel_size, s_out[i], s_in[j], :)

        # Expand subkernel κ to full basis kernel k
        # Resulting size: (kernel_size^2, a, b, d)
        # where a is the size of ρ_out and b is the size of ρ_in
        @tensor k[x, a, b, d] := A[a, α] * κ[x, α, β, d] * B[β, b]
        k
    end

    stack(b; dims = 5)
end

GaussianRadialProfile(r0, σ) = r -> exp(-(r - r0)^2 / 2σ^2)

function get_basis_block(ψm, ψn, ϕ; maximum_frequency)
    @assert ψm.group == ψn.group
    get_basis_block(ψm.group, ψm, ψm, maximum_frequency)
end

function get_basis_block(group::CyclicGroup, ψm, ψn, ϕ; maximum_frequency)
    N = order(group)
    M = maximum_frequency
    m, n = ψm.freq, ψn.freq
    mscalar = m == 0 || 2 * m == N
    nscalar = n == 0 || 2 * n == N
    get_a(i, j, t) =
        if 2 * i == N ⊻ 2 * j == N
            (t + 1) / 2
        else
            t
        end
    if mscalar && nscalar
        nt = if 2 * i == N ⊻ 2 * j == N
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
        b = cat(stack(b1), stack(b2); dims = 2)
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
            stack(-s, c)
        end
        b2 = map(start:stop) do t
            a = get_a(i, j, t)
            s = @. sin((j + a * N) * ϕ)
            c = @. cos((j + a * N) * ϕ)
            stack(c, s)
        end
        b = cat(stack(b1), stack(b2); dims = 3)
    else
        b1 = map(-div(M + m - n, N):div(M - m + n, N)) do t
            s = @. sin((m - n + t * N) * ϕ)
            c = @. cos((m - n + t * N) * ϕ)
            stack(c, s, -s, c)
        end
        b2 = map(-div(M + m - n, N):div(M - m + n, N)) do t
            s = @. sin((m - n + t * N) * ϕ)
            c = @. cos((m - n + t * N) * ϕ)
            stack(-s, c, -c, -s)
        end
        b3 = map(-div(M + m + n, N):div(M - m - n, N)) do t
            s = @. sin((m + n + t * N) * ϕ)
            c = @. cos((m + n + t * N) * ϕ)
            stack(c, s, s, -c)
        end
        b4 = map(-div(M + m + n, N):div(M - m - n, N)) do t
            s = @. sin((m + n + t * N) * ϕ)
            c = @. cos((m + n + t * N) * ϕ)
            stack(-s, c, -c, -s)
        end
        b = cat(stack(b1), stack(b2), stack(b3), stack(b4); dims = 3)
    end
    sm, sn = size(ψm, 1), size(ψn, 1)
    reshape(b, length(ϕ), sm, sn, :)
end

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
        i = @. (bounds[1] < Agrid[1, :] < bounds[2]) &&
           (bounds[1] < Agrid[2, :] < bounds[2])
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
