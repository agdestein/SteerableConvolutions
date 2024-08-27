"Steerable convolutional layer on ``\\mathbb{R}^2``."
struct R2Conv{G,RI,RO,IW,IB,C} <: Lux.AbstractExplicitLayer
    gspace::G
    ρ_in::RI
    ρ_out::RO
    c_in::Int
    c_out::Int
    init_weight::IW
    init_bias::IB
    conv::C
end

function R2Conv(
    gspace,
    representations,
    channels;
    kernel_size,
    activation = identity,
    use_bias = true,
    init_weight = glorot_uniform,
    init_bias = zeros32,
    kwargs...,
)
    ρ_in, ρ_out = representations
    c_in, c_out = channels
    n_in, n_out = size(ρ_in, 1), size(ρ_out, 1)
    conv = Conv(
        (kernel_size, kernel_size),
        c_in * n_in => c_out * n_out,
        activation;
        use_bias,
        kwargs...,
    )
    R2Conv(gspace, ρ_in, ρ_out, c_in, c_out, init_weight, init_bias, conv)
end

function (c::R2Conv)(x, params, states)
    (; conv) = c
    @tensor weight[kx, ky, nin, cin, nout, cout] :=
        states.weightbasis[kx, ky, nin, nout, d] * params.weight[d, cin, cout]
    kx, ky, n_in, c_in, n_out, c_out = size(weight)
    weight = reshape(weight, kx, ky, n_in * c_in, n_out * c_out)
    if uses_bias(c)
        @tensor bias[nout, cout] := states.biasbasis[nout, b] * params.bias[b, cout]
        bias = reshape(bias, :)
        params = (; weight, bias)
    else
        params = (; weight)
    end
    x = reshape(x, size(x, 1), size(x, 2), n_in * c_in, :)
    x, convstates = conv(x, params, states.convstates)
    x = reshape(x, size(x, 1), size(x, 2), n_out, c_out, :)
    x, (; states..., convstates)
end

uses_bias(c::R2Conv) = uses_bias(c.conv)
uses_bias(c::Conv{N,use_bias}) where {N,use_bias} = use_bias

function Lux.initialparameters(rng::AbstractRNG, c::R2Conv)
    (; gspace, ρ_in, ρ_out, c_in, c_out, init_weight, init_bias) = c
    (; group) = gspace
    weightbasis = build_weightbasis(c)
    d = size(weightbasis, 5)
    weight = init_weight(rng, d, c_in, c_out)
    if uses_bias(c)
        # One bias degree of freedom for each trivial output representation
        nbias = sum(f -> istrivial(Irrep(group, f)), frequencies(ρ_out))
        bias = init_bias(rng, nbias, c_out)
        (; weight, bias)
    else
        (; weight)
    end
end

function Lux.initialstates(rng::AbstractRNG, c::R2Conv)
    convstates = Lux.initialstates(rng, c.conv)
    states = (; convstates, weightbasis = build_weightbasis(c))
    if uses_bias(c)
        states = (; states..., biasbasis = build_biasbasis(c))
    end
    states
end

function build_biasbasis(c::R2Conv)
    (; gspace, ρ_out) = c
    (; group) = gspace
    irreps = map(f -> Irrep(group, f), frequencies(ρ_out))
    trivials = map(ψ -> istrivial(ψ), irreps)
    s = map(ψ -> size(ψ, 1), irreps)
    stops = cumsum(s)
    # Trivials are 1x1
    basis(ρ_out)[:, stops[trivials]]
end

function build_weightbasis(c::R2Conv)
    (; gspace, ρ_in, ρ_out, conv) = c
    (; kernel_size, dilation) = conv
    kernel_size, dilation = kernel_size[1], dilation[1]

    (; basis_filter, rings, σ, cutoffs) = compute_basis_params(; kernel_size, dilation)

    # Build coordinates
    points = get_grid_coords(2, kernel_size, dilation)
    x, y = points[1, :], points[2, :]
    radii = @. sqrt(x^2 + y^2)
    angles = @. atan(y, x)

    # Build radial basis
    b = map(rings, σ, cutoffs) do ring, σ, cutoff
        radial_profile = GaussianRadialProfile(ring, σ).(radii)
        angular_basis = build_angular_basis(;
            gspace,
            angles,
            ρ_in,
            ρ_out,
            σ,
            maximum_frequency = cutoff,
        )
        radial_profile .* angular_basis
    end
    b = cat(b...; dims = 4)
    s = size(b)
    reshape(b, kernel_size, kernel_size, s[2], s[3], :)
end

get_grid_coords(d, kernel_size, dilation) =
    broadcast(1:d, (1:kernel_size^d)') do i, j
        p = div(j - 1, kernel_size^(i - 1)) % kernel_size
        p = dilation * (p - (kernel_size - 1) / 2)
    end

function compute_basis_params(;
    kernel_size,
    dilation,
    rings = nothing,
    cutoffs = nothing,
    σ = nothing,
    basis_filter = nothing,
)
    width = dilation * (kernel_size - 1) / 2
    rmax = width * sqrt(2)

    if isnothing(rings)
        n_rings = ceil(Int, kernel_size / 2)
        rings = range(0, div((kernel_size - 1), 2), n_rings) * dilation
    end

    @assert all(r -> 0 ≤ r ≤ rmax, rings)

    if isnothing(cutoffs)
        cutoffs = map(rings) do r
            if r ≈ 0
                0
            else
                round(Int, min(2 * r, r ≈ div(kernel_size, 2) ? 2 : 2 * r - (r + 1) % 2))
            end
        end
    end

    if isnothing(σ)
        σ = vcat(fill(0.6, length(rings) - 1), 0.4)
        @. σ[rings≈0] = 0.005
    end

    maximum_frequency = maximum(cutoffs)

    # fco_filter = bandlimiter(cutoff)
    # if isnothing(custom_basis_filter)
    #     basis_filter = fco_filter
    # else
    #     basis_filter = d -> custom_basis_filter(d) && fco_filter(d)
    # end

    (; basis_filter, rings, σ, cutoffs, maximum_frequency)
end

# bandlimiter(cutoff) = a -> abs(a.ψ.freq) ≤ cutoff(a.radius)
