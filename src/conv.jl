"Steerable convolutional layer on ``\\mathbb{R}^2``."
struct R2Conv{FI,FO,IW,IB,C} <: Lux.AbstractExplicitLayer
    field_in::FI
    field_out::FO
    nin::Int
    nout::Int
    cin::Int
    cout::Int
    init_weight::IW
    init_bias::IB
    conv::C
end

function R2Conv(
    width,
    fieldtypes,
    channels;
    activation = identity,
    use_bias = true,
    init_weight = glorot_uniform,
    init_bias = zeros32,
    kwargs...,
)
    field_in, field_out = fieldtypes
    cin, cout = channels
    @assert field_in.gspace == field_out.gspace
    nin = size(field_in.representation, 1)
    nout = size(field_out.representation, 1)
    conv = Conv((width, width), cin * nin => cout * nout, activation; use_bias, kwargs...)
    R2Conv(field_in, field_out, nin, nout, cin, cout, init_weight, init_bias, conv)
end

function (c::R2Conv)(x, params, state)
    (; nin, nout, cin, cout, conv) = c
    @tensor weight[x, y, nin, cin, nout, cout] =
        state.weightbasis[x, y, nin, nout] * params.weight[d, nin, cin, nout, cout]
    if uses_bias(c)
        @tensor bias[nout, cout] = state.biasbasis[b, nout] * params.bias[b, nout, cout]
        params = (; weight, bias)
    else
        params = (; weight)
    end
    y, convstate = conv(x, params, state.convstate)
    y, (; state..., convstate)
end

uses_bias(c::R2Conv) = uses_bias(c.conv)
uses_bias(c::Conv{N,use_bias}) where {N,use_bias} = use_bias

function Lux.initialparameters(rng::AbstractRNG, c::R2Conv)
    (; nin, cin, nout, cout, field_out, init_weight, init_bias) = c
    (; gspace, representation) = field_out
    (; group) = gspace
    weightbasis = build_weightbasis(c)
    b = size(weightbasis, 5)
    weight = init_weight(rng, b, nin, cin, nout, cout)
    if uses_bias(c)
        # One bias degree of freedom for each trivial output representation
        ψ = trivial_representation(group)
        nbias = sum(==(ψ.freq), frequencies(representation))
        bias = init_bias(rng, nbias, nout, cout)
        (; weight, bias)
    else
        (; weight)
    end
end

function Lux.initialstates(rng::AbstractRNG, c::R2Conv)
    conv = Lux.initialstates(rng, c.conv)
    state = (; conv, weigthbasis = build_weightbasis(c))
    if uses_bias(c)
        state = (; state..., biasbasis = build_biasbasis(c))
    end
    state
end

function build_biasbasis(c::R2Conv)
    (; field_out) = c
    (; gspace, representation) = field_out
    (; group) = gspace
    ψ = trivial_representation(group)
    trivials = map(==(ψ.freq), frequencies(representation))
    basis(representation)[:, trivials]
end

function build_weightbasis(c::R2Conv)
    (; field_in, field_out, conv) = c
    (; kernel_size, dilation) = conv
    kernel_size, dilation = kernel_size[1], dilation[1]

    (; basis_filter, rings, σ, cutoffs) = compute_basis_params(; kernel_size, dilation)

    # Build coordinates
    points = get_grid_coords(2, kernel_size, dilation)
    x, y = points[1, :], points[2, :]
    radii = @. sqrt(x^2 + y^2)
    angles = @. atan(y, x)

    # Build radial basis
    basis = map(rings, σ, cutoffs) do ring, σ, cutoff
        radial_profile = GaussianRadialProfile(ring, σ).(radii)
        angular_basis = build_angular_basis(;
            gspace,
            angles,
            ρ_in = field_in.representation,
            ρ_out = field_out.representation,
            σ,
            maximum_frequency = cutoff,
        )
        radial_profile .* angular_basis
    end
    basis = stack(basis; dims = 4)
    _, sz... = size(basis)
    reshape(basis, kernel_size, kernel_size, sz...)
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
