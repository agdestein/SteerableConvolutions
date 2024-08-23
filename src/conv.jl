"Steerable convolutional layer on ``\\mathbb{R}^2``."
struct R2Conv{FI,FO,IW,IB,C} <: Lux.AbstractExplicitLayer
    field_in::FI
    field_out::FO
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
    activation,
    use_bias = true,
    init_weight = glorot_uniform,
    init_bias = zeros32,
    kwargs...,
)
    field_in, field_out = fields
    cin, cout = channels
    @assert field_in.gspace == field_out.gspace
    nin = size(field_in.representation, 1)
    nout = size(field_out.representation, 1)
    conv = Conv((width, width), cin * nin => cout * nout, activation; use_bias, kwargs...)
    R2Conv(field_in, field_out, cin, cout, init_weight, init_bias, conv)
end

uses_bias(c::R2Conv) = uses_bias(c.conv)
uses_bias(c::Conv{N,use_bias}) where {N,use_bias} = use_bias

function Lux.initialparameters(rng::AbstractRNG, c::R2Conv)
    (; field_out, init_weight, init_bias) = c
    (; gspace, representation) = field_out
    (; group) = gspace
    b = basisdimension(c)
    weight = init_weight(rng, b)
    if uses_bias(c)
        # One bias degree of freedom for each trivial output representation
        ψ = trivial_representation(group)
        nbias = sum(==(ψ.freq), frequencies(representation))
        bias = init_bias(rng, nbias)
        (; weight, bias)
    else
        (; weight)
    end
end

function Lux.initialstates(rng::AbstractRNG, c::R2Conv)
    state = Lux.initialstates(rng, c.conv)
    state = (; state..., weigthbasis = build_weightbasis(c))
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
    (; basis_filter, rings, σ, maximum_frequency) =
        compute_basis_params(; kernel_size = kernel_size[1], dilation = dilation[1])
    grid = get_grid_coords(2, kernel_size[1], dilation[1])
    basis = build_kernel_basis(
        gspace,
        field_in.representation,
        field_out.representation,
        σ,
        rings,
        maximum_frequency,
    )
    block_expansion = block_basisexpansion(; basis, grid, basis_filter)
end

get_grid_coords(d, kernel_size, dilation) =
    broadcast(1:kernel_size^d, (1:d)') do i, j
        p = div(i - 1, kernel_size^(j - 1)) % kernel_size
        p = dilation * (p - (kernel_size - 1) / 2)
        j > 1 ? -p : p
    end

function compute_basis_params(;
    kernel_size,
    dilation,
    cutoff = nothing,
    rings = nothing,
    σ = nothing,
    custom_basis_filter = nothing,
)
    width = dilation * (kernel_size - 1) / 2
    rmax = width * sqrt(2)

    if isnothing(rings)
        n_rings = ceil(Int, kernel_size / 2)
        rings = range(0, div((kernel_size - 1), 2), n_rings) * dilation
    end

    @assert all(r -> 0 ≤ r ≤ rmax, rings)

    if isnothing(σ)
        σ = vcat(fill(0.6, length(rings) - 1), 0.4)
        @. σ[rings≈0] = 0.005
    end

    if isnothing(cutoff)
        cutoff = r -> if r ≈ 0
            0.0
        else
            min(2 * r, r ≈ div(kernel_size, 2) ? 2 : 2 * r - (r + 1) % 2)
        end
    end

    maximum_frequency = round(Int, maximum(cutoff, rings))

    fco_filter = bandlimiter(cutoff)
    if isnothing(custom_basis_filter)
        basis_filter = fco_filter
    else
        basis_filter = d -> custom_basis_filter(d) && fco_filter(d)
    end

    (; basis_filter, rings, σ, maximum_frequency)
end

bandlimiter(cutoff) = a -> abs(a.irrep_frequency) ≤ cutoff(a.radius)

function block_basisexpansion(basis, grid, basis_filter)
    mask = basis_filter.(basis)
    @assert any(mask)
end
