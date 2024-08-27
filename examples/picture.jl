using SteerableConvolutions
using SteerableConvolutions.Lux
using Random
using CairoMakie
using Images, FileIO
using Statistics

raw = load(joinpath(@__DIR__, "SC2.jpg"))
size(raw)

img = raw[1:3000, 1:3000]

i = 1:3000

mask(x, y) = (x - 1500)^2 + (y - 1500)^2 < 1500^2

img = @. img * mask(i, i')

img = imresize(img, ratio=1/6)
size(img)
k = size(img, 1)
img
G = CyclicGroup(12)
gspace = GSpace(G, 2)
ψ = Irrep(G, 0)
ρ = regular_representation(G)

channels = channelview(img)
x = cat(
    Float32.(channels[1, :, :]),
    Float32.(channels[2, :, :]),
    Float32.(channels[3, :, :]);
    dims = 4,
)
x = reshape(x, size(x)..., 1)

function rgb2fiberfield(img)
    channels = channelview(img)
    x = cat(
        Float32.(channels[1, :, :]),
        Float32.(channels[2, :, :]),
        Float32.(channels[3, :, :]);
        dims = 4,
    )
    x = reshape(x, size(x)..., 1)
    FiberField(gspace, ψ, x)
end

function fiberfield2rgb(f)
    (; x) = f
    kx, ky, _ = size(x)
    x = reshape(f.x, kx, ky, 3)
    RGB.(x[:, :, 1], x[:, :, 2], x[:, :, 3])
end

f = rgb2fiberfield(img)

g = G(6)

gf = g * f |> f32

gimg = fiberfield2rgb(gf)

round.(Int, ψ(G(3)))
round.(Int, ρ(G(7)))

SteerNN = let
    kernel_size = 9
    # chain = Chain(
    #     R2Conv(gspace, ψ => ψ, 3 => 3; kernel_size, use_bias = false, pad = SamePad()),
    #     # R2Conv(gspace, ρ => ψ, 5 => 3; kernel_size, use_bias = false, pad = SamePad()),
    # )
    chain = Chain(
        R2Conv(gspace, ψ => ρ, 3 => 5; activation = relu, kernel_size, use_bias = true, pad = SamePad()),
        # R2Conv(gspace, ρ => ψ, 5 => 5; activation = relu, kernel_size, use_bias = true, pad = SamePad()),
        R2Conv(gspace, ρ => ψ, 5 => 3; kernel_size, use_bias = false, pad = SamePad()),
    )
    rng = Random.Xoshiro(6)
    params, states = Lux.setup(rng, chain) |> f32 |> gpu

    function SteerNN(f)
        x = f.x |> f32 |> gpu
        y, _ = chain(x, params, states)
        y ./ std(y)
        FiberField(gspace, ψ, y) |> cpu
    end
end

CNN = let
    kernel_size = 9
    chain = Chain(
        Conv((kernel_size, kernel_size), 3 => 5, relu; pad = SamePad()),
        # Conv((kernel_size, kernel_size), 20 => 15, relu; pad = SamePad()),
        Conv((kernel_size, kernel_size), 5 => 3; use_bias = false, pad = SamePad()),
    )
    rng = Random.Xoshiro(0)
    params, states = Lux.setup(rng, chain) |> f32 |> gpu

    function CNN(f)
        x = f.x |> f32 |> gpu
        x = reshape(x, size(x, 1), size(x, 2), 3, 1)
        y, _ = chain(x, params, states)
        y ./ std(y)
        y = reshape(y, size(y, 1), size(y, 2), 1, 3, 1) |> cpu
        FiberField(gspace, ψ, y)
    end
end

function compare(NN, f, g)
    NNf = NN(f)
    NNgf = NN(g * f)
    gNNf = g * NNf
    hcat(
        NNf |> fiberfield2rgb,
        gNNf |> fiberfield2rgb,
        NNgf |> fiberfield2rgb,
    )
end

compare(SteerNN, f, g)
compare(CNN, f, g)
