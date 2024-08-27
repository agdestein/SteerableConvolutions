# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Scratch
#
# Scratchpad test script.

using SteerableConvolutions
using SteerableConvolutions.Lux
using Random
using GLMakie

set_theme!(; colormap = :seaborn_icefire_gradient)

lines(cumsum(randn(100)))

σ = 0.1
r0 = 0.6
# r0 = 0
gp(r) = exp(-(r - r0)^2 / 2σ^2)

x = LinRange(-1, 1, 100)
y = x'

getr(x, y) = sqrt(x^2 + y^2)
getϕ(x, y) = atan(y, x)

r = getr.(x, y)
ϕ = getϕ.(x, y)

N = 4
ψ = @. sin(1 * N * ϕ)
ψ = @. cos(1 * N * ϕ)
# mat = gp.(r)
# mat = ψ
mat = @. gp(r) * ψ
heatmap(mat; colorrange = (-1, 1))

using SteerableConvolutions
G = CyclicGroup(8)
gspace = GSpace(G, 2)
ρ_in = Irrep(G, 0)
ρ_out = regular_representation(G)

kernel_size = 9
c = R2Conv(gspace, ρ_in => ρ_out, 1 => 1; kernel_size, use_bias = false);
rng = Random.Xoshiro()
params, states = Lux.setup(rng, c);

states.weightbasis
params.weight

w = states.weightbasis
w = reshape(w, kernel_size, kernel_size, :)

o = Observable(w[:, :, 1])
heatmap(o)
for i in 1:size(w, 3)
    o[] = w[:, :, i]
    sleep(0.1)
end

size(states.weightbasis)
size(params.weight)

size(states.biasbasis)
size(params.bias)

x = randn(50, 50, 1, 1, 1)
f = FiberField(gspace, ρ_in, x)

g = G(1)
cf, _ = c(f, params, states)
gcf = g * cf
gf = g * f
cgf, _ = c(gf, params, states)


gcf.x
cgf.x

n = 50
grid = LinRange(0, 1, n)
mask = @. (1 - grid^2) * (1 - (grid')^2)
f = @. sinpi(5 * grid) * mask
f = reshape(f, n, n, 1, 1, 1)
f = FiberField(gspace, ρ_in, f)

heatmap(f.x[:, :])

gf = g * f

heatmap(gf.x[:, :])

cf, _ = c(f, params, states)
gcf = g * cf
gf = g * f
cgf, _ = c(gf, params, states)

heatmap(cgf.x[:, :, 4])
heatmap(gcf.x[:, :, 4])

cgf.x[:, :, 4]
gcf.x[:, :, 4]

gcf.x[25, 24, 4]
cgf.x[25, 25, 4]

using LinearAlgebra
norm(gcf.x - cgf.x) / norm(gcf.x)


ρ = Irrep(G, 0)
ρ_hidden = regular_representation(G)

c = Chain(
    R2Conv(gspace, ρ => ρ_hidden, 1 => 1; activation = relu, kernel_size = 5, use_bias = true),
    R2Conv(gspace, ρ_hidden => ρ, 1 => 1; kernel_size = 5, use_bias = false),
)

params, states = Lux.setup(rng, c)

cf, _ = c(f, params, states)
gcf = g * cf
gf = g * f
cgf, _ = c(gf, params, states)

heatmap(cgf.x[:, :])
heatmap(gcf.x[:, :])
