using SteerableConvolutions
using SteerableConvolutions.Lux
using Random
using CairoMakie

G = CyclicGroup(8)
gspace = GSpace(G, 2)
ρ_in = Irrep(G, 0)
# ρ_out = Irrep(G, 0) ⊕ regular_representation(G)
ρ_out = regular_representation(G)
kernel_size = 9
c = R2Conv(gspace, ρ_in => ρ_out, 1 => 1; kernel_size, use_bias = false);
rng = Random.Xoshiro()
params, states = Lux.setup(rng, c);

states.weightbasis |> size
params.weight |> size

w = reshape(states.weightbasis, kernel_size, kernel_size, :)
o = Observable(w[:, :, 1])
fig = heatmap(o; figure = (; size = (600, 600)), axis = (; aspect = DataAspect()))
# for i in 1:size(w, 3)
anim = record(fig, joinpath(@__DIR__, "basis.mp4"), 1:size(w, 3); framerate = 20) do i
    o[] = w[:, :, i]
end

hm(x) = heatmap(x; figure = (; size = (600, 600)), axis = (; aspect = DataAspect()))

let
    i = 9
    j = 23
    states.weightbasis[:, :, 1, i, j] |> hm
end
