# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Scratch
#
# Scratchpad test script.

using SteerableConvolutions

G = CyclicGroup(4)
gspace = GSpace(G, 2)
fieldtype = FieldType(gspace, [irrep(G, 1)])
x = zeros(5, 5, 2)
x[2, 1, 1] = 2.2
x[1, 3, 2] = 8.3
f = FiberField(fieldtype, x)
f.x

g = Element(G, 2)

fnew = g * f

f.x
fnew.x
