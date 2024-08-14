# LSP hack                                     #src
if false                                       #src
    include("../src/SteerableConvolutions.jl") #src
    using .SteerableConvolutions               #src
end                                            #src

# # Equivariance
#
# Introduction to GSpaces, FieldTypes, and FiberFields
# using the SteerableConvolutions.jl package.

using SteerableConvolutions

# ## GSpaces
#
# Create 2-dimensional GSpace with cyclic group of order 4.
G = CyclicGroup(4)
gspace = GSpace(G, 2)

# ## Field types
#
# Create a field type from a list of representations.
representations = [irrep(G, 0), irrep(G, 0), irrep(G, 1)]
fieldtype = FieldType(gspace, representations)

# ## Fiber fields
#
# Create a fiber field from a field type and data.
x = zeros(5, 5, 4)
x[1, 1, 1] = 4.7
x[2, 1, 3] = 2.2
x[1, 3, 4] = 8.3
f = FiberField(fieldtype, x)
f.x

# We now apply a group transform (rotation by 90 degrees):

g = G(1)
fnew = g * f

# Compare the original and transformed fields:

f.x

#-
fnew.x

# Note how all the 2D fields have been rotated by 90 degrees counter-clockwise.
# In addition, field 3 and 4 have been swapped, with a change of sign.
