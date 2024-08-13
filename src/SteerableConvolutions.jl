"""
Steerable convolutional layers.

## Exports

The following symbols are exported by SteerableConvolutions:

$(EXPORTS)
"""
module SteerableConvolutions

using DocStringExtensions
using LinearAlgebra
using Lux

# Docstring templates
@template MODULES = """
                   $(DOCSTRING)

                   ## Exports

                   $(EXPORTS)
                   """
@template (FUNCTIONS, METHODS) = """
                                 $TYPEDSIGNATURES

                                 $DOCSTRING
                                 """
@template TYPES = """
                  $TYPEDEF

                  $DOCSTRING

                  ## Fields

                  $FIELDS
                  """

"$LICENSE"
license = "MIT"

include("group.jl")

export AbstractGroup,
    AbstractFiniteGroup,
    Rotation,
    CyclicGroup,
    Element,
    IrreducibleRepresentation,
    Representation,
    elements,
    irrepmat,
    irrep,
    irreps,
    regular_representation,
    ⊕,
    directsum

end
