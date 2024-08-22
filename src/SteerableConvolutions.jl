"""
Steerable convolutional layers.

## Exports

The following symbols are exported by SteerableConvolutions:

$(EXPORTS)
"""
module SteerableConvolutions

using DocStringExtensions
using EnumX
using LinearAlgebra
using Lux
using TensorOperations

# Docstring templates
@template MODULES = """
    $(DOCSTRING)

    # Exports

    $(EXPORTS)
    """
@template (FUNCTIONS, METHODS, MACROS) = """
    $TYPEDSIGNATURES

    $DOCSTRING

    # Methods

    $METHODLIST
    """
@template TYPES = """
    $TYPEDEF

    $DOCSTRING

    # Fields

    $FIELDS
    """

"$LICENSE"
license = "MIT"

include("utils.jl")
include("group.jl")
include("field.jl")

export AbstractGroup,
    AbstractFiniteGroup,
    AbstractRepresentation,
    RotationGroup,
    CyclicGroup,
    DihedralGroup,
    Element,
    FieldType,
    FiberField,
    GSpace,
    Irrep,
    IrrepType,
    Representation,
    basis,
    elements,
    frequencies,
    istrivial,
    regular_representation,
    ⊕,
    directsum

end
