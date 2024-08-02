"""
Steerable convolutional layers.

## Exports

The following symbols are exported by SteerableConvolutions:

$(EXPORTS)
"""
module SteerableConvolutions

using DocStringExtensions
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

end
