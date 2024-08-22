struct R2Conv{F1,F2,C} <: Lux.AbstractExplicitLayer
    field_in::F1
    field_out::F2
    cin::Int
    cout::Int
    conv::C
end

function R2Conv(width, fieldtypes, channels; activation, kwargs...)
    field_in, field_out = fields
    cin, cout = channels
    @assert field_in.gspace == field_out.gspace
    nin = size(field_in.representation, 1)
    nout = size(field_out.representation, 1)
    conv = Conv((width, width), cin * nin => cout * nout, activation; kwargs...)
    R2Conv(field_in, field_out, cin, cout, conv)
end

function Lux.initialparameters(rng::AbstractRNG, c::R2Conv)
    b = basisdimension(c)
end
