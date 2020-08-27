include("NeuronModel.jl")

"""
    GatingVariable(n0, α, β)

A structure to represent a gating variable of a Hodgkin–Huxley neuron, for use within a [`Chan`](@ref) struct.
`n0` is the initial value of the gating variable. `α` and `β` are functions of voltage, of the form

    function a(v :: typeof(1.0u"V")) :: typeof(1.0u"1/s")

and must be written with awareness of the Unitful units. They are applied in [`gv_function`](@ref)
for the general Hodgkin–Huxley equation for a gating variable:

    dn/nt = α(V) . (1 - n) - β(V) . n

# Examples
```julia-repl
julia> α = v -> 1000 * 0.07u"1/s" * exp(-v/0.02u"V")
#48 (generic function with 1 method)

julia> β = v -> 1000 * 0.125u"1/s" * exp(-v/0.08u"V")
#52 (generic function with 1 method)

julia> gv = GatingVariable(1.0, α, β)
GatingVariable(1.0, var"#50#51"(), var"#52#53"())
```

See also: [`Chan`](@ref), [`gv_function`](@ref)
"""
mutable struct GatingVariable
    n0 :: Float64
    α :: Function
    β :: Function
    GatingVariable(n0, α, β) = new(n0, α, β)
    GatingVariable() = GatingVariable(0.0, x -> 0.0, x -> 0.0)
end

"""
    Chan(Vi, g, gating_vars :: Array{GatingVariable, 1}, gv_mult)

A structure to represent a channel of a Hodgkin–Huxley neuron, for use within an [`HHNeuronModel`](@ref) struct.
`Vi` is the stable voltage for the channel, given in Unitful volts.
`g` is the conductance of the channel, given in Unitful siemens.
`gating_vars` is a (possibly empty) array of [`GatingVariable`](@ref) structs.
`gv_mult` is a function of the form

    function gv_mult(gvs :: Array{Float64, 1})

that describes how the gating variables will be applied when calculating the current. See the example below.
In the [`channel_current`](@ref) function, the gating variables are applied as a multiplier

    chan.gv_mult(gvs)

where `gvs` is an array of the present values of the gating variables, given in the same order as `gating_vars`.

# Examples
With defined `GatingVariable`s gv1 and gv2:
```julia-repl
julia> c = Chan(0.055u"V", 0.04u"S", [gv1, gv2], gvs -> gvs[1]^3 * gvs[2])
Chan(0.055 V, 0.04 S, GatingVariable[GatingVariable(...), GatingVariable(...)], var"#99#100"())
```
This will define a channel whose current at voltage v is given by

    c.g * gv1^3 * gv2 * (c.Vi - v)

See also: [`HHNeuronModel`](@ref), [`GatingVariable`](@ref), [`channel_current`](@ref)
"""
mutable struct Chan
    Vi :: V
    g :: S
    gating_vars :: Array{GatingVariable, 1}
    gv_mult :: Function
    Chan(vi, g, gvs, mult) = new(vi, g, gvs, mult)
    Chan() = Chan(0.0, 0.0, [], gvs -> 1)
end

"""
    InitValues(V0, C, I)

A structure of values for intialising a [`HHNeuronModel`](@ref).
`V0` is the starting voltage, given in Unitful volts.
`C` is the membrane capacitance, given in Unitful Farads.
`I` is the input current, which should be a function of time, of the form

    function input_current(t :: typeof(1.0u"s")) :: typeof(1.0u"A")

and must be written with awareness of the Unitful units.
"""
mutable struct InitValues
    V0 :: V # Starting voltage
    C :: F # Capacitance
    I :: Function # Input current, now a function of time
    InitValues(v, c, i) = new(v, c, i)
    InitValues() = InitValues(0.0, 0.0, 0.0)
end

"""
    HHNeuronModel(chans, inits)

A structure to represent a Hodgkin–Huxley neuron, for simulation with the [`simulate`](@ref) function.
Holds an array of channels as [`Chan`](@ref) structs, and the initialisation values as an [`InitValues`](@ref) struct.

See also: [`simulate`](@ref), [`Chan`](@ref), [`InitValues`](@ref)
"""
mutable struct HHNeuronModel <: NeuronModel
    chans :: Array{Chan, 1}
    init_values :: InitValues
    HHNeuronModel(arr, inits) = new(arr, inits)
    HHNeuronModel() = HHNeuronModel([], InitValues())
end

"""
    get_gvs(model :: HHNeuronModel) :: Tuple{Array{GatingVariable, 1}, Dict{Chan, Int}}

Creates a tuple `(gvs, indices)` from an [`HHNeuronModel`](@ref).
Mainly intended for internal use within [`voltage_function`](@ref) and [`get_variables`](@ref).
`gvs` is an array of all the [`GatingVariable`](@ref)s from all the channels of the model.
indices is a [`Dict`](@ref) which acts as a lookup table to find the [`GatingVariable`](@ref)s within `u` in [`voltage_function`](@ref).

See also: [`HHNeuronModel`](@ref), [`voltage_function`](@ref), [`get_variables`](@ref), [`GatingVariable`](@ref)
"""
function get_gvs(model :: HHNeuronModel) :: Tuple{Array{GatingVariable, 1}, Dict{Chan, Int}}
    gvs = []
    x = 2 # First index is 2 since 1 will be V.
    indices = Dict{Chan, Int}()
    for chan in model.chans
        append!(gvs, chan.gating_vars)
        if length(chan.gating_vars) > 0
            indices[chan] = x
            x += length(chan.gating_vars)
        end
    end
    return (gvs, indices)
end

Float64

# Given a channel, the gvs, the model and the values of the variables, return the current through that channel.
"""
    channel_current(chan :: Chan, u :: Array, indices :: Dict{Chan, Int}) :: typeof(1.0u"A")

Finds the current through a given [`Chan`](@ref).
`u` is the variables array for the overall differential equation.
`indices` is a [`Dict`](@ref), from [`get_gvs`](@ref), which acts as a lookup table to find the [`GatingVariable`](@ref)s within `u`.

# Examples
```julia-repl
julia> channel_current(model.chans[3], [-0.06u"V", 0.0, 1.0, 0.0], get_gvs(model)[2])
-1.500000000000001e-6 A
```

See also: [`voltage_function`](@ref), [`Chan`](@ref), [`get_gvs`](@ref), [`GatingVariable`](@ref)
"""
function channel_current(chan :: Chan, u :: Array, indices :: Dict{Chan, Int}) :: A
    gv_mult = 1.0 # The muliplier for the gvs, eg m^3 * h
    if length(chan.gating_vars) > 0
        gvs = u[indices[chan] : indices[chan] + length(chan.gating_vars) - 1]
    else
        gvs = []
    end
    mult = chan.gv_mult(gvs)
    # V is u[1]
    return uconvert(u"A", chan.g * mult * (chan.Vi - u[1]))
end

# Given a gv and the voltage, acts as the in-place DE for the gv.
"""
    gv_function(gv :: GatingVariable, v :: typeof(1.0u"V"), n) :: typeof(1.0u"1/s")

Differential equation function to find the rate of change of a [`GatingVariable`](@ref) `gv`, given
its current value `n` and membrane voltage `v`.

See also: [`voltage_function`](@ref), [`GatingVariable`](@ref)
"""
function gv_function(gv :: GatingVariable, v :: V, n) :: per_s
    return (gv.α(v) * (1 - n) - gv.β(v) * n)
end

function get_variables(model :: HHNeuronModel) :: Array
    u0 :: Array{Quantity{Float64,D,U} where U where D,1} = [model.init_values.V0]
    for gv in get_gvs(model)[1]
        push!(u0, gv.n0)
    end
    return u0
end

"""
    voltage_function(model :: HHNeuronModel)

Implementation for the [`HHNeuronModel`](@ref).
Applies the Hodgkin–Huxley equations to the voltage and the gating variables.

# Examples
```julia-repl
julia> voltage_function(hh_model)
(::var"#f!#18"{HHNeuronModel,Array{GatingVariable,1},Dict{Chan,Int64}}) (generic function with 1 method)
```

See also: [`HHNeuronModel`](@ref), [`simulate`](@ref)
"""
function voltage_function(model :: HHNeuronModel) # Type is f for f(V, p, t) = ...
    (gvs, indices) = get_gvs(model)
    function f!(du :: Array{}, u, p, t) :: Array{} # du is array of diff eqs. u is array of vars: V and gvs.
        # 1u"1/s" * ustript(u"V/s", # Forcefully remove the volts
        du[1] :: typeof(1.0u"V/s") = uconvert(u"V/s",
            (model.init_values.I(t) :: A
            + sum(map(chan -> channel_current(chan, u, indices), model.chans))) / model.init_values.C)
        for i = 2:length(du)
            du[i] :: per_s = gv_function(gvs[i-1], uconvert(u"V", u[1]), u[i])
        end
        return du
    end
    return f!
end

"Convenience function to plot all the values of an H-H simulation. Strips the voltage units."
function plot_all(sol :: ODESolution{}, labs="")
    arr = map(xs -> pushfirst!(xs[2:end], ustrip(u"V", xs[1])), sol.u)
    arr2 = permutedims(reshape(hcat(arr...), (length(arr[1]), length(arr))))
    plot(sol.t, arr2, labels=labs, xaxis="Time", yaxis="Value")
end

"Convenience function to plot all the gating variables of an H-H simulation without the voltage."
function plot_gvs(sol :: ODESolution{}, labs="")
    arr = map(xs -> xs[2:end], sol.u)
    arr2 = permutedims(reshape(hcat(arr...), (length(arr[1]), length(arr))))
    plot(sol.t, arr2, labels="V", xaxis="Time", yaxis="Value")
end
