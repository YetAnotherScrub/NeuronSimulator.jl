# Bring some units into a convenient namespace.
V = typeof(1.0u"V")
mV = typeof(1.0u"mV")
A = typeof(1.0u"A")
S = typeof(1.0u"S")
F = typeof(1.0u"F")
s = typeof(1.0u"s")
per_s = typeof(1.0u"1/s")

"Blank parent type for NeuronModels. Implemented model types should inherit from this."
abstract type NeuronModel end

"""
    voltage_function(model :: t where t <: NeuronModel)

Gives the differential equation function for a model, to be used in [`simulate`](@ref).
This must return a function of type suitable for an [`ODEProblem`](@ref) call.
The correct implementation depends on the model, and needs to be defined separately for each model type.

See Also: [`simulate`](@ref), [`ODEProblem`](@ref)
"""
function voltage_function(model :: t where t <: NeuronModel) :: Function
    throw(ArgumentError("voltage_function is not defined for this type of NeuronModel.\n
        Ensure that you are calling it on a concrete type for which voltage_function is defined."))
    return nothing
end

# Builds an array of the inital values of the variables for the DE.
# Should be defined for the specific NeuronModel subtype separately.
"""
    get_variables(model :: t where t <: NeuronModel)

Gives an array of the initial values of the variables for the differential equation of the given model, to be used in [`simulate`](@ref).
The correct implementation depends on the model, and needs to be defined separately for each model type.
Care should be taken to ensure that the order is correct, matching what is expected by the associated [`voltage_function`](@ref) implementation.

# Examples
```julia-repl
julia> get_variables(model)
4-element Array{Float64,1}:
 -0.06
  0.0
  1.0
  0.0
```

#See Also: [`simulate`](@ref), [`voltage_function`](@ref)
"""
function get_variables(model :: t where t <: NeuronModel) :: Array
    return []
end

"""
The input parameters for the simulation. Contains the timespan, a tuple in Unitful seconds,
and the input current, which is a function of time, of the form

    function input_current(t :: typeof(1.0u"s")) :: typeof(1.0u"A")

and must be written with awareness of the Unitful units.
Other arguments for the solver may be added in the future.

# Examples
```julia-repl
julia> Params((0.0u"s", 1.0u"s"), t -> 2.7e-5u"A")
Params((0.0 s, 1.0 s), var"#26#27"())
```

See also: [`simulate`](@ref)
"""
mutable struct Params
    tspan :: Tuple{s, s}
    input_current :: Function
    # Other arguments could be added here?
end

"""
    simulate(model, params :: Params) :: ODESolution{}

Simulates a given instance of a model, according to the parameters given, using the [`Tsit5`](@ref) algorithm.
Requires [`voltage_function`](@ref) and [`get_variables`](@ref) to be defined first.

# Examples
An example using a Hodgkin–Huxley implementation, output abbreviated.
```julia-repl
julia> sol = simulate(model, param)
retcode: Success
Interpolation: specialized 4th order "free" interpolation
t: 93-element Array{...}
    0.0 s
    ...
u: 93-element Array{Array{Float64,1},1}:
    [-0.06, 0.0, 1.0, 0.0]
    ...
```

See also: [`voltage_function`](@ref), [`get_variables`](@ref), [`ODESolution`](@ref)
"""
function simulate(model, param :: Params) :: ODESolution{}
    if !(model isa NeuronModel) throw(ArgumentError("The model must be a child of the abstract NeuronModel type.")) end
    f = voltage_function(model)
    u0 = get_variables(model)
    prob = ODEProblem(f, u0, param.tspan, param.input_current)
    sol = solve(prob, Tsit5())
    return sol
end

"Convenience function to plot the voltage of a simulation."
function plot_voltage(sol :: ODESolution{})
    plot(sol.t, map(first, sol.u), labels="V", xaxis="Time", yaxis="Voltage")
end

import Base.zero
# Hack to make DifferentialEquations work with an array with some Unitful and some raw values.
function zero(arr ::Array{Quantity{Float64,D,U} where U where D}) :: Array{Quantity{Float64,D,U} where U where D}
    return map(x -> x * 0.0, arr)
end
