abstract type Variable end

mutable struct GatingVariable <: Variable
    n0 :: Float64
    # dn/nt = a . V . (1 - n) - b . V . n
    a :: Function # A function of V
    b :: Function
    GatingVariable(n0, a, b) = new(n0, a, b)
    GatingVariable() = GatingVariable(0.0, x -> 0.0, x -> 0.0)
end

function foo() :: Function
    return (x -> throw(ArgumentError("voltage_function is not defined for this type of NeuronModel.\n
    Ensure that you are calling it on a concrete type for which voltage_function is defined.")))
end
