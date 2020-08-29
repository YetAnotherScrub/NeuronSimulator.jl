# An example to test the Hodgkin-Huxley model implementation.
using NeuronSimulator
using BenchmarkTools
# Need to explicitly include these
using Unitful
using UnitfulRecipes

# Bring some units into a convenient namespace.
V = typeof(1.0u"V")
mV = typeof(1.0u"mV")
A = typeof(1.0u"A")
S = typeof(1.0u"S")
F = typeof(1.0u"F")
s = typeof(1.0u"s")
per_s = typeof(1.0u"1/s")

function hh_channels() :: Array{Chan, 1}
    # Original formulation from Wikipedia
    ma = v -> 1000 * 100u"1/s" * (0.025u"V" - v) / (1u"V" * (exp((0.025u"V" - v) / 0.01u"V") -1))
    mb = v -> 1000 * 4u"1/s" * exp(-v/0.018u"V")
    ha = v -> 1000 * 0.07u"1/s" * exp(-v/0.02u"V")
    hb = v -> 1000 * 1u"1/s" / (exp((0.03u"V" - v) / 0.01u"V") + 1)
    na = v -> 1000 * 10u"1/s" * (0.01u"V" - v) / (1u"V" * (exp((0.01u"V" - v) / 0.01u"V") - 1))
    nb = v -> 1000 * 0.125u"1/s" * exp(-v/0.08u"V")
    m = GatingVariable(0.0, ma, mb) # m
    h = GatingVariable(1.0, ha, hb) # h
    n = GatingVariable(0.0, na, nb) # n
    sodium = Chan(0.055u"V", 0.04u"S", [m, h], gvs -> gvs[1]^3 * gvs[2]) # Na, uses m^3 * h
    potassium = Chan(-0.077u"V", 0.035u"S", [n], gvs -> gvs[1]^4) # K, uses n ^4
    leak = Chan(-0.065u"V", 0.0003u"S", [], gvs -> 1) # Leak, no gvs

    return [sodium, potassium, leak]
end

function input_current(t :: s)
    return 2.7e-5u"A"
end

function vary_input_current(t :: s)
    if mod(ustrip(u"s", t), 0.05) < 0.0075
        return 3.0e-5u"A"
    else
        return 5.0e-6u"A"
    end
end

function test_constant()
    model = HHNeuronModel(hh_channels(), InitValues(-0.06u"V", 1.0e-6u"F"))
    param = Params((0.0u"s", 0.02u"s"), input_current)
    sol = simulate(model, param)
    plot_voltage(sol)
end

function test_vary()
    model = HHNeuronModel(hh_channels(), InitValues(-0.06u"V", 1.0e-6u"F"))
    param = Params((0.0u"s", 0.2u"s"), vary_input_current)
    sol = simulate(model, param)
    plot_voltage(sol)
end

# For timing with the @benchmark macro
function time(maxtime)
    model = HHNeuronModel(hh_channels(), InitValues(-0.06u"V", 1.0e-6u"F"))
    param = Params((0.0u"s", maxtime * 1u"s"), vary_input_current)
    sol = simulate(model, param)
end
