module NeuronSimulator

using DifferentialEquations
using Plots
using Unitful
using UnitfulRecipes

include("NeuronModel.jl")
include("HH.jl")

export NeuronModel, simulate, voltage_function, get_variables, Params, plot_voltage,
    Chan, GatingVariable, InitValues, HHNeuronModel, plot_all, plot_gvs

end
