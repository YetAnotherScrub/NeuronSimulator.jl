# NeuronSimulator.jl
This package provides a framework to facilitate the simulation of neuron models in Julia.
It is currently only at a "proof of concept" level, supporting only single neurons of the Hodgekin-Huxley model.
In future development, many more model types may be added, as well as support for networks of multiple neurons.
The package is dependent on several supporting packages: [**DifferentialEquations**](https://github.com/SciML/DifferentialEquations.jl), [**Unitful**](https://github.com/PainterQubits/Unitful.jl) and [**Plots**](https://github.com/JuliaPlots/Plots.jl).

## Core functionality
The base functionality of the system is defined in **NeuronModel.jl**. This outlines the core types and functions necessary for any neuron model to be simulated within this package. There are 3 steps to using this: *define* a model type, *instantiate* a model object, and *simulate* the resulting neuron.

*Defining* the model type is done using the components outlined in **NeuronModel.jl**. The suggested order for this is as follows:

1. Create a data type. The type should inherit from `NeuronModel`. The data type should represent the structure of the model in code. For example, if the model of the neuron has ion channels, then the model data type should have attached components to represent these channels. All physical units within the model should have attached **Unitful** units. To make a model data type more useful, try to keep it as flexible as possible, whilst still representing the core concepts of the specific model. For example with ion channels, it is preferable to store them as an array, whose size can be varied based on the implementation, rather than hard-coding a fixed number of channels.

2. Define the `get_variables` function for the new data type. Given an instance of the data type, this function will return an array of the starting values of all the variables needed for the differential equation that models the voltage across the neuron membrane. The first variable should be the voltage itself. Any other related variables whose behaviour follows a differential equation should also be in here. Again, remember to include **Unitful** units for the voltage, and any other variables which represent physical units. Be aware that if you want to create an array of mixed unit types, you may need to explicitly force a mixed Unitful type, such as `Array{Quantity{Float64,D,U} where U where D,1}`.

3. Define the `voltage_function` function for the new data type. This takes an instance of the model data type as input, and returns a function that represents the differential equation for the model. That differential equation function should be appropriate for creating an [`ODEProblem`](https://diffeq.sciml.ai/latest/types/ode_types/). The input `**u**` will be the array of the values of the variables, in the same order as was defined in `get_variables`. The return value should be an array of the rates of change of the variables. Again, these should have attached **Unitful** units, such a `u"V/s"`.

A specific example of this can be seen in **HH.jl**. Currently, this is the only model type that has been defined. However once a model data type is defined, if the design is suitably flexible, it can be used for many different implementations. Hence in the future, by defining the remaining standard models in common use, I hope to cover most of the likely use cases, whilst still supporting custom models for those that need them.

*Instantiating* the model object should then be relatively straightforward, given an understanding of the model data type. Simply give the constructors the correct physical values (with attached **Unitful** units as always) and the appropriate equations if applicable, and build the model object like any other `struct`.

Finally, *simulating* the neuron is then very straightforward. This is done using the `simulate` function from **NeuronModel.jl**. This function uses `get_variables` and `voltage_function` to create the differential equation, and solves it using the **DifferentialEquations** package. As well as the model object, it takes some parameters from a `Params` struct, namely the timespan of the simulation and the function for input current. It then runs the simulation and returns the `ODESolution` object. The results can be displayed using the `plot_voltage` function, which is just a convenience function to call the appropriate **Plots** package function - if you wish to display the results in another format, follow the appropriate documentation for that package.

## Developer
My name is Rainbow-Anthony Lilico, I'm a student at the University of Bristol studying computer science. This project is the primary subject of my final year thesis. Although it is currently in a very primitive stage, I hope to significantly expand the functionality in the future. If you are interesting in using or contributing to this software, feel free to contact me here.
