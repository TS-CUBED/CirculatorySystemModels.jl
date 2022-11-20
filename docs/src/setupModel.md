```@meta
CurrentModule = CirculationModels
```

## Setting up a simple model

CirculationModels.jl is an acausal modelling system. Models are composed from _elements_, which are _connected_.

Each _connection_ will be set up to fulfill the Kirchoff laws:

- all pressures at a connection are equal.
- all flows in and out of the connection add up to zero.
