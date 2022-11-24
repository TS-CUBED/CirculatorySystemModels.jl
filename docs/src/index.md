```@meta
-
CurrentModule = CirculationModels
```

# CirculationModels

Documentation for [CirculationModels](https://github.com/TS-CUBED/CirculationModels.jl).

## An acausal modelling library for Circulation Models

CirculationModels.jl is an acausal modelling library for zero-dimensional, or _lumped parameter_ modelling of the circulatory system for [ModelingToolkit.jl (MTK)](https://github.com/SciML/ModelingToolkit.jl).

Integration into MTK and the [SciML](https://docs.sciml.ai/Overview/stable/) environment enables efficient solver methods ([DifferentialEquations.jl](https://diffeq.sciml.ai/latest/)), local and global parameter optimisation ([Optimization.jl](https://optimization.sciml.ai/stable/)), Sensitivity Analysis ([GlobalSensitivity.jl](https://gsa.sciml.ai/stable/)), and many other features in the SciML framework.

The acausal modelling approach in MTK was chosen to support a "if you can draw it you can model it" paradigm.

The main model states that are modelled throughout these models are:

- _volume flow rate_ (at nodes and through elements): $q\ [\mathrm{ml/s}]\ [\mathrm{cm^{3}/s}]$. Flow into an element is positive, flow out of an element is negative.
- _pressure_ (at nodes): $p\ [\mathrm{mm_{{Hg}}}]$.
- _pressure difference_ (over elements): $\Delta p\ [\mathrm{mm_{{Hg}}}]$. The pressure difference is following the usual fluid mechanical definition $\Delta p = p_{out} - p_{in}$. and is usually negative in flow direction!


### Units

There are many unit systems that are used in circulation models.
This modelling system uses the most common one, which uses $\mathrm{mm_{Hg}}$ for pressures and $\mathrm{ml/s}$ or $\mathrm{cm^3/s}$ for flow rates.

This is a variation of the $\mathrm{[g, cm, s]}$ system, which could be called $\mathrm{[g, cm, s, mm_{Hg}]}$ system.

Different model components are developed based on publications that use different unit systems. In those cases we attempted to keep the equations in the published system and do unit conversions transparently within the component function, while the outside API stays in the $\mathrm{[g, cm, s, mm_{Hg}]}$ system.

All model parameters are to be given in the $\mathrm{[g, cm, s, mm_{Hg}]}$  system unless otherwise specified in the component documentation.

### Main variables

The flow is modelled in terms of pressure $p$ and flow $q$. In many 0D
models of the circulation system these are replaced by the electrical
terms of voltage $v$ and current $i$. For this model we want to use the
physiologically relevant parameters. This also avoids the confusion of
using the same symbol $v$ to denote both, potential at a connection, and
the difference in potential over a component, as is commonly done in
electrical analogon models. To denote the pressure drop over a
component, this model uses the symbol $\Delta p$.

Time is a parameter in all the symbolic operations, so needs to be
defined as such (do not use `t` as a variable name elsewhere!). Parameters can be set using the `@parameters` macro.

```julia
@parameters t
```
