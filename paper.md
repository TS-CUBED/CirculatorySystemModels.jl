---
title: CirculationModels.jl - A ModelingToolkit Library for 0D-Lumped-Parameter Models of the Cardiovascular Circulation
tags:
  - Julia
  - ModelingToolkit
  - cardio-vascular
  - lumped parameter
  - circulation
  - acausal
  - patient-specific
authors:
  - name: Torsten Schenkel
    orcid: 0000-0001-5560-1872 
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Harry Saxton
    orcid: 0000-0001-7433-6154

    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Department of Engineering and Mathematics, Sheffield Hallam University, UK
   index: 1
 - name: Materials and Engineering Research Institute MERI, Sheffield Hallam University, UK
   index: 2
bibliography: paper.bib

---

# Summary

Within the realm of circulatory mechanics, lumped parameter (0D)
modelling [@shi2011review] offers the unique ability to examine both
cardiac function and global hemodynamics within the context of a single
model. Due to the interconnected nature of the cardiovascular system,
being able to quantify both is crucial for the effective prediction and
diagnosis of cardiovascular diseases. Lumped parameter modelling derives
one of its main strengths from the minimal computation time required to
solve ODEs and algebraic equations. Furthermore, the relatively simple
structure of the model allows most personalized simulations to be
automated. Meaning the ability to embed these lumped parameter models
into a clinical workflow could one day become trivial
[@bozkurt2022patient; @holmes2018clinical].

_CirculationModels.jl_ is a [Julia](https://www.julialang.org) [@bezanson2017julia] 
package, built on the acausal modelling framework provided by _ModelingToolkit.jl_
[@ma2021modelingtoolkit], containing all the common elements plus more
needed for effective and realistic lumped parameter modelling. Currently
_CirculationModels.jl_ supports common elements such as a capacitor,
resistor, inductance and diodes [@westerhof2010snapshots], which act as
simple valve functions. We also make extensions to the common elements
to include constant compliance chambers, non-linear and Poiseuille
resistances [@pfitzner1976poiseuille]. Plus the Double-Hill, and Shi
activation functions which are used as the cardiac driving chamber
elastance's
[@stergiopulos1996determinants; @korakianitis2006numerical].
We also include non-linear valve functions from Shi, and Mynard
[@korakianitis2006numerical; @mynard2012simple].
Alongisde individual components we also have created a collection of sub
compartments including, full circulatory, systemic, pulmonary and heart
models. We then also break down these full systems into collections of
elements such as the famous Windkessel models [@westerhof2009arterial]
to give the user full control over their modelling.

Users can easily add new elements to _CirculationModels.jl_ 
using _ModelingToolkit.jl_ functions.

# Statement of need

Lumped parameter modelling has become an essential part of contributing
strongly to our understanding of circulatory physiology. Lumped
parameter models have been used in many different contexts across
cardiovascular medicine such as aortic valve stenosis, neonatal
physiology and detection of coronary artery disease
[@laubscher2022dynamic; @sepulveda2022openmodelica; @dash2022non]. There
already exists some popular packages within both Simulink and Modelica
such as the "Cardiovascular library for lumped-parameter modeling"
[@rosalia2021object] and "Physiolibary" [@matejak2014physiolibrary]
respectively. These languages operate a block orientated "drag and
drop" approach with Modelica been the common choice due to it's acausal
modelling approach [@schweiger2020modeling]. Other languages, based on
XML, exist for lumped parameter modelling such as CellML
[@cuellar2003overview] and SBML [@hucka2003systems], while these are
great for exchanging models they are often difficult to implement and
model analysis is limited. A common theme within all current lumped 
parameter modelling software is the systems inability to deal with
complex event handling and non-linear components. 
Being based on _ModelingToolkit.jl_, _CirculationModels.jl_ overcomes
these limitations by leveraging the wider _SciML_ framework.

_CirculationModels.jl_
provides the Julia community with a quick and effective way to perform
lumped parameter modelling, being the first within the field to leverage
both multiple dispatch and JIT compilation. As a result of Julia's
architecture, as the complexity of the model increases the model
analysis time does not become unreasonable. 

Other packages exist which
allow users to import models from other frameworks, CellMLToolkit.jl
[@CellMLToolKit] and OpenModelica.jl [@tinnerholm2022modular]. 
_CirculationModels.jl_ goes beyond these by providing a lumped parameter modelling library
with seamless integration with the SciML framework 
[@Dixit2022; @rackauckas2020universal; @rackauckas2017differentialequations]
which allows for extensive and efficient model analysis.
Since both the modelling library and the framework it is built on are pure Julia,
new components can be developed in a transparent and consistent manner.

Using Julia, _CirculationModels.jl_ models compute significantly faster
than models implemented in Matlab or Python, and at the same speed as 
specialised C-code (\autoref{tbl:benchmarks}). This allows the models to run in real-time,
and opens the possibility of global parameter optimisation, global sensitivity analysis.

The modular, acausal approach also allows quick and straightforward changes to the
topology of the model, which can be automated using meta-programming techniques.



# Example

Validation and benchmarking was performed on a full, 4-chamber, model of the circulation system
proposed by [@korakianitis2006numerical] (\autoref{fig-shi-diagram}). This model was previously implemented in CellML [@Shi2018],
which makes it an ideal candidate for validation of the new modeling library^[Note that CellML does not allow the callbacks which are required for the
extended valve model, so only the simplified model can be compared. The _CirculationModels.jl_
implementation of [@korakianitis2006numerical] includes the extended model as well.].


Model results from _CirculationModels.jl_ model (\autoref{fig:shi-results}) are a perfect match for the CellML model.
The CellML model was run in three versions: (1) imported into _ModelingToolkit.jl_ using _CellMLToolKit.jl_,
(2) Matlab code exported from CellML, (3) Python/SciPy code exported from CellML. 
Speedup against Matlab and Python is 2 and 3 orders of magnitude, respectively (\autoref{tbl:benchmarks}).

![4-chamber, full-circulationmodel from [@korakianitis2006numerical]. Groupings in dashed rectangles are implemented as compound subsystems, which in turn have been composed from individual resistor, compliance, and inertance elements. Ventricles and Atria implemented as time-variable elastances. Both simplified and dynamic, non-linear valve models are implemented. \label{fig-shi-diagram}](fig-shi-diagram.pdf){width=100%}

| CirculationModels.jl | CellMLToolkit.jl | Matlab (ode45) | Python (scipy.solve) |
|:--------------------:|:----------------:|:--------------:|:--------------------:|
| 1x                   | 1.6x             | 272x           | 963x                 |
: Simulation time comparison for a single run of [@korakianitis2006numerical]. CirculationModels.jl model was implemented from scratch. CellML model was imported into ModelingToolkit.jl using CellMLToolkit.jl, Matlab and Python models were created from the CellML code and downloaded from the [CellML Model Repository](http://models.cellml.org/exposure/c49d416ae3a5132882e6ea7479ba50f5/ModelMain.cellml/view).  \label{tbl:benchmarks}

![(a-f) Results for simplified model [@korakianitis2006numerical] implemented in _CirculationModels.jl_. These results match the results from the CellML models (not shown). \label{fig:shi-results}](fig-shi-results.pdf){width=100%}


# Acknowledgements

Harry Saxton is supported by a Sheffield Hallam University Graduate Teaching Assistant PhD scholarship.

# References {-}
