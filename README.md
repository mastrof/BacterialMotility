# BacterialMotility

BacterialMotility is a framework for numerical simulations of individual bacterial behavior (motility, chemokinesis, chemotaxis...) in arbitrary external fields.
The library implements some standard models with highly customizable traits and provides a friendly interface for the integration of user-defined models.

The library currently includes:
- different motile patterns (all with customizable reorientation statistics via `Distributions`)
    * run-and-tumble
    * run-reverse-flick
- agent-based models of chemotaxis
    * Brumley et al. 2019 PNAS
- common external fields
    * steady-state concentration field from a leaking spherical source


![A simple simulation with different motile patterns](https://github.com/mastrof/BacterialMotility/blob/main/scripts/simple.gif)
