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


## Examples from the library
### Different motile patterns
![A simple simulation with different motile patterns](https://github.com/mastrof/BacterialMotility/blob/main/scripts/simple.gif)
### Rotational diffusion
![Run-tumble motility with and without rotational diffusion](https://github.com/mastrof/BacterialMotility/blob/main/scripts/rotational-diffusion.gif)
### Chemotaxis (Brumley et al. 2019 PNAS)
![Random walk vs perfect chemotaxis vs noisy chemotaxis](https://github.com/mastrof/BacterialMotility/blob/main/scripts/brumley.gif)
