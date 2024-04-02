# Converging Radial Flow to a Borehole in a Graded Rock
This repository contains both python (using mpmath) and fortran programs for evaluating the flow solution to a borehole in a gradent material. Both the permeability and the porosity are graded (increasing towards the excavation, as would be caused by damage). The solution includes single-porosity, a Warren-Root style double porosity (i.e., a mass-function coefficient thin-film type relationship between the fracture and matrix), and a Kazemi style double porosity matrix diffusion solution.

The python-based solution is slower to execute but easier to use and modify, and is based on the mpmath extende-precition library. The fortran version runs faster (i.e., for parameter estimation), but it limited to double precision variables.

This solution has been published in Mathematical Geosciences as "Generalized Solution for Double-Porosity Flow through a Graded Excavation Damaged Zone" (April 2024).

The author's publications and related software can be found at https://kris.kuhlmans.net
