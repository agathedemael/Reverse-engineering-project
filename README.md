# blade_FEM:
Simple FEM analysis of a rotating blender for an immersion blender. Many parameters have been simplified or neglected to compute an initial blade geometry optimization. The conclusion of this analysis is a trend for optimization rather than a specific shape that performs best in real-world blending scenarios.

# blade_guard_CFD:
Simulation using a simplified incompressible CFD-style solver based on semi-Lagrangian advection, viscous diffusion, and pressure projection. The numerical structure is inspired by standard projection methods for incompressible flow and stable fluids. The blender blade and exchange openings are represented using simplified forcing regions rather than fully resolved rotating geometry. Results must be interpreted qualitatively to compare flow patterns, suction tendency, and relative mixing behavior rather than as verified quantitative predictions.
# Warning: 
If the stl files of the three different geometries are not downloaded in the same folder, it is necessary to first run blade_guard_geometry.m, then blade_guard_CFD. Both the CFD and the optimisation algorithms are heavy to run on MATLAB, except around 5 minutes and 15 minutes of run time respectively.
