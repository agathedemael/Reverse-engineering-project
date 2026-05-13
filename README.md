# blade_FEM:
Simple FEM analysis of a rotating blender for an immersion blender. Many parameters have been simplified or neglected to compute an initial blade geometry optimization. The conclusion of this analysis is a trend for optimization rather than a specific shape that performs best in real-world blending scenarios.

# blade_guard_CFD:
Simulation using a simplified incompressible CFD-style solver based on semi-Lagrangian advection, viscous diffusion, and pressure projection. The numerical structure is inspired by standard projection methods for incompressible flow and stable fluids. The blender blade and exchange openings are represented using simplified forcing regions rather than fully resolved rotating geometry. Results must be interpreted qualitatively to compare flow patterns, suction tendency, and relative mixing behavior rather than as verified quantitative predictions.
