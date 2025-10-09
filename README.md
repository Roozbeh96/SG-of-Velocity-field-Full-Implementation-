# Stochastic Generation of 2D-TBL Velocity Field

Tools to stochastically generate turbulent boundary-layer (TBL) velocity fields—from synthetic profile generation through reorganization and vortex addition. Available for MATLAB and Python. Find details about the model in the following paper:

1. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E

2. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-modal-velocity-field-in-roughwall-turbulence/800D97E60609F4FE4D92847CCDB1C7A0

3. coming soon

4. https://journals.ametsoc.org/view/journals/bams/105/1/BAMS-D-23-0066.1.xml

## Contents
- [MATLAB package](#matlab-package)
- [Python package](#python-package)
- [Sample of stochatically generated velocity field](#Stochatically-generated-velocity-field)

## MATLAB package (R2025b)
<!-- brief blurb or link to docs/install/usage -->
<!-- e.g., Installation, Quickstart, API, Examples -->

<details>
  <summary>Click to see the documentation for MATLAB package</summary>

### MATLAB files
All the code files are in `Matlab/src` folder:
`Matlab/src`
- main.m: defining initial conditions. The documentation for the code is provided in the files as comments.
-stochastic_generation: defining an object and transfers to different functions
`Matlab/src/SG_VelProf`
- SGVP: codes to generate velocity profiles
`Matlab/src/SG_VelField`
- SGVF: codes to reorganize the generated profiles into correlated field
`Matlab/src/SG_VorX`
- SGVorX: codes to add vortex cores.


</details>

## Python package
<!-- brief blurb or link to docs/install/usage -->
<details>
  <summary>Click to see the documentation for Python files</summary>

### Python files
Full docs here...
- Installation
- Quickstart
- API
- Examples

</details>

## Stochatically Generated Velocity Field
<img src = "Gen.gif">