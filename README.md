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
Image processing toolbox should be installed.<br>
All the code files are in `Matlab/src` folder:

`Matlab/src`
- main.m: defining initial conditions. The documentation for the code is provided in the files as comments.
- stochastic_generation: defining an object and transfers to different functions

`Matlab/src/utils/SG_VelProf`
- SGVP: codes to generate velocity profiles.

`Matlab/src/utils/SG_VelField`
- SGVF: codes to reorganize the generated profiles into correlated field.

`Matlab/src/utils/SG_VorX`
- SGVorX: codes to add vortex cores.

`Matlab/src/utils/plots`
- Plots: After generating the object, plot the statistics of the generated field.

</details>

## Python package
<!-- brief blurb or link to docs/install/usage -->
<details>
<summary> Installation:
</summary>

1. 🛠️ Installing Poetry

To install [Poetry](https://python-poetry.org/) (Python dependency management and packaging tool), run the following command in your terminal:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

After installation, make sure Poetry is in your `PATH`
- macOS/Linux
```bash
export PATH="$HOME/.local/bin:$PATH
```

Verify installation:

```bash
poetry version
```

Keep venv inside the project (works great with VS Code) poetry 

```bash
config virtualenvs.in-project true
```

2. Installing environment: 

This environment is set with Python 3.13. Change the requires-python = ">=3.13" in pyproject.tmol file if you have other versions on your PC.

run:
```bash
poetry install
```

In case if you want to make environment from scratch, run:(Do not recommended)

```bash
poetry new project_name
```
</details>

<details>
  <summary>Click to see the documentation for Python files</summary>

### Python files
- Explanation:
</details>

## Stochatically Generated Velocity Field
<img src = "Gen.gif">