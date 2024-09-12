# stochKBE
Code for Stochastic Known Boundary Emulations

Accompanying code to go with the corresponding paper. Split into four files:

- `baseFunctions.R`: Functions that are (in principle) independent of model.
- `modelFunctions.R`: Functions that depend in some way on the particular model (for example specific analytic mean/standard deviation calculations).
- `plotting.R`: Helper functions annd definitions to minimise lines of code dedicated to plotting output.
- `script.R`: The main script, from which emulators and relevant plots can be generated.

Only `script.R` need be run: it sources the plotting functions and the model functions, which in turn source the base functions.
