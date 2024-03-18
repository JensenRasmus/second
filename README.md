# secondorder_spst

# About

These files contain the files, which are either to replace standard files in the Manopt package for Matlab, or used by ```numerical_examples.m```

# List of files
- ```numerical_experiments.m```: Run the experiments by uncommenting selected parts of the code
- ```optim_sd.m``` and ```optim_cg.m```: The steepest descend and conjugate gradient method
- ```trustregions.m```: Replace this with the corresponding file in the manopt package to enable the use of  ```symplecticStiefelfactory.m```
- ```symplecticStiefelfactory.m```: Manifold structure for the Manopt--package. Note that this is not directly compatible with the standard Manopt--package, but this can be done editing a few lines in the solvers one wants to use. One can check ```stiefelFactory.m``` for a comparison. 