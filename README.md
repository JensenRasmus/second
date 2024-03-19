# secondorder_spst

## About

These files contain the files, which are either to replace standard files in the Manopt package for Matlab, or used by ```numerical_examples.m```
 
## List of files
- ```numerical_experiments.m```: Run the experiments by uncommenting selected parts of the code. One defines a problem structure similar to how it is done in Manopt [1], and data is subsequently captured. 
- ```optim_sd.m``` and ```optim_cg.m```: The steepest descend- and conjugate gradient method.
- ```trustregions.m```: Replace this with the corresponding file Manopt to enable the use of  ```symplecticStiefelfactory.m```.
- ```symplecticStiefelfactory.m```: Manifold structure for the Manopt--package. Note that this is not directly compatible with the standard Manopt--package, but this can be done editing a few lines in the solvers one wants to use. One can check ```stiefelFactory.m``` for a comparison. 
- ```williasondiag.m```: Computes the symplectic eigenvalue decomposition. This file is part of the work [2] 

## References
[1]: Manopt, a Matlab Toolbox for Optimization on Manifolds, Nicolas Boumal et.al., Journal of Machine Learning Research, 2014. [Website](https://www.manopt.org)
[2]: Riemannian Optimization on the Symplectic Stiefel Manifold, Bin Gao, Nguyen Thanh Son, P.-A. Absil,Tatjana Stykel, SIAM Journal on Matrix Analysis and Applications, 2021. Source code is found [Here](https://github.com/opt-gaobin/spopt)

## Copyright

Copyright (C) 2024, Rasmus Jensen and Ralf Zimmermann 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/