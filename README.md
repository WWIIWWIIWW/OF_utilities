# OF_utilities

## Some utilities and solvers used by Kai Zhang - collected and modified based on online sources.:

(a) CalcMixtureFraction: Postprocess mixture fraction.

(b) reactingLMFoam - low Mach reactingFoam solver with runTime mixture fraction calculation and some extra outputs.

(c) reactingLMFoam_diff - low Mach reactingFoam solver with diffusion term changed: mixture averaged method and non-Unity Schmidt number methhod. 

(d) RemoveVortons: Post-process to Remove Vortons for old LEMOS BC condition - because paraview fails with vortons in U, new version of LEMOS is ok.

(e) runTimeReconstructPar - runTimeReconstrct cases.

Email: Kai.Zhang.1@city.ac.uk; kaizhang@kth.se


## Citations of lowMach solver:

Duwig, Christophe, Sébastien Ducruix, and Denis Veynante. "Studying the Stabilization Dynamics of Swirling Partially Premixed Flames by Proper Orthogonal Decomposition." Journal of Engineering for Gas Turbines and Power 134 (2012): 101501.

Duwig, Christophe, et al. "Large Eddy Simulations of a piloted lean premix jet flame using finite-rate chemistry." Combustion Theory and Modelling 15.4 (2011): 537-568.
