# OF_utilities

## Some utilities and solvers used by Kai Zhang - collected and modified based on online sources.:

(a) CalcMixtureFraction: Postprocess mixture fraction.

(b) reactingLMFoam_Z - low Mach reactingFoam solver with runTime mixture fraction calculation and some extra outputs.

(c) reactingLMFoam_diff - low Mach reactingFoam solver with diffusion term changed: mixture averaged method and non-Unity Schmidt number methhod. 

(d) RemoveVortons: Post-process to Remove Vortons for old LEMOS BC condition - because paraview fails with vortons in U, new version of LEMOS is ok.

(e) runTimeReconstructPar - runTimeReconstrct cases.

(f) reactingLMFoam_TFM - low Mach reactingFoam solver with Thickend Flame Model implemented.

Email: Kai.Zhang.1@city.ac.uk; kaizhang@kth.se


## We kindly ask you to cite the following paper if you are using lowMach solver:

Zhang, K., Shen, Y.Z., Duwig, C., 2021. Finite rate simulations and analyses of wet/distributed flame structure in swirl-stabilized combustion. Fuel, 289, pp.119922

Zhang, K., Dybe, S., Shen, Y., Schimek, S., Paschereit, C.O. and Duwig, C., 2020. Experimental and Numerical Investigation of Ultra-Wet Methane Combustion Technique for Power Generation. Journal of Engineering for Gas Turbines and Power.

Duwig, C., Ducruix, S. and Veynante, D., 2012. Studying the stabilization dynamics of swirling partially premixed flames by proper orthogonal decomposition. Journal of Engineering for Gas Turbines and Power, 134(10).

Duwig, C., Nogenmyr, K.J., Chan, C.K. and Dunn, M.J., 2011. Large eddy simulations of a piloted lean premix jet flame using finite-rate chemistry. Combustion Theory and Modelling, 15(4), pp.537-568.
