import cantera as ct
import numpy as np
gas = ct.Solution('BFER_mech.cti')
# setup parameters
Lx=0.02
tol_ss      = [1.0e-6, 1.0e-14]        # [rtol atol] for steady-state problem
tol_ts      = [1.0e-5, 1.0e-13]        # [rtol atol] for time stepping
loglevel    = 0                        # amount of diagnostic output (0
refine_grid = True                     # True to enable refinement
dt          = 1.0e-5


#H2O = 0.57
#composition in {fuel_dry_pct_ini = 1-H2O_pct_ini)
#H2 = 0.083947
CO = 0.157781
CO2 = 0.201272
CH4 = 0.183947
N2 = 0.457

phi = 1.0
fuel = 'CO:%f, CO2:%f, CH4:%f, N2:%f'  %( CO, CO2, CH4, N2)
gas.TP = 300, 1*ct.one_atm
gas.set_equivalence_ratio(phi, fuel , 'O2:0.21, N2:0.79')
print (gas())
f = ct.FreeFlame(gas, width=Lx)
f.transport_model = 'Multi'
f.soret_enabled=True

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)
f.set_refine_criteria(ratio=3, slope=0.01, curve=0.01)

f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)

#print(gas.multi_diff_coeffs)
print (gas.binary_diff_coeffs)
#f.write_csv("a.csv", species='X', quiet=False)
#np.savetxt("a.txt", np.transpose(f.X))
#print (gas())
#gas.TPX = 780, 1*ct.one_atm, f.X[-1,:]
#print(gas.multi_diff_coeffs)
