"""
A module containing functions to calculate uncertainty in QSPR models
"""
import hybmdl

dat_Plp = np.loadtxt("../qspr/Kow_Plp.txt")
dat_Kcc = np.loadtxt("../qspr/Kow_Kcc.txt")
dat_Ksc = np.loadtxt("../qspr/Kow_Ksc.txt")


paras0 = np.log( np.array([4.2, 0.31, 0.69]) )
bnds = ((-10, 10), (-10, 10), (-10, 10))