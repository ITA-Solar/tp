#-------------------------------------------------------------------------------
# Created 19.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Constants.jl
#
#-------------------------------------------------------------------------------
# Module containing physical constants
#-------------------------------------------------------------------------------


module Constants

k_B      = 1.380649e-23 # J/K..................... Boltzmann's constant
e        = 1.602176634e-19 # C....................... Elementary charge
c        = 299792458 # m/s.................... Speed of light in vacuum
G        = 6.67430e−11 # m^3 kg^-1 s^-2........ Constant of gravitation
h        = 6.62607015e−34 # J Hz^-1.................. Planck's constant
μ_0      = 1.25663706212e−6 # N A^-2.......Vacuum magnetic permeability
ϵ_0      = 8.8541878128e−12 # F m^-1.......Vacuum electric permittivity
m_e      = 9.1093837015e-31 # Kg ........................ Electron mass
m_p      = 1.67262192369e-27 # Kg ......................... Proton mass
cSqrdInv = 1/c^2 # s^2 m^-2 ........ Inverse of the light speed squared

# Scaling factors
cgs2SI_u = 0.01 #...........................cgs-velocity to SI-velocity
cgs2SI_b = 1e-4 #.....................cgs-magnetic field to SI-velocity


end # module Constants
