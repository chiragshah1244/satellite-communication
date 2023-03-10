import numpy as np
import math
import pandas as pd
import sys


#initilazition 
mu = 3.986005e14;                  # WGS 84 value of the earth’s gravitational constant for GPS user
omega_dot_earth = 7.2921151467e-5; # omega.e WGS 84 value of the earth’s rotation rate

#Semi-major axis
a=0.515375480270E+04

#beginning of end of week crossover
t=383760.0
toe=388800.0              

#mean motion correction
dn=0.583845748090E-08 

#Mean anomaly
Mo=-0.286954703389e+01

#eccentric anomal
e=0.167867515702E-01

#Argument of Latitude
w=np.pi

#Argument of Latitude Correction
Cus= 0.277347862720E-05 
Cuc= -0.379979610443E-06 

# Corrected Radius
Crs= -0.965625000000E+01
Crc= 0.293218750000E+03

#Inclination Correction
Cis= 0.173225998878E-06
Cic= 0.199303030968E-06
idot=0.789318592573E-10
io=0.903782727230E+00

# Corrected longitude of ascending node
omg_o= -0.657960408566E+00
omg_dot= -0.868929051526E-08
omg_dot_e= 7.2921151467E-5



A= a**2                        # Semi-major axis
print('A=' ,A)

no=np.sqrt(mu/A**3)               # Computed mean motion (rad/sec)
print('no=',no)

tk=t-toe                          # Time from ephemeris reference epoch
print('tk=',tk)

# account for beginning of end of week crossover
if (tk > 302400):
	tk = tk-604800;
	print('sub 604800',tk)
	
elif (tk < -302400):
	tk = tk+604800;
	print('add 604800',tk)

else :
	tk=tk
	print('no neeed 604800')
	
# apply mean motion correction
n=no+dn;
print('n=',n)
	  
# Mean anomaly
Mk=Mo+(n*tk)
print('Mk=',Mk)

# solve for eccentric anomal
E_k=Mk;

for j in range(0,3):
  E_k=E_k+(Mk - E_k + e*np.sin(E_k))/(1-e*np.cos(E_k))
print('3rd value (E3) E_k=',E_k)

# True anomaly (unambiguous quadrant)
V_k=2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E_k/2))
print('V_k=',V_k)

# Argument of Latitude
ph_k=V_k+w
print('ph_k=',ph_k)

# Argument of Latitude Correction
del_uk=Cus*np.sin(2*ph_k) + Cuc*np.cos(2*ph_k)
print(del_uk)

# Corrected Argument of Latitude
uk= ph_k + del_uk
print('Corrected Latitude uk=',uk)
	 
# Radius Correction
del_rk=Crs*np.sin(2*ph_k)+Crc*np.cos(2*ph_k)
print(del_rk)

#corrected radius
rk = A*(1-e*np.cos(E_k))+ del_rk
print('Corrected radius rk=',rk)

# Inclination Correction
del_ik=Cis*np.sin(2*ph_k) 
print(del_ik)

# Corrected Inclination
ik=io+del_ik+idot*tk
print('Corrected radius ik=',ik)

# Positions in Orbital Plane
X_k= rk*np.cos(uk)
print('X_k=',X_k)
Y_k= rk*np.sin(uk)
print('Y_k=',Y_k)

print('Positions in Orbital Plane',X_k,Y_k)

# Corrected longitude of ascending node
omg_k=omg_o + (omg_dot-omg_dot_e)*tk - omg_dot_e*toe
print ('Corrected longitude of ascending node', omg_k)

# Earth-fixed coordinates
x_k= X_k*np.cos(omg_k)- Y_k*np.cos(ik)*np.sin(omg_k)
print('x_k=',x_k)
y_k= X_k*np.sin(omg_k)+ Y_k*np.cos(ik)*np.cos(omg_k)
print('y_k=',y_k)
z_k= Y_k*np.sin(ik)
print('z_k=',z_k)

pos=[]
pos.append(x_k)
pos.append(y_k)
pos.append(z_k)

print(pos)





