"""
@ A program for analytical and numerical solutions for one-dimensional shock tube problem (Riemann problem), and the corresponding OpenFOAM case
@ author ZimoJupiter
@ date 08 Mar 2024
@ license MIT License
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

#%% Initial Parameters
rho1 = 1.0
rho4 = 2.0
p1 = 1e5
p4 = 2e5
gamma = 1.4
t_end = 0.001
R = 287.1

#%% Exact Solution
a1 = np.sqrt(gamma*p1/rho1)
a4 = np.sqrt(gamma*p4/rho4)

def f(p2):
    return ((p2/p1)*(1-((gamma-1)*(a1/a4)*((p2/p1)-1))/ \
                 np.sqrt(2*gamma*(2*gamma+(gamma+1)*((p2/p1)-1))))** \
                    (-2*gamma/(gamma-1))) - p4/p1

p2 = fsolve(f, p1)

g = ((gamma*p2)/(p2-p1))-((gamma-1)/2)
W = np.sqrt(g*(p2-p1)/(rho1))
u2 = W/g
rho2 = g/(g-1)

p3 = p2
u3 = u2
rho3 = rho4*(p3/p4)**(1/gamma)
a3 = np.sqrt(gamma*p3/rho3)

t = t_end
x1 = 1+W*t
x2 = 1+u2*t
x3 = 1-(a3-u3)*t
x4 = 1-a4*t

x43 = np.arange(x4,x3,0.0001)
u43 = (u3*(a4*t+x43-1))/((a4-a3+u3)*t)
rho43 = rho4*(1-((gamma-1)*u43)/(2*a4))**(2/(gamma-1))
p43 = p4*(1-((gamma-1)*u43)/(2*a4))**((2*gamma)/(gamma-1))

x04 = np.arange(0,x4,0.0001)
x04 = x04[:-1]
rho04 = rho4*np.ones(len(x04))
p04 = p4*np.ones(len(x04))
u04 = np.zeros(len(x04))

x32 = np.arange(x3,x2,0.0001)
x32 = x32[:-1]
rho32 = rho3*np.ones(len(x32))
p32 = p3*np.ones(len(x32))
u32 = u3*np.ones(len(x32))

x21 = np.arange(x2,x1,0.0001)
x21 = x21[:-1]
rho21 = rho2*np.ones(len(x21))
p21 = p2*np.ones(len(x21))
u21 = u2*np.ones(len(x21))

x10 = np.arange(x1,2,0.0001)
rho10 = rho1*np.ones(len(x10))
p10 = p1*np.ones(len(x10))
u10 = np.zeros(len(x10))

x = np.concatenate([x04,x43,x32,x21,x10])
rho = np.concatenate([rho04,rho43,rho32,rho21,rho10])
p = np.concatenate([p04,p43,p32,p21,p10])
u = np.concatenate([u04,u43,u32,u21,u10])
T = p/rho/R

plt.figure(figsize=(8, 6), dpi=300)
plt.suptitle('Exact Solution /t=0.001s')
plt.subplot(4, 1, 1)
plt.plot(x, rho)
plt.ylim([0.9, 2.1])
plt.xlabel('x')
plt.ylabel(r'$\rho (kg/m^3)$')
plt.subplot(4, 1, 2)
plt.plot(x, p)
plt.ylim([0.9e5, 2.1e5])
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('x')
plt.ylabel('p (Pa)')
plt.subplot(4, 1, 3)
plt.plot(x, u)
plt.ylim([-5, 100])
plt.xlabel('x')
plt.ylabel('u (m/s)')
plt.tight_layout()
plt.subplot(4, 1, 4)
plt.plot(x, T)
plt.ylim([310, 390])
plt.xlabel('x')
plt.ylabel('T (K)')
plt.tight_layout()
plt.savefig("Exact Solution.png")
# plt.show()
rho_exc = rho.copy()
p_exc = p.copy()
u_exc = u.copy()
x_exc = x.copy()
T_exc = T.copy()

#%% Numerical Solution (Original flux)
# Meshing
x = np.linspace(0, 2, 100)
Delta_x = (x[1]-x[0])
Delta_t = 1e-5
t = np.arange(0, t_end+Delta_t, Delta_t)

# Initialization
p = np.array([p4]*int(len(x)/2)+[p1]*int(len(x)/2))
F_plus = np.array([[0.0]*len(x),[0.0]*len(x),[0.0]*len(x)])
F_minus = np.array([[0.0]*len(x),[0.0]*len(x),[0.0]*len(x)])
rho = np.array([rho4]*int(len(x)/2)+[rho1]*int(len(x)/2))
u = np.array([0.0]*len(x))
rho_u = rho*u
pFplus = np.array([[0.0]*len(x),[0.0]*len(x),[0.0]*len(x)])
pFminus = np.array([[0.0]*len(x),[0.0]*len(x),[0.0]*len(x)])
E = p/(gamma-1)+0.5*u**2*rho

# Calculation
for j in range(len(t)):

   #Steger-Warming FVS
   F = [rho_u, rho*u**2+p, u*(E+p)]
   C = np.sqrt(abs(gamma*p/rho))
   lambdaa = np.array([u,u-C,u+C])
   lambdaplus = (abs(lambdaa)+lambdaa)/2
   lambdaminus = (-abs(lambdaa)+lambdaa)/2
   for i in range(len(x)):
      F_plus[0][i] = rho[i]/(2*gamma)* \
         (2*(gamma-1)*lambdaplus[0][i]+lambdaplus[1][i]+lambdaplus[2][i])
      F_plus[1][i] = rho[i]/(2*gamma)* \
         (2*(gamma-1)*lambdaplus[0][i]*u[i]+lambdaplus[1][i]*(u[i]-C[i])+lambdaplus[2][i]*(u[i]+C[i]))
      F_plus[2][i] = rho[i]/(2*gamma)* \
         ((gamma-1)*lambdaplus[0][i]*u[i]**2+1/2*lambdaplus[1][i]*(u[i]-C[i])**2+ \
          1/2*lambdaplus[2][i]*(u[i]+C[i])**2+(3-gamma)/(2*gamma-2)*(lambdaplus[1][i]+lambdaplus[2][i])*(C[i])**2)
      F_minus[0][i] = rho[i]/(2*gamma)* \
         (2*(gamma-1)*lambdaminus[0][i]+lambdaminus[1][i]+lambdaminus[2][i])
      F_minus[1][i] = rho[i]/(2*gamma)* \
         (2*(gamma-1)*lambdaminus[0][i]*u[i]+lambdaminus[1][i]*(u[i]-C[i])+lambdaminus[2][i]*(u[i]+C[i]))
      F_minus[2][i] = rho[i]/(2*gamma)* \
         ((gamma-1)*lambdaminus[0][i]*u[i]**2+1/2*lambdaminus[1][i]*(u[i]-C[i])**2+ \
          1/2*lambdaminus[2][i]*(u[i]+C[i])**2+(3-gamma)/(2*gamma-2)*(lambdaminus[1][i]+lambdaminus[2][i])*(C[i])**2)
   
   # 5-order WENO
   for i in range(len(x)-4):
      IS0plus = 0
      IS1plus = 0
      IS2plus = 0
      IS0minus = 0
      IS1minus = 0
      IS2minus = 0
      for k in range(3):
            IS0plus += 1/4*((F_plus[k][i]-4*F_plus[k][i+1]+3*F_plus[k][i+2])**2)+ \
               13/12*((F_plus[k][i]-2*F_plus[k][i+1]+F_plus[k][i+2])**2)
            IS1plus += 1/4*((F_plus[k][i+1]-F_plus[k][i+3])**2)+ \
               13/12*((F_plus[k][i+1]-2*F_plus[k][i+2]+F_plus[k][i+3])**2)
            IS2plus += 1/4*((3*F_plus[k][i+2]-4*F_plus[k][i+3]+F_plus[k][i+4])**2)+ \
               13/12*((F_plus[k][i+2]-2*F_plus[k][i+3]+F_plus[k][i+4])**2)
            IS0minus += 1/4*((F_minus[k][i+4]-4*F_minus[k][i+3]+3*F_minus[k][i+2])**2)+ \
               13/12*((F_minus[k][i+4]-2*F_minus[k][i+3]+F_minus[k][i+2])**2)
            IS1minus += 1/4*((F_minus[k][i+3]-F_minus[k][i+1])**2)+ \
               13/12*((F_minus[k][i+3]-2*F_minus[k][i+2]+F_minus[k][i+1])**2)
            IS2minus += 1/4*((3*F_minus[k][i+2]-4*F_minus[k][i+1]+F_minus[k][i])**2)+ \
               13/12*((F_minus[k][i+2]-2*F_minus[k][i+1]+F_minus[k][i])  **2)
            
      omega0plus = (0.1/(IS0plus+1e-6)**2/(0.1/(IS0plus+1e-6)**2+ \
                                           0.6/(IS1plus+1e-6)**2+0.3/(IS2plus+1e-6)**2))
      omega1plus = (0.6/(IS1plus+1e-6)**2/(0.1/(IS0plus+1e-6)**2+ \
                                           0.6/(IS1plus+1e-6)**2+0.3/(IS2plus+1e-6)**2))
      omega2plus = 1-omega1plus-omega0plus
      omega0minus = (0.1/(IS0minus+1e-6)**2/(0.1/(IS0minus+1e-6)**2+ \
                                             0.6/(IS1minus+1e-6)**2+0.3/(IS2minus+1e-6)**2))
      omega1minus = (0.6/(IS1minus+1e-6)**2/(0.1/(IS0minus+1e-6)**2+ \
                                             0.6/(IS1minus+1e-6)**2+0.3/(IS2minus+1e-6)**2))
      omega2minus = 1-omega1minus-omega0minus
      for k in range(3):
            pFplus[k][i+2] = omega0plus*((1/3)*F_plus[k][i]-(7/6)*F_plus[k][i+1]+(11/6)*F_plus[k][i+2])+ \
            omega1plus*(-(1/6)*F_plus[k][i+1] +(5/6)*F_plus[k][i+2]+(1/3)*F_plus[k][i+3])+ \
            omega2plus*((1/3)*F_plus[k][i+2]+(5/6)*F_plus[k][i+3]-(1/6)*F_plus[k][i+4])
            pFminus[k][i+2] = omega0minus*((1/3)*F_minus[k][i+4]-(7/6)*F_minus[k][i+3]+(11/6)*F_minus[k][i+2])+ \
            omega1minus*(-(1/6)*F_minus[k][i+3]+(5/6)*F_minus[k][i+2]+(1/3)*F_minus[k][i+1])+ \
            omega2minus*((1/3)*F_minus[k][i+2]+(5/6)*F_minus[k][i+1]-(1/6)*F_minus[k][i])
   for i in range(len(x)-4):
      if(i>=2):
            rho[i+1] -= (Delta_t/Delta_x)*(pFplus[0][i+1]-pFplus[0][i]+pFminus[0][i+2]-pFminus[0][i+1])
            rho_u[i+1] -= (Delta_t/Delta_x)*(pFplus[1][i+1]-pFplus[1][i]+pFminus[1][i+2]-pFminus[1][i+1])
            E[i+1] -= (Delta_t/Delta_x)*(pFplus[2][i+1]-pFplus[2][i]+pFminus[2][i+2]-pFminus[2][i+1])
   u = rho_u/rho
   p = (gamma-1)*(E-0.5*rho*u**2)
T = p/rho/R

plt.figure(figsize=(8, 6), dpi=300)
plt.suptitle('Numerical Solution /t=0.001s')
plt.subplot(4, 1, 1)
plt.plot(x, rho)
plt.ylim([0.9, 2.1])
plt.xlabel('x')
plt.ylabel(r'$\rho (kg/m^3)$')
plt.subplot(4, 1, 2)
plt.plot(x, p)
plt.ylim([0.9e5, 2.1e5])
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('x')
plt.ylabel('p (Pa)')
plt.subplot(4, 1, 3)
plt.plot(x, u)
plt.ylim([-5, 100])
plt.xlabel('x')
plt.ylabel('u (m/s)')
plt.subplot(4, 1, 4)
plt.plot(x, T)
plt.ylim([310, 390])
plt.xlabel('x')
plt.ylabel('T (K)')
plt.tight_layout()
plt.savefig("Numerical Solution OF.png")
# plt.show()
rho_num_of = rho.copy()
p_num_of = p.copy()
u_num_of = u.copy()
x_num_of = x.copy()
T_num_of = T.copy()

#%% Data of openFoam
p_openfoam = arr = np.genfromtxt("openFoam_p.txt", delimiter=",", dtype=float)
u_openfoam = arr = np.genfromtxt("openFoam_u.txt", delimiter=",", dtype=float)
rho_openfoam = arr = np.genfromtxt("openFoam_rho.txt", delimiter=",", dtype=float)
T_openfoam = arr = np.genfromtxt("openFoam_T.txt", delimiter=",", dtype=float)

#%% Results comparison
plt.figure(figsize=(8, 6), dpi=300)
# plt.suptitle('Solutions /t=0.001s')
plt.subplot(4, 1, 1)
plt.plot(x_exc, rho_exc, color='blue', alpha=0.5, linestyle='--', label='Exact solution', linewidth=2)
plt.plot(x_num_of, rho_num_of, color='red', alpha=0.8, linestyle='-', label='Finite difference method', \
         markersize='5', markevery=10, linewidth=1)
plt.plot(x_num_of, rho_openfoam, color='green', alpha=0.8, linestyle='-', label='openFoam', \
         markersize='5', markevery=10, linewidth=1)
plt.legend(fontsize=8, loc = 'upper right')
plt.ylim([0.9, 2.1])
plt.xlabel('x')
plt.ylabel(r'$\rho (kg/m^3)$')
plt.subplot(4, 1, 2)
plt.plot(x_exc, p_exc, color='blue', alpha=0.5, linestyle='--', label='Exact solution', linewidth=2)
plt.plot(x_num_of, p_num_of, color='red', alpha=0.8, linestyle='-', label='Finite difference method', \
         markersize='3', markevery=10, linewidth=1)
plt.plot(x_num_of, p_openfoam, color='green', alpha=0.8, linestyle='-', label='openFoam', \
         markersize='3', markevery=10, linewidth=1)
plt.ylim([0.9e5, 2.1e5])
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('x')
plt.ylabel('p (Pa)')
plt.subplot(4, 1, 3)
plt.plot(x_exc, u_exc, color='blue', alpha=0.5, linestyle='--', label='Exact solution', linewidth=2)
plt.plot(x_num_of, u_num_of, color='red', alpha=0.8, linestyle='-', label='Finite difference method', \
         markersize='3', markevery=10, linewidth=1)
plt.plot(x_num_of, u_openfoam, color='green', alpha=0.8, linestyle='-', label='openFoam', \
         markersize='3', markevery=10, linewidth=1)
plt.ylim([-5, 100])
plt.xlabel('x')
plt.ylabel('u (m/s)')
plt.subplot(4, 1, 4)
plt.plot(x_exc, T_exc, color='blue', alpha=0.5, linestyle='--', label='Exact solution', linewidth=2)
plt.plot(x_num_of, T_num_of, color='red', alpha=0.8, linestyle='-', label='Finite difference method', \
         markersize='3', markevery=10, linewidth=1)
plt.plot(x_num_of, T_openfoam, color='green', alpha=0.8, linestyle='-', label='openFoam', \
         markersize='3', markevery=10, linewidth=1)
plt.ylim([310, 390])
plt.xlabel('x')
plt.ylabel('T (K)')
plt.tight_layout()
plt.savefig("Results comparison.png")
plt.show()