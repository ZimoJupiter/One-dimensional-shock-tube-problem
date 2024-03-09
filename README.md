# One-dimensional-shock-tube-problem
A solver for analytical and numerical solutions of one-dimensional shock tube problem (Riemann problem), and the corresponding modified OpenFOAM case.

## Initial conditions
![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Initial%20conditions.png)

## Wave propagation patterns
![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/After%20diaphragm%20broke.png)

## Exact solutions
Conservation of mass, momentum and energy flux is satisfied on the shock wave. The equation for parameters relation on both sides of shock wave, in other words, region 2 and region 1 can be written as following form:
$$\rho_2(u_2 - u_s) = \rho_1(u_1 - u_s)$$
$$\rho_2 u_2(u_2 - u_s) + p_2 = \rho_1 u_1(u_1 - u_s) + p_1$$
$$\rho_2 E_2 (u_2 - u_s) + u_2 p_2 = \rho_1 E_1 (u_1 - u_1) + u_1 p_1$$

Where: $u_s$ is velocity of normal shock wave, $E_1 = \frac{RT_1}{\gamma -1} + \frac{u_1^2}{2}$ and $E_2 = \frac{RT_2}{\gamma -1} + \frac{u_2^2}{2}$ are functions of other parameters, $R=287.1 J/kg\cdot K$ is a constant, $u_s$ is velocity of shock wave. Meanwhile, relations for expansion wave need to satisfy the isentropy relation and  the equality of the Riemann invariants:
$$\displaystyle\frac{p_3}{\rho_3^\gamma} =  \displaystyle\frac{p_4}{\rho_4^\gamma} $$
$$u_3 +  \displaystyle\frac{2a_3}{\gamma-1} = u_4 +  \displaystyle\frac{2a_4}{\gamma-1} $$

Where: $a_3 = \sqrt{\frac{\gamma p_3}{\rho_3}}$, $a_4 = \sqrt{\frac{\gamma p_4}{\rho_4}}$ represent local sound speed in region 3 and region4. The following derivation of the analytical solution can be given by following the partitioning within the shock tube in Fig. \ref{fig:Development of shock tube problem}.Initial condition of region 1 and region 4 is considered as known: $P_1 = 1\times10^5 Pa$, $P_4 = 2\times10^5 Pa$, $\rho_1 = 1 kg/m^3$, $\rho_2 = 2 kg/m^3$, $T_1 = T_4 = 300K$. For air, specific heat capacity ratio $\gamma = 1.4$. The calculation range is $0 \leq x \leq 2$

Firstly, pressure, velocity, density and temperature crossing shock wave in region 2 can be derived as:
$$\displaystyle\frac{p_2}{p_1} =  \displaystyle\frac{p_4}{p_1}\left[ 1 -  \displaystyle\frac{(\gamma - 1)(\frac{a_1}{a_4})(\frac{p_2}{p_1} - 1)}{\sqrt{2\gamma\left(2\gamma + (\gamma + 1)\left(\frac{p_2}{p_1} - 1\right)\right)}} \right]^ {\frac{2\gamma}{\gamma - 1}}$$

Where: $a_1 = \sqrt{\frac{\gamma p_1}{\rho_1}}$ represents local sound speed in region 1.
$$u_2 = \displaystyle\frac{a_1}{\gamma} \left( \displaystyle\frac{p_2}{p_1}-1 \right) \left(\displaystyle\frac{\frac{2\gamma}{\gamma+1}}{\frac{p_2}{p_1} + \frac{\gamma-1}{\gamma+1}}\right)^{0.5}$$
$$\rho_2 = p_2 R T_2$$
$$\displaystyle\frac{T_1}{T_2} = \displaystyle\frac{p_2}{p_1}\left(\displaystyle\frac{\frac{\gamma+1}{\gamma-1}+\frac{p_2}{p_1}}{1+\left(\frac{p_2}{p_1}\right)\left(\frac{\gamma+1}{\gamma-1}\right)} \right)$$

Then the properties in region 2 is obtained. Only density in both sides of contact surface changes suddenly, however, other values do not vary. Density can be calculated by the relation of expansion wave. Parameters in region 3 can be get as following procedures.
$$p_3 = p_2$$
$$u_3 = u_2$$
$$\rho_3 = \rho_4 \left(\displaystyle\frac{p_3}{p_4}\right)^{\frac{1}{\gamma}}$$
$$T_3 = T_2$$

After getting unknown values of region 2 and region 3. The expansion wave can be studied. Movement of surfaces contacted with region 3 and region 4 can be derived as:
$$x_3 = 1-(a_3-u_3)t$$
$$x_4 = 1-a_4 t$$

Where: $a_3 = \sqrt{\frac{\gamma p_3}{\rho_3}}$ is local sound speed in region 3, $x_3$ and $x_4$ are displacement of surfaces contacted with region 3 and region 4, $t$ is time. Then the velocity, pressure, density and temperature distribution within expansion wave can be derived as:
$$u_e = \displaystyle\frac{u_3(a_4 t + x_e -1)}{(a_4 -a_3 +u_3)t}$$
$$\rho_e = \rho_4\left[1-\displaystyle\frac{(\gamma-1)u_e}{2a_4}\right]^\frac{2}{\gamma -1}$$
$$p_e = p_4\left[1-\displaystyle\frac{(\gamma-1)u_e}{2a_4}\right]^\frac{2\gamma}{\gamma -1}$$
$$T_e = \displaystyle\frac{p_e}{\rho_e R}$$

Where: $u_e$, $x_e$, $p_e$, $\rho_e$, $T_e$ is velocity, displacement, pressure and density distribution within in expansion wave, correspondingly. Finally, the displacement of contact surface and shock wave can be written as:
$$x_s = 1+ \sqrt{g\displaystyle\frac{p_2-p_1}{\rho_1}}t$$
$$x_c = 1+u_2 t$$

Where: $x_s$ and $x_c$ corresponding to displacement of shock wave and contact surface. Above is the process of obtaining an analytical solution to the one-dimensional surge tube problem.

## Numerical solutions
This project investigates one-dimensional inviscid flow, so the control equations are chosen as one-dimensional Euler equation:

$$
\left\{\begin{matrix}
 \displaystyle\frac{\partial\rho}{\partial t} + \displaystyle\frac{\partial(\rho u)}{\partial x} = 0 \\
 \displaystyle\frac{\partial(\rho u)}{\partial t} + \displaystyle\frac{\partial(\rho u^2 +p)}{\partial x} = 0 \\
 \displaystyle\frac{\partial(\rho E)}{\partial t} + \displaystyle\frac{\partial(\rho E u + pu)}{\partial x} = 0 
\end{matrix}\right.
$$

Form two new matrices:

$$
\textbf{U} = \left[\begin{matrix}
    \rho \\ 
    \rho u \\ 
    \rho E    
\end{matrix} \right]
$$

$$
\textbf{F} = \left[\begin{matrix}
    \rho u \\ 
    \rho u^2 +p \\ 
    \rho u E + pu    
\end{matrix} \right]
$$

Then the equation can be rewritten as:
$$\frac{\partial \textbf{U}}{\partial t} + \frac{\partial \textbf{F}}{\partial t} = 0$$

Jacobian matrix $\textbf{A} = \displaystyle\frac{\partial\textbf{F}}{\partial\textbf{U}}$ has 3 eigenvalues: $u$, $u+a$ and $u-a$.

Split flux vector with Steger-Warming scheme:

$$
\left\{\begin{matrix}
        \lambda^+_i = \displaystyle\frac{\lambda_i + |\lambda_i|}{2}\\
        \lambda^-_i = \displaystyle\frac{\lambda_i - |\lambda_i|}{2}
    \end{matrix}\right.
$$

$$
\begin{matrix}
        \textbf{F}^+ = \frac{\rho}{2\gamma}\left[\begin{array}{l}
        \lambda^+_1 + 2(\gamma-1)\lambda^+_2 + \lambda^+_3\\
        (u-a)\lambda^+_1 + 2(\gamma-1)u\lambda^+_2 + (u+a)\lambda^+_3\\
        (H-ua)\lambda^+_1 + 2(\gamma-1)u^2\lambda^+_2 + (H+ua)\lambda^+_3\\
        \end{array} \right]
    \end{matrix}
$$

$$
    \begin{matrix}
        \textbf{F}^- = \frac{\rho}{2\gamma}\left[\begin{array}{l}
        \lambda^-_1 + 2(\gamma-1)\lambda^-_2 + \lambda^-_3\\
        (u-a)\lambda^-_1 + 2(\gamma-1)u\lambda^-_2 + (u+a)\lambda^-_3\\
        (H-ua)\lambda^-_1 + 2(\gamma-1)u^2\lambda^-_2 + (H+ua)\lambda^-_3\\
        \end{array} \right]
    \end{matrix}
$$

Where $H = \displaystyle\frac{\gamma RT}{\gamma-1}+\frac{u^2}{2}$. After splitting, Carry out 5-order WENO scheme for space derivative.

Smoothness indicators:

$$ IS_0^+ = \frac{1}{4}(F_{i-2}^{+}-4F_{i-1}^{+}+3F_{i}^{+})^2+\frac{13}{12}(F_{i-2}^{+}-2F_{i-1}^{+}+F_{i}^{+})^2 $$
        
$$ IS_1^+ = \frac{1}{4}(F_{i-1}^+-F_{i+1}^+)^2+\frac{13}{12}(F_{i-1}^+-2F_{i}^++F_{i+1}^+)^2 $$
        
$$ IS_2^+ = \frac{1}{4}(3F_{i}^+-4F_{i+1}^++F_{i+2}^+)^2+\frac{13}{12}(F_{i}^+-2F_{i+1}^++F_{i+2}^+)^2 $$
        
$$ IS_0^- = \frac{1}{4}(F_{i+2}^--4F_{i+1}^-+3F_{i}^-)^2+\frac{13}{12}(F_{i+2}^--2F_{i+1}^-+F_{i}^-)^2 $$

$$ IS_1^- = \frac{1}{4}(F_{i+1}^--F_{i-1}^-)^2+\frac{13}{12}(F_{i+1}^--2F_{i}^-+F_{i-1}^-)^2 $$

$$ IS_2^- = \frac{1}{4}(3F_{i}^--4F_{i-1}^-+F_{i-2}^-)^2+\frac{13}{12}(F_{i}^--2F_{i-1}^-+F_{i-2}^-)^2 $$

Modified weight:

$$
    \alpha_0 = \displaystyle\frac{\gamma_0}{(\epsilon+IS_0)^2},
    \alpha_1 = \displaystyle\frac{\gamma_1}{(\epsilon+IS_1)^2},
    \alpha_2 = \displaystyle\frac{\gamma_2}{(\epsilon+IS_2)^2}
$$

$$
    \omega_0 = \displaystyle\frac{\alpha_0}{\alpha_0+\alpha_1+\alpha_2},
    \omega_1 = \displaystyle\frac{\alpha_0}{\alpha_1+\alpha_1+\alpha_2},
    \omega_2 = \displaystyle\frac{\alpha_0}{\alpha_2+\alpha_1+\alpha_2}
$$

Where: $\epsilon = 10^{-6}; \gamma_0 = 0.1; \gamma_1 = 0.6, \gamma_2 = 0.3$

Then the fluxes can be derived as:

$$ 
f_{i}^+ = \omega_0^+ ( \frac{1}{3} F_{i-2}^+ - \frac{7}{6} F_{i-1}^+ + \frac{11}{6} F_{i}^+ ) + \omega_1^+ ( -\frac{1}{6} F_{i-1}^+ + \frac{5}{6} F_{i}^+ +  \frac{1}{3} F_{i+1}^+ ) + \omega_2^+ ( \frac{1}{3} F_{i}^+ + \frac{5}{6} F_{i+1}^+ - \frac{1}{6} F_{i+2}^+ )
$$

$$
f_{i}^- = \omega_0^- ( \frac{1}{3} F_{i+2}^- - \frac{7}{6} F_{i+1}^- + \frac{11}{6} F_{i}^- ) + \omega_1^- ( -\frac{1}{6} F_{i+1}^- + \frac{5}{6} F_{i}^- + \frac{1}{3} F_{i-1}^- ) + \omega_2^- ( \frac{1}{3} F_{i}^- + \frac{5}{6} F_{i-1}^- - \frac{1}{6} F_{i-2}^- )
$$

Then the flux item can be written as:

$$
    \displaystyle\frac{\partial\textbf{F}}{\partial x} = \displaystyle\frac{f^+-f^-}{\Delta x}
$$

Use first-order difference formula for time derivative:

$$
    \textbf{U}^{n+1}_i = \textbf{U}^{n}_i - \displaystyle\frac{\Delta t}{\Delta x} \textbf{F}^n_i
$$

## OpenFOAM case
blockMesh

![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Mesh.png)

setFields

rhoCentralFoam

## Results
### Results comparision ($t = 0.001s$)
![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Results%20comparison.png)

### Post-process of OpenFOAM case ($t = 0.001s$)
![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/OpenFOAM%20results.png)

### Animation of OpenFOAM case ($t < 0.03s $)

![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Animation_rho.gif)

![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Animation_p.gif)

![image](https://github.com/ZimoJupiter/One-dimensional-shock-tube-problem/blob/main/Pictures%20and%20results/Animation_u.gif)
