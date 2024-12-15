import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

b = 0.1681 # m
h = 0.1366 # m
t = 5e-04 # m
L = 5.47 * 9.81 # N
I = 2*(b*(t**3)/12 + b*t*((h/2)**2) + (h**3)*t/12) # m^4
M_max_fus = 4 # Nm

s1 = np.linspace(0, b/2, 20)

q1 = []
for i in range(len(s1)):
    q1.append(L*h*s1[i]*t/(2*I))


s2 = np.linspace(0, h, 20)
q2 = []

for i in range(len(s2)):
    q2.append(q1[-1] + L*(h-s2[i])*s2[i]*t/(2*I))

q2 = np.array(q2)

#I = 1
A = 152.1e-06
E = 69e09
G = 26e09
F_qx = L*h*t*(s1**2)/(4*I)
F_qy = L*h*b*s2*t/(4*I) + L*t*(h*(s2**2)/2 - s2**3/3)/(2*I)
F_qx_b = L*h*b*t*s1/(4*I) - L*h*(s1**2)*t/(4*I)
M_0 = 14.315

def top_member(x):
    F_A, V_A, M_A = x
    V = V_A
    H = F_A - F_qx
    H = np.full_like(s1, H)
    M = M_A + V*s1
    M = np.full_like(s1, M)

    dU_MA = np.trapz(M/(E*I), s1) 
    dU_FA = np.trapz(H/(E*A), s1)
    dU_VA = np.trapz(M*s1/(E*I), s1) + np.trapz(np.ones_like(s1), s1)*(V/(G*A))

    return [dU_FA, dU_VA, dU_MA]

def vertical_member(x):
    F_A, V_A, M_A = x
    V = F_A - F_qx[-1]
    H = L/2 - F_qy + V_A
    H = np.full_like(s2, H)
    M = M_A + M_0 - V*s2 + V_A*b/2
    M = np.full_like(s2, M) 

    dU_MA = np.trapz(M/(E*I) , s2)
    dU_FA = np.trapz((-M*s2)/(E * I), s2) + V/(G*A)*np.trapz(np.ones_like(s2),s2)
    dU_VA = np.trapz(-(M*b/2)/(E*I), s2) + np.trapz(H/(E*A), s2)
    return [dU_FA, dU_VA, dU_MA]

def bottom_member(x):
    F_A, V_A, M_A = x
    V = L/2 - F_qy[-1] + V_A
    H = L*h*t/(4*I) * (b**2/4 - s1**2) - F_A
    H = np.full_like(s1, H)
    M = M_A + M_0 + F_qx_b*h + H*h - V*s1 + V_A*b/2
    M = np.full_like(s1, M)
    dU_MA = np.trapz(M/(E*I) , s1)
    dU_FA = np.trapz(-h*M/(E*I), s1) + np.trapz(-H/(E*A), s1)
    dU_VA = np.trapz(M*s1/(E*A), s1) + np.trapz(np.ones_like(s1), s1)*(V/(G*A)) 
    return [dU_FA, dU_VA, dU_MA]

# Initial guess for F_A and M_A
initial_guess = [1e-05, 1e-05, 1e-05]  # Initial guess for [F_A, M_A]

# Solve the system of equations
sol_top = fsolve(top_member, initial_guess) 
sol_vert = fsolve(vertical_member, initial_guess)
sol_bottom = fsolve(bottom_member, initial_guess)

# Sum the results to get the final solution
F_A_sol = sol_top[0] + sol_vert[0] + sol_bottom[0]
V_A_sol = sol_top[1] + sol_vert[1] + sol_bottom[1]
M_A_sol = sol_top[2] + sol_vert[2] + sol_bottom[2]

print("F_A = ", F_A_sol)
print("M_A = ", M_A_sol)
print("V_A = ", V_A_sol)

def M_top(F_A,V_A, M_A):
    V = V_A
    H = F_A - F_qx
    H = np.full_like(s1, H)
    M = M_A + V*s1
    M = np.full_like(s1, M)

    return M

def M_vertical(F_A, V_A, M_A):
    V = F_A - F_qx[-1]
    H = L/2 - F_qy + V_A
    H = np.full_like(s2, H)
    M = M_A + M_0 - V*s2 + V_A*b/2
    M = np.full_like(s2, M) 

    return M

def M_bottom(F_A, V_A, M_A):
    V = L/2 - F_qy[-1] + V_A
    H = L*h*t/(4*I) * (b**2/4 - s1**2) - F_A
    H = np.full_like(s1, H)
    M = M_A + M_0 + F_qx_b*h + H*h - V*s1 + V_A*b/2
    M = np.full_like(s1, M)

    return M

M1 = M_top(F_A_sol, V_A_sol, M_A_sol)
M2 = M_vertical(F_A_sol, V_A_sol, M_A_sol)
M3 = M_bottom(F_A_sol, V_A_sol, M_A_sol)
print(min(M1))

plt.plot(s1, M1, 'b', label = 'Top Member')
plt.plot(s2+s1[-1], M2, 'r', label = 'Vertical Member')
plt.plot(s2[-1]+s1[-1]+s1, M3, 'g', label = 'Bottom Member')
plt.title("Bending Moment Diagram of Bulkhead")
plt.xlabel("Distance along the half section (m)")
plt.ylabel('Bending Moment  (Nm)')
plt.grid(True)
plt.legend()
plt.show()

plt.plot(s1, q1)
plt.title('Variation of shear flow in the horizontal part from the center to the edge')
plt.grid(True)
plt.xlabel('Distance, $s_{1}$')
plt.ylabel('Shear Flow, $q_{1}$')
plt.show()

plt.plot(s2, q2)
plt.title('Variation of shear flow in the Vertical part')
plt.grid(True)
plt.xlabel('Distance, $s_{2}$')
plt.ylabel('Shear Flow, $q_{2}$')
plt.show()

F1_max = M_max_fus*h/(2*I) + q1[-1]/t
F2 = []
for i in range(len(s2)):
    F2.append(q2[i]/t + M_max_fus*(h-s2[i])*s2[i]*t/(2*I))
F2_max = max(F2)
print(F1_max)
print(F2_max)
print(I)