import numpy as np
import matplotlib.pyplot as plt

b1 = 0.114967 # m
b2 = 0.23069 # m
a = 0.032 # m
K1 = 5.34 + 4*((b1/a)**2)
K2 = 5.34 + 4*((b2/a)**2)
nu = 0.33
E = 69e09 #Pa
t = 5e-04 # m
N_cr1 = (K1*((np.pi)**2)*E/(12*(1 - nu**2)))*((t/b1)**2)  # Pa
N_cr2 = (K2*((np.pi)**2)*E/(12*(1 - nu**2)))*((t/b2)**2)  # Pa
print('N_cr1',N_cr1)
print('N_cr2',N_cr2)
M_wing_max = 14.315 #Nm
I_skin = 24352.84e-12 # m^4

y1 = [-0.13, -3.37, -6.52, -9.39, -11.84, -13.68, -14.79, -15.12, 6.05, 8.48, 10.38, 11.73, 12.53, 12.84, 12.72, 12.29]
q1 = [-0.63, -0.62, -0.59, -0.54, -0.46, -0.35, -0.22, 0.00, 0.60, 0.57, 0.52, 0.46, 0.37, 0.28, 0.17, -0.01]
y2 = [-14.64, -13.55, -11.99, -10.08, -7.94, -5.70, -3.49, -1.44, 0.31, 1.67, 2.52, 11.65, 10.80, 9.80, 8.74, 7.67, 6.66, 5.75, 4.95, 4.29, 3.80, 3.49]
q2 = [0.13, 0.26, 0.38, 0.47, 0.54, 0.58, 0.61, 0.62, 0.61, 0.61, 0.61, -0.11, -0.22, -0.31, -0.39, -0.46, -0.51, -0.55, -0.57, -0.59, -0.61, -0.61]
x1 = [1.11, 5.07, 11.80, 21.16, 32.91, 46.76, 62.35, 79.28, 1.72, 6.19, 13.26, 22.77, 34.46, 48.05, 63.23, 79.65]
x2 = [97.04, 115.14, 133.21, 150.80, 167.49, 182.87, 196.55, 208.22, 217.59, 224.44, 228.62, 96.98, 114.86, 132.77, 150.27, 166.93, 182.32, 196.08, 207.86, 217.35, 224.30, 228.55]
y1_up = [-0.13, -3.37, -6.52, -9.39, -11.84, -13.68, -14.79, -15.12]
y2_up = [-14.64, -13.55, -11.99, -10.08, -7.94, -5.70, -3.49, -1.44, 0.31, 1.67, 2.52]
q1_up = [-0.63, -0.62, -0.59, -0.54, -0.46, -0.35, -0.22, 0.00]
q2_up = [0.13, 0.26, 0.38, 0.47, 0.54, 0.58, 0.61, 0.62, 0.61, 0.61, 0.61]
q_up = q1_up + q2_up
x_up = [1.11, 5.07, 11.80, 21.16, 32.91, 46.76, 62.35, 79.28, 97.04, 115.14, 133.21, 150.80, 167.49, 182.87, 196.55, 208.22, 217.59, 224.44, 228.62]
x_down = [1.72, 6.19, 13.26, 22.77, 34.46, 48.05, 63.23, 79.65, 96.98, 114.86, 132.77, 150.27, 166.93, 182.32, 196.08, 207.86, 217.35, 224.30, 228.55]
q_down = [0.60, 0.57, 0.52, 0.46, 0.37, 0.28, 0.17, -0.01, -0.11, -0.22, -0.31, -0.39, -0.46, -0.51, -0.55, -0.57, -0.59, -0.61, -0.61]

F1 = []
F2 = []

for i in range(len(q1)):
    q1[i] = abs(q1[i])*1000
    F1.append(q1[i]/t + 0.001*M_wing_max*y1[i]/I_skin)
print('max_F1',max(F1))

for i in range(len(q2)):
    q2[i] = abs(q2[i])*1000
    F2.append(q2[i]/t + 0.001*M_wing_max*y2[i]/I_skin)
print('max_F2',max(F2))


plt.plot(x_up, q_up)
plt.grid(True)
plt.title('Shear Flow on upper surface')
plt.xlabel('X (mm)')
plt.ylabel('Shear Flow (N/mm)')
plt.show()

plt.plot(x_down, q_down)
plt.grid(True)
plt.title('Shear Flow on lower surface')
plt.xlabel('X (mm)')
plt.ylabel('Shear Flow (N/mm)')
plt.show()