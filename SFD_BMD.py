# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
plt.figure(figsize=(15,6))
span = 1.84/2
area = 0.42
chord = 0.23
ceqroot = chord * 4 / 3.14159
cellipse = []
crect = []
cschrenk = []
xaxis = []
lx=[]
Lift = []
Liftdist=[]
SF = []
BM = []
Width = []
Iskin = 0
t = 0.12 * chord

ww = 0.95*9.81
wwf = ww/100
Weight = 5.47*9.81
for i in range(101):
    crect.append(chord)
    x = i*span/100
    xaxis.append(x)
    ce = ceqroot*(1-(x/span)**2)**0.5
    cellipse.append(ce)
    cschrenk.append((ce+chord)/2)

plt.plot(xaxis,cschrenk,'r--',label= "Schrenk's Chord")
plt.plot(xaxis,cellipse,'g--', label = "Equivalent Elliptic Wing")
plt.plot(xaxis,crect,'b--', label="Actual Wing")
plt.xlabel("Span(m)")
plt.grid(True)
plt.ylabel("Chord (m)")
plt.legend()
plt.show()
AreaUC = np.trapz(cschrenk,xaxis)
k = Weight / (2*AreaUC) 
for i in cschrenk:
    Lift.append(i*k/100)
    Liftdist.append(i*k)

for x in xaxis:
    lx.append(x)
    lx.append(-x)
lx.sort()
Liftdist.sort()
for i in cschrenk:
    Liftdist.append(i*k)
for i in range(len(Lift)):
    if(i==0):
        SF.append(-Lift[i]+(Weight)/2)
        BM.append(0)
    else:
        SF.append(SF[i-1] - Lift[i] )
        product = [a * b for a,b in zip(Lift[0:i][::-1],xaxis[0:i])]
        BM.append(sum(product))
plt.figure(figsize=(15,6))
plt.xlabel("Span (m)")
plt.ylabel("Shear Force (N)")
plt.grid(True)
plt.plot(xaxis,SF,'r--')
plt.show()
BM = BM[::-1]
plt.figure(figsize=(15,6))
plt.grid(True)
plt.xlabel("Span (m)")
plt.ylabel("Bending Moment (N)")
plt.plot(xaxis,BM,'g--')
plt.show()
max(BM)

# %%
plt.figure(figsize=(15,6))
plt.xlabel("Span (m)")
plt.ylabel("Lift distribution (N/m)")
plt.grid(True)
plt.legend()
plt.plot(lx,Liftdist,'c--', label = 'Lift Distribution')

# %%



print('HII')

sol = 7/36 - 1/54 - 49/(27*12) - 1/(20*27) + 7/(8*27)
print(sol)