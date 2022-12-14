import time
start=time.time()
import numpy as np
import math
from scipy import optimize
import matplotlib.pyplot as plt
from itertools import repeat
z=[0.3,0.3,0.4]   #z[i] the fraction of component i in the feed stream of stage 3 with i=0 for pentane , i=1 for n-butane , i=2 for n-pentane
U=[50,0,0,0,0]   #U[j] is the liquid side flow rate outputting stage j
P=100 #feed pressure(psia)
Ph=5171.5   #feed pressure(mmhg)
T=[65,90,115,140,165]  #estimated temperature
F=[0,0,100,0,0]  #F[j]the flow rate of the feed stream to stage j
V = [0, 150, 150, 150, 150,0]  #V[j]the vapor flow rate outputting stage j
a = [15.726, 15.6782, 15.8333]  #a[i]the ANTOINE coefficient for element i
b = [1872.46, 2154.9, 2477.07]  # b[i] the ANTOINE coefficient for i
c = [-25.16, -34.42, -39.94]  # c[i] the ANTOINE coefficient for i
t1=[]
n=36
E=list(repeat([0],n))
X1T=0.1  #XjT :initialize the summation of fractions in stage j
X2T=0.1
X3T=0.1
X4T=0.1
X5T=0.1
O=0
while abs(X1T-1)>10**(-10) and abs(X2T-1)>10**(-10) and abs(X3T-1)>10**(-10) and abs(X4T-1)>10**(-10) and abs(X5T-1)>10**(-10):
    for i in range(len(T)):
        t1.append([0])
        t1[i] = (T[i] - 32) * (5 / 9) + 273.15  # converting T to Kelvin
    P1 = []   #saturation vapor pressure of C3
    P2 = []   #saturation vapor pressure of n-C4
    P3 = []   #P.s of n-C5
    K1 = []   #phase  equilibrium  ratio of (C3)
    K2 = []   #phase equilibrium ratio of (n-C4)
    K3 = []   #phase equilibrium ratio of(n-C5)
    for i in range(5):  # to calculate the saturation vapor pressure by antoine equation
        P1.append([0])
        P1[i] = math.exp(a[0] - (b[0] / (t1[i] + c[0])))
        P2.append([0])
        P2[i] = math.exp(a[1] - (b[1] / (t1[i] + c[1])))
        P3.append([0])
        P3[i] = math.exp(a[2] - (b[2] / (t1[i] + c[2])))
    # calculate phase equilibrium ratio kij of element i of stage j
    for i in range(5):
        K1.append([0])
        K1[i] = (P1[i]) / Ph
        K2.append([0])
        K2[i] = (P2[i]) / Ph
        K3.append([0])
        K3[i] = (P3[i]) / Ph
    # the first matrix
    A1 = [0, 0, 0, 0, 0]
    s = 0
    s1 = 0
    for i in range(1, 5): #calculation of the element A of the matrix(C.2)
        for j in range(i):
            s = s + F[j] - U[j]
        A1[i] = V[i] + s
        s = 0
    B1 = [0, 0, 0, 0, 0]
    for i in range(5): #calculation of the element B of the matrix(C.3)
        for j in range(i + 1):
            s1 = s1 + F[j] - U[j]
        B1[i] = -(V[i + 1] + s1 + U[i] + (V[i] * K1[i]))
        s1 = 0
    C1 = [0, 0, 0, 0, 0]
    for i in range(4): #calculation of the element C of the matrix (C.4)
        C1[i] = V[i + 1] * K1[i + 1]
    D1 = [0, 0, 0, 0, 0]
    for i in range(5): #calculation of D (C.5)
        D1[i] = -F[i] * z[1]
    M1 = [B1[0], C1[0], 0, 0, 0], [A1[1], B1[1], C1[1], 0, 0], [0, A1[2], B1[2], C1[2], 0], [0, 0, A1[3], B1[3],C1[3]], [0, 0, 0, A1[4],B1[4]]
    X1 = np.linalg.solve(M1, D1) #matrix resolution
    # the second matrix
    A2 = [0, 0, 0, 0, 0]
    s2 = 0
    s12 = 0
    for i in range(1, 5):
        A2[0] = 0
        for j in range(i):
            s2 = s2 + F[j] - U[j]
        A2[i] = V[i] + s2
        s2 = 0
    B2 = [0, 0, 0, 0, 0]
    for i in range(5):
        for j in range(i + 1):
            s12 = s12 + F[j] - U[j]
        B2[i] = -(V[i + 1] + s12 + U[i] + (V[i] * K2[i]))
        s12 = 0
    C2 = [0, 0, 0, 0, 0]
    for i in range(4):
        C2[i] = V[i + 1] * K2[i + 1]
    D2 = [0, 0, 0, 0, 0]
    for i in range(5):
        D2[i] = -F[i] * z[1]
    M2 = [B2[0], C2[0], 0, 0, 0], [A2[1], B2[1], C2[1], 0, 0], [0, A2[2], B2[2], C2[2], 0], [0, 0, A2[3], B2[3],C2[3]], [0, 0, 0, A2[4],B2[4]]
    X2 = np.linalg.solve(M2, D2)
    # the third matrix
    A3 = [0, 0, 0, 0, 0]
    s3 = 0
    s13 = 0
    for i in range(1, 5):
        A3[0] = 0
        for j in range(i):
            s3 = s3 + F[j] - U[j]
        A3[i] = V[i] + s3
        s3 = 0
    B3 = [0, 0, 0, 0, 0]
    for i in range(5):
        for j in range(i + 1):
            s13 = s13 + F[j] - U[j]
        B3[i] = -(V[i + 1] + s13 + U[i] + (V[i] * K3[i]))
        s13 = 0
    C3 = [0, 0, 0, 0, 0]
    for i in range(4):
        C3[i] = V[i + 1] * K3[i + 1]
    D3 = [0, 0, 0, 0, 0]
    for i in range(5):
        D3[i] = -F[i] * z[2]
    M3 = [B3[0], C3[0], 0, 0, 0], [A3[1], B3[1], C3[1], 0, 0], [0, A3[2], B3[2], C3[2], 0], [0, 0, A3[3], B3[3], C3[3]], [0, 0, 0, A3[4], B3[4]]
    X3 = np.linalg.solve(M3, D3)
    X1T = X1[0] + X2[0] + X3[0]
    X2T = X1[1] + X2[1] + X3[1]
    X3T = X1[2] + X2[2] + X3[2]
    X4T = X1[3] + X2[3] + X3[3]
    X5T = X1[4] + X2[4] + X3[4]
    # normalization
    X1N = [X1[0] / X1T, X2[0] / X1T, X3[0] / X1T]
    X2N = [X1[1] / X2T, X2[1] / X2T, X3[1] / X2T]
    X3N = [X1[2] / X3T, X2[2] / X3T, X3[2] / X3T]
    X4N = [X1[3] / X4T, X2[3] / X4T, X3[3] / X4T]
    X5N = [X1[4] / X5T, X2[4] / X5T, X3[4] / X5T]
    # Bubble point method to calculate T
    def f1(T1):
        return ((((math.exp(a[0] - (b[0] / (T1 + c[0])))) / Ph) * X1N[0]) + (((math.exp(a[1] - (b[1] / (T1 + c[1])))) / Ph) * X1N[1]) + (((math.exp(a[2] - (b[2] / (T1 + c[2])))) / Ph) * X1N[2])) - 1
    T0 = 291
    res1 = optimize.leastsq(f1, T0)
    T1 = ((float(res1[0] - 273.15)) * (9 / 5) + 32)
    def f2(T2):
        return ((((math.exp(a[0] - (b[0] / (T2 + c[0])))) / Ph) * X2N[0]) + (((math.exp(a[1] - (b[1] / (T2 + c[1])))) / Ph) * X2N[1]) + (((math.exp(a[2] - (b[2] / (T2 + c[2])))) / Ph) * X2N[2])) - 1
    res2 = optimize.leastsq(f2, T0)
    T2 = ((float(res2[0] - 273.15)) * (9 / 5) + 32)
    def f3(T3):
        return ((((math.exp(a[0] - (b[0] / (T3 + c[0])))) / Ph) * X3N[0]) + (((math.exp(a[1] - (b[1] / (T3 + c[1])))) / Ph) * X3N[1]) + (((math.exp(a[2] - (b[2] / (T3 + c[2])))) / Ph) * X3N[2])) - 1
    res3 = optimize.leastsq(f3, T0)
    T3 = ((float(res3[0] - 273.15)) * (9 / 5) + 32)
    def f4(T4):
        return ((((math.exp(a[0] - (b[0] / (T4 + c[0])))) / Ph) * X4N[0]) + ( ((math.exp(a[1] - (b[1] / (T4 + c[1])))) / Ph) * X4N[1]) + (((math.exp(a[2] - (b[2] / (T4 + c[2])))) / Ph) * X4N[2])) - 1
    res4 = optimize.leastsq(f4, T0)
    T4 = ((float(res4[0] - 273.15)) * (9 / 5) + 32)
    def f5(T5):
        return ((((math.exp(a[0] - (b[0] / (T5 + c[0])))) / Ph) * X5N[0]) + (((math.exp(a[1] - (b[1] / (T5 + c[1])))) / Ph) * X5N[1]) + (((math.exp(a[2] - (b[2] / (T5 + c[2])))) / Ph) * X5N[2])) - 1
    res5 = optimize.leastsq(f5, T0)
    T5 = ((float(res5[0] - 273.15)) * (9 / 5) + 32)
    T = [T1, T2, T3, T4, T5]
    O=O+1  #iteration count calculator
    E[O-1] = abs(X1T - 1) #to register the error at each iteration


#results:
print("---------------------the molar fraction of propane ----------")
for i in range(5):
    print("---------------------------------")
    print(" X1",i+1,"|", X1[i])

print("---------------------the molar fraction of n-butane ----------")
for i in range(5):
    print("---------------------------------")
    print("X2",i+1,"|", X2[i])
print("---------------------the molar fraction of  n-pentane ----------")
for i in range(5):
    print("---------------------------------")
    print("X3",i+1,"|", X3[i])
print("----------iteration number  ----------")
print("-------------------",O,"-------------------")
print("----------temperature as function stage  ----------")
for i in range(5):
    print("T",i+1,"=",T[i])
print("----------fraction as function stage ----------")
print("sum of fraction in stage 1:",X1T,"\n"
      "sum of fraction in stage 2:", X2T,"\n"
      "sum of fraction in stage 3:",X3T,"\n"
      "sum of fraction in stage 4:",X4T,"\n"
      "sum of fraction in stage 5:",X5T)
#Plot T as fynction of stage by python and DWSIM
x = np.array([1,2,3,4,5])
t = np.array([T[0], T[1], T[2], T[3], T[4]])
td=[79.9771,117.204,147.461,173.673,194.472]
tDWSIM=np.array(td)
fig, ax = plt.subplots()
ax.plot(t, x,"o-",label="temperature using python")
ax.plot(tDWSIM, x,"v--",label="temperature using DWSIM")
ax.set_title("temperature as function stage")
ax.set_xlabel("temperature")
ax.set_ylabel("stage")
ax.legend()
#Plot x as function of stage by python and DWSIM
x1=np.array([X1[0], X1[1], X1[2], X1[3], X1[4]])
x2=np.array([X2[0], X2[1], X2[2], X2[3], X2[4]])
x3=np.array([X3[0], X3[1], X3[2], X3[3], X3[4]])
x1D=[0.592116,0.251307,0.111844,0.0337883,0.00791146]
x2D=[0.363269,0.539759,0.463989,0.38137,0.236737]
x3D=[0.0446157,0.208933,0.424166,0.584842,0.755351]
x1D=np.array(x1D)
x2D=np.array(x2D)
x3D=np.array(x3D)
fig, ax1 = plt.subplots()
ax1.plot(x1,x,"o-", label="propane fraction using python")
ax1.plot(x2,x,"o-", label="n-butane fraction using python")
ax1.plot(x3,x,"o-", label="n-pentane fraction using python")
ax1.plot(x1D,x,"v--", label="propane fraction using DWSIM")
ax1.plot(x2D,x,"v--", label="n-butane fraction using DWSIM")
ax1.plot(x3D,x,"v--", label="n-pentane fraction using DWSIM")
ax1.legend()
ax1.set_title("fraction molar as function of stage")
ax1.set_xlabel("X")
ax1.set_ylabel("stage")
#Plot k as function of stage by python and DWSIM
k1=np.array([K1[0], K1[1], K1[2], K1[3], K1[4]])
k2=np.array([K2[0], K2[1], K2[2], K2[3], K2[4]])
k3=np.array([K3[0], K3[1], K3[2], K3[3], K3[4]])
k1D=[1.4519,2.35598,3.34312,4.41188,5.41944]
k2D=[0.373171,0.672974,1.02636,1.43004,1.82545]
k3D=[0.104986,0.213526,0.353159,0.522392,0.694954]
k1D=np.array(k1D)
k2D=np.array(k2D)
k3D=np.array(k3D)
fig, ax2 = plt.subplots()
ax2.plot(k1,x,"o-", label="Kpropane using python")
ax2.plot(k2,x,"o-", label="Kn-butane using python")
ax2.plot(k3,x,"o-", label="Kn-pentane using python")
ax2.plot(k1D,x,"*--", label="Kpropane using DWSIM")
ax2.plot(k2D,x,"*--", label="Kn-butane using DWSIM")
ax2.plot(k3D,x,"*--", label="Kn-pentane using DWSIM")
ax2.legend()
ax2.set_title("k as function of stage ")
ax2.set_xlabel("K")
ax2.set_ylabel("etage")
#Plot error as function of iterations
erreur=np.array(E)
iteration=np.arange(1, O+1, 1)
fig, axe = plt.subplots()
axe.plot(iteration,erreur,"o-", label="error")
axe.legend()
axe.set_title("error as function of iteration")
axe.set_xlabel("iteration")
axe.set_ylabel("error")
print("excution time",time.time()-start,"second")
plt.show()

