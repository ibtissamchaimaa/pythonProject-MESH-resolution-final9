import numpy as np
import math
from scipy import optimize
import matplotlib.pyplot as plt
from itertools import repeat
z=[0.3,0.3,0.4]   #La fraction de feed
U=[50,0,0,0,0]   #Le flux de liquide extrait à  chaque étage
P=100 #La pression de feed en psia
Ph=5171.5   #La pression de feed en mmhg
T=[65,90,115,140,165]  #la temperature estimee
F=[0,0,100,0,0]  #Le feed à chaque étage
V = [0, 150, 150, 150, 150,0]  #Le flux de vapeur
a = [15.726, 15.6782, 15.8333]  # coefficient d'ANTOINE A pour les 3 constituants
b = [1872.46, 2154.9, 2477.07]  # coefficient d'ANTOINE B pour les 3 constituants
c = [-25.16, -34.42, -39.94]  # coefficient d'ANTOINE C pour les 3 constituants
t1=[]
X1T=0.1  #La somme de la composition à l'étage 1 estime
X2T=0.1  #La somme de la composition à l'étage 2 estime
X3T=0.1  #La somme de la composition à l'étage 3 estime
X4T=0.1  #La somme de la composition à l'étage 4 estime
X5T=0.1  #La somme de la composition à l'étage 5 estime
O=0
n=37
E=list(repeat([0],n))
while abs(X1T-1)>10**(-10) and abs(X2T-1)>10**(-10) and abs(X3T-1)>10**(-10) and abs(X4T-1)>10**(-10) and abs(X5T-1)>10**(-10):
    for i in range(len(T)):
        t1.append([0])
        t1[i] = (T[i] - 32) * (5 / 9) + 273.15  # convertir la temperature de F en Kelvin
    P1 = []   #Pression de vapeur saturante de C3
    P2 = []   #Pression de vapeur saturante de n-C4
    P3 = []   #Pression de vapeur saturante de n-C5
    K1 = []   # volatilite de element 1 (C3)
    K2 = []   # volatilite de element 2 (n-C4)
    K3 = []   # volatilite de element 3 (n-C5)
    for i in range(5):  # pour calculer la pression de vapeur saturante Pij de l'element i dans etage j
        P1.append([0])
        P1[i] = math.exp(a[0] - (b[0] / (t1[i] + c[0])))
        P2.append([0])
        P2[i] = math.exp(a[1] - (b[1] / (t1[i] + c[1])))
        P3.append([0])
        P3[i] = math.exp(a[2] - (b[2] / (t1[i] + c[2])))
    # calculer Kij valeur de la volatilite absolue de l'element i dans l'etage j
    for i in range(5):
        K1.append([0])
        K1[i] = (P1[i]) / Ph  # la loi de raout/dalton
        K2.append([0])
        K2[i] = (P2[i]) / Ph
        K3.append([0])
        K3[i] = (P3[i]) / Ph
    # la premiere matrice
    A1 = [0, 0, 0, 0, 0]
    s = 0
    s1 = 0
    for i in range(1, 5): #Pour calculer L'élèment de matrice A
        for j in range(i):
            s = s + F[j] - U[j]
        A1[i] = V[i] + s
        s = 0
    B1 = [0, 0, 0, 0, 0]
    for i in range(5): #Pour calculer l'élèment de matrice B
        for j in range(i + 1):
            s1 = s1 + F[j] - U[j]
        B1[i] = -(V[i + 1] + s1 + U[i] + (V[i] * K1[i]))
        s1 = 0
    C1 = [0, 0, 0, 0, 0]
    for i in range(4): #Pour calculer l'élèment de matrice C
        C1[i] = V[i + 1] * K1[i + 1]
    D1 = [0, 0, 0, 0, 0]
    for i in range(5): #Pour calculer l'élèment de matrice D
        D1[i] = -F[i] * z[1]
    M1 = [B1[0], C1[0], 0, 0, 0], [A1[1], B1[1], C1[1], 0, 0], [0, A1[2], B1[2], C1[2], 0], [0, 0, A1[3], B1[3],C1[3]], [0, 0, 0, A1[4],B1[4]]
    X1 = np.linalg.solve(M1, D1) #La resolution matricielle par la methode de Gauss
    # la deuxieme matrice
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
    # la troisieme matrice
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
    # la normalisation
    X1N = [X1[0] / X1T, X2[0] / X1T, X3[0] / X1T]
    X2N = [X1[1] / X2T, X2[1] / X2T, X3[1] / X2T]
    X3N = [X1[2] / X3T, X2[2] / X3T, X3[2] / X3T]
    X4N = [X1[3] / X4T, X2[3] / X4T, X3[3] / X4T]
    X5N = [X1[4] / X5T, X2[4] / X5T, X3[4] / X5T]
    # estimation de T
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
    E[O] = abs(X1T - 1)
    O=O+1  #Pour calculer le nombre itération

#Les problèmes des moindres carrés, minimisant la norme d'une fonction vectorielle, ont une structure spécifique qui peut être utilisée dans l'algorithme de Levenberg – Marquardt implémenté dans scipy.optimize.leastsq().
print("---------------------la fraction molaire de propane ----------")
for i in range(5):
    print("---------------------------------")
    print(" X1",i+1,"|", X1[i])

print("---------------------la fraction molaire de n-butane ----------")
for i in range(5):
    print("---------------------------------")
    print("X2",i+1,"|", X2[i])
print("---------------------la fraction molaire de n-pentane ----------")
for i in range(5):
    print("---------------------------------")
    print("X3",i+1,"|", X3[i])
print("----------le nombre d'eteration  ----------")
print("-------------------",O,"-------------------")
print("----------la temperature dans chaque etage est  ----------")
for i in range(5):
    print("T",i+1,"=",T[i])
print("----------la somme des x dans chaque etage ----------")
print("la somme de la composition dans l'etage 1 est:",X1T,"\n"
      ,"la somme de la composition dans l'etage 2 est:", X2T,"\n"
      "la somme de la composition dans l'etage 3 est:",X3T,"\n"
      "la somme de la composition dans l'etage 4 est:",X4T,"\n"
      "la somme de la composition dans l'etage 5 est:",X5T)
x = np.array([1,2,3,4,5])
y = np.array([T[0], T[1], T[2], T[3], T[4]])
fig, ax = plt.subplots()
ax.plot(y, x,"o-")
ax.set_title("LA TEMPERATURE DANS CHAQUE ETAGE")
ax.set_xlabel("temperature")
ax.set_ylabel("etage")

y2 = np.array([1,2,3,4,5])
x1=np.array([X1[0], X1[1], X1[2], X1[3], X1[4]])
x2=np.array([X2[0], X2[1], X2[2], X2[3], X2[4]])
x3=np.array([X3[0], X3[1], X3[2], X3[3], X3[4]])
fig, ax1 = plt.subplots()
ax1.plot(x1,y2,"o-", label="propane")
ax1.plot(x2,y2,"o-", label="n-butane")
ax1.plot(x3,y2,"o-", label="n-pentane")
ax1.legend()
ax1.set_title("LA COMPOSITION EN FONCTION DE L ETAGE")
ax1.set_xlabel("X")
ax1.set_ylabel("etage")

 # affiche la figure a l'ecran
e2 = np.array([1,2,3,4,5])
y4=np.array([K1[0], K1[1], K1[2], K1[3], K1[4]])
y5=np.array([K2[0], K2[1], K2[2], K2[3], K2[4]])
y6=np.array([K3[0], K3[1], K3[2], K3[3], K3[4]])
fig, ax2 = plt.subplots()
ax2.plot(y4,e2,"o-", label="Kpropane")
ax2.plot(y5,e2,"o-", label="Kn-butane")
ax2.plot(y6,e2,"o-", label="Kn-pentane")
ax2.legend()
ax2.set_title("LA VOLATILITE ABSOLUE EN FONCTION DE L ETAGE")
ax2.set_xlabel("K")
ax2.set_ylabel("etage")
plt.show()
