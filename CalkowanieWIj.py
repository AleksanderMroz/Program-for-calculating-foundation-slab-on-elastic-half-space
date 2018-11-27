import unittest
import numpy as np
from sympy import *
import Kernel as KE
import Zemoczkin as ZE
import time

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt








def calkowanieWIJ(X,Y,EP,N,DX,DY,DEP,DP):
    x=Symbol("x")
    y=Symbol("y")
    ep=Symbol("ep")
    n=Symbol("n")
    E0=Symbol("E0")
    v=Symbol("v")

    B=1/((((ep-x)**2+(n-y)**2))**0.5)
    Di=1/(DX*DY)
    Dj=1/(DEP*DP)

    ### CALKOWANIE q(x,y) po dx i dy daje nam R
    C = B.subs([(x, X), (y, Y), (ep, EP), (n,N)])
    print(C)

    CX = integrate(C, (x, 0, DX))
    CX = integrate(CX, (y, 0, DY))
    CX = integrate(CX, (ep, 0, DEP))
    CX = integrate(CX, (n, 0, DP))

    OUT=Di*Dj*CX
    print(X,Y,EP,N)
    return OUT

def TablicaFki(A,H,m,n):
    Fki=np.zeros((n,m))
    Fki[0][0]=SamNaSiebie(A,H)

    for i in range (n):
        for j in range(m):
            if(i==0 and j==0):
                pass
            if(Fki[i][j]!=0):
                pass

            else:
                if(i+j<8):
                    #print();print()
                    #print("OK", i, j)
                    Fki[i][j]=JedenWij(A,H,j*A,i*H)
                else:
                    #print("OK", i, j)
                    Fki[i][j] = DrugiWij(A, H, i * A, j * H)
    return Fki

def Calka(CX,GR,GB):

    return CX*GR*GR*GB*GB

def JedenWij(A,H,x,y):
    X=x
    Y=y
    Sum = 0
    DP = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            #print(X, Y)
            D = TablicaXYEN(A, H, X, Y)
            DP[i][j] = D
            Sum += D
            X = X - A/5
        X = x
        Y -= H/5
    return Sum/25

def DrugiWij(A,H,X,Y):
    Di = 1 / ((A) * (H))
    Dj = 1 / ((A) * (H ))

    DX = X
    DY = Y

    B = 1 / ((DX ** 2 + DY ** 2) ** 0.5)
    CX = Calka(B, A , H )
    return CX*Di*Dj



def TablicaXYEN(A,H,X,Y):
    x = Symbol("x")
    y = Symbol("y")
    ep = Symbol("ep")
    n = Symbol("n")

    D=np.zeros((5,5))



    Di = 1 / ((A/5)*(H/5))
    Dj = 1 / ((A/5)*(H/5))

    DX=X
    DY=Y
    for i in range(5):
        for j in range(5):
            #print((DX,DY,i,j))
            B=1/((DX**2+DY**2)**0.5)
            CX = Calka(B,A/5,H/5)
            D[i][j]=Di * Dj * CX
            DX=DX+A/5
        DX=X
        DY=DY+H/5

        SUM=0

    for i in range(5):
        for j in range(5):
            SUM+=D[i][j]

    return SUM / 25

def SamNaSiebie(A,H):
    c=A
    b=H
    Z=ln(1+sqrt((c/b)**2+1))
    Y=ln((c/b)+sqrt((c/b)**2 +1 ))
    X=ln(b/c)+ (b/c)*Y + Z
    P=2*c/b * X
    P=P/c
    return P/1.1697

def SamNaSiebieAleObok(A,H,x,y):
    c=A
    b=H

    P1=2*ln(b/c)
    P2=ln(2*(x/c)-1)
    P3U=2*x/c+1
    P3D=2*x/c-1
    P3=-2*x/c * ln(P3U/P3D)

    T1=(2*x/b +c/b)
    T2=(2*x/b -c/b)
    T3=sqrt((2*x/b +c/b)**2 +1)
    T4 = sqrt((2 * x / b - c / b) ** 2 + 1)

    P4=b/c *ln((T1+T3)/(T2+T4))

    P5=2*x/c *ln((1+T3)/(1+T4))

    P6= ln((1+T3)*(1+T4))
    P=c/b*(P1-P2-P3+P4+P5+P6)
    return P1,-P2,-P3,P4,P5,P6




def Dob():
    # some 3-dim points

    x = []
    y = []
    z =[]
    X=0
    Y=0
    for p in range(10):
        for j in range(10):
            y.append(Y)
        Y += 1
        for i in range(10):
            x.append(i + X)

    print(x)
    print(y)
    H=0
    for i in range(100):
        z.append((0.3)*H**2)
        H+=1



    data = np.c_[x, y, z]

    # regular grid covering the domain of the data
    X, Y = np.meshgrid(np.arange(-3.0, 3.0, 0.5), np.arange(-3.0, 3.0, 0.5))
    XX = X.flatten()
    YY = Y.flatten()

    order = 1  # 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
        C, _, _, _ = np.linalg.lstsq(A, data[:, 2])  # coefficients

        # evaluate it on grid
        Z = C[0] * X + C[1] * Y + C[2]

        # or expressed using matrix/vector product
        # Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:, :2], np.prod(data[:, :2], axis=1), data[:, :2] ** 2]
        C, _, _, _ = np.linalg.lstsq(A, data[:, 2])

        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX * YY, XX ** 2, YY ** 2], C).reshape(X.shape)

    # plot points and fitted surface
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='r', s=50)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    ax.axis('tight')
    plt.show()


def Dob2(Mac):


    N,M=np.shape(Mac)
    n=N*2-1
    m=M*2-1
    X=np.zeros((n,m))

    for i in range(N):
        for j in range(M):
            X[i*2][j*2]=Mac[i][j]



    for i in range(int((n-1)/2)):
        for j in range(int((m-1)/2)):
            X[2*i+1][2*j+1]=(X[2*i][2*j]+X[2*i][2*j+2]+X[2*i+2][2*j+2]+X[2*i+2][2*j])/4


    for i in range(int((m-1)/2)):
        X[0][2 * i + 1] = (X[0][2 * i] + X[0][2 * i + 2]) / 2
        X[n-1][2 * i + 1] = (X[n-1][2 * i] + X[n-1][2 * i + 2]) / 2

    for i in range(int((n-1)/2)):
        X[2 * i + 1][0] = (X[2 * i][0] + X[2 * i + 2][0]) / 2
        X[2 * i + 1][m-1] = (X[2 * i][m-1] + X[2 * i + 2][m-1]) / 2


    for i in range(n):
        for j in range(m):
            if(X[i][j]==0):
                X[i][j]=(X[i+1][j]+X[i-1][j]+X[i][j+1]+X[i][j-1])/4

    return X




if __name__ == '__main__':
    print(TablicaFki(0.33333,1,5,5))
    print("X")
    print(TablicaFki(1,0.33333,5,5))