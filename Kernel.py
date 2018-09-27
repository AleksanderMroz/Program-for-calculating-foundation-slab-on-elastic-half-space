
import numpy as np
from sympy import *
import time as time
init_printing(use_unicode=True)
import sys
import Zemoczkin as ZE
import AproxFunkKsztal as AFK

import CalkowanieWIj as CAL
def fK0():
    ### Macierz Bezwładnosci ###
    B0 = [["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"],
          ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"],
          ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]]

    x = Symbol("x")
    y = Symbol("y")
    b = Symbol("b")
    a = Symbol("a")


    B0[0][0] = (6 * (2 * x - 1) * (y - 1)) / (a ** 2)
    B0[0][1] = (0) / (a ** 2)
    B0[0][2] = (-2 * a * (3 * x - 2) * (y - 1)) / (a ** 2)
    B0[0][3] = (-6 * (2 * x - 1) * (y - 1)) / (a ** 2)
    B0[0][4] = (0) / (a ** 2)
    B0[0][5] = (-2 * a * (3 * x - 1) * (y - 1)) / (a ** 2)
    B0[0][6] = (6 * y * (2 * x - 1)) / (a ** 2)
    B0[0][7] = (0) / (a ** 2)
    B0[0][8] = (2 * a * y * (3 * x - 1)) / (a ** 2)
    B0[0][9] = (6 * y * (-2 * x + 1)) / (a ** 2)
    B0[0][10] = (0) / (a ** 2)
    B0[0][11] = (2 * a * y * (3 * x - 2)) / (a ** 2)

    B0[1][0] = (6 * (x - 1) * (2 * y - 1)) / (b ** 2)
    B0[1][1] = (2 * b * (x - 1) * (3 * y - 2)) / (b ** 2)
    B0[1][2] = (0) / (b ** 2)
    B0[1][3] = (6 * x * (-2 * y + 1)) / (b ** 2)
    B0[1][4] = (2 * b * x * (-3 * y + 2)) / (b ** 2)
    B0[1][5] = (0) / (b ** 2)
    B0[1][6] = (6 * x * (2 * y - 1)) / (b ** 2)
    B0[1][7] = (2 * b * x * (-3 * y + 1)) / (b ** 2)
    B0[1][8] = (0) / (b ** 2)
    B0[1][9] = (-6 * (x - 1) * (2 * y - 1)) / (b ** 2)
    B0[1][10] = (2 * b * (x - 1) * (3 * y - 1)) / (b ** 2)
    B0[1][11] = (0) / (b ** 2)

    B0[2][0] = (12 * x ** 2 - 12 * x + 12 * y ** 2 - 12 * y + 2) / a / b
    B0[2][1] = (2 * b * (y - 1) * (3 * y - 1)) / a / b
    B0[2][2] = (-2 * a * (x - 1) * (3 * x - 1)) / a / b
    B0[2][3] = (-12 * x ** 2 + 12 * x - 12 * y ** 2 + 12 * y - 2) / a / b
    B0[2][4] = (-2 * b * (y - 1) * (3 * y - 1)) / a / b
    B0[2][5] = (2 * a * x * (-3 * x + 2)) / a / b
    B0[2][6] = (12 * x ** 2 - 12 * x + 12 * y ** 2 - 12 * y + 2) / a / b
    B0[2][7] = (2 * b * y * (-3 * y + 2)) / a / b
    B0[2][8] = (2 * a * x * (3 * x - 2)) / a / b
    B0[2][9] = (-12 * x ** 2 + 12 * x - 12 * y ** 2 + 12 * y - 2) / a / b
    B0[2][10] = (2 * b * y * (3 * y - 2)) / a / b
    B0[2][11] = (2 * a * (x - 1) * (3 * x - 1)) / a / b

    ### Transponowanie macierzy B ###
    X = Matrix(B0)
    Y = X.T

    ### MACIERZ ODKSZTALCEN ###
    v = Symbol("v")
    E0 = Symbol("E0")
    h0 = Symbol("h0")

    D1 = [["0", "0", "0"], ["0", "0", "0"], ["0", "0", "0"]]
    # D1 = [[1, v, 0], [v, 1, 0], [0, 0, ((1 - v) / 2)]]

    D1[0][0] = E0 * h0 ** 3 / 12 / (1 - v ** 2)
    D1[0][1] = v * E0 * h0 ** 3 / 12 / (1 - v ** 2)
    D1[0][2] = 0

    D1[1][0] = v * E0 * h0 ** 3 / 12 / (1 - v ** 2)
    D1[1][1] = E0 * h0 ** 3 / 12 / (1 - v ** 2)
    D1[1][2] = 0

    D1[2][0] = 0
    D1[2][1] = 0
    D1[2][2] = ((1 - v) / 2) * E0 * h0 ** 3 / 12 / (1 - v ** 2)

    Z = Matrix(D1)

    #M=Macierz K0e niepocałkowana i pomnożona przaz a i b
    M = Y * Z * X
    return M




### FUNKCJA CALKOWANIA SYMBOLICZNEGO ###

def SymIntegrate(M,A,B,V,e0,H0):
    x = Symbol("x")
    y = Symbol("y")
    b = Symbol("b")
    a = Symbol("a")
    v = Symbol("v")
    E0 = Symbol("E0")
    h0 = Symbol("h0")

    K0e=np.zeros((12,12))
    for i in range(12):
        for j in range(12):
            #PAMIĘTAJ, ŻE CALKA JEST A*B INTEGRATE _ WZORY DZIALAJA
            N = M.col(i).row(j)
            C = N.subs([(a, A), (b, B), (v, V), (E0, e0), (h0, H0)])
            CX = integrate(C, (x, 0, 1))
            CX1=CX*A
            CY = integrate(CX1, (y, 0, 1))
            CY1=CY*B

            if (CY1 != 0):
                K0e[i][j] = CY1[0,0]
            else:
                K0e[i][j] = 0

    return K0e

def MatWspAndMp(le,lw,n,m,a,b):
    mp = np.zeros((le, 4))
    for i in range(n):
        for j in range(m):
            k = (i) * m + j
            # U nas wzory inne niż u Kazia - bo i,j leci od 0 do x, A NIE OD 0
            mp[k, 0] = (i) * m + i + 1 + j
            mp[k, 1] = (i) * m + i + j + 2
            mp[k, 2] = (i + 1) * m + i + j + 3
            mp[k, 3] = (i + 1) * m + i + j + 2

    wsp = np.zeros((lw, 2))
    for i in range(n + 1):
        for j in range(m + 1):
            k = (i) * (m + 1) + j
            wsp[k, 0] = (j) * a
            wsp[k, 1] = (i) * b

    mp = np.int_(mp)



    return mp,wsp

def MatK0(lw,le,mp,K0e):

    K0 = np.zeros((lw * 3, lw * 3))
    for k in range(le):
        for i in range(4):
            for j in range(4):
                i1 = mp[k][i]
                j1 = mp[k][j]
                #print(i1, j1)
                for r in range(3):
                    for s in range(3):
                        K01 = (i1 - 1) * 3 + r
                        K02 = (j1 - 1) * 3 + s
                        K0e1 = (i + 1 - 1) * 3 + r
                        K0e2 = (j + 1 - 1) * 3 + s
                        #print(i1, j1, "X")

                        K0[K01][K02] = K0[K01][K02] + K0e[K0e1][K0e2]
                        # print(K01,K02,K0e1,K0e2)
                        # Pamiętaj- że normalnie leci od 1do 3 a u Ciebie od 0 do 2

    return K0


def MatK0UTW(K0,m,n):

    for i in range(m + 1):
        K0[3 * (i)][3 * (i)] = 10 ** 6
        K0[3 * (i) + 2][3 * (i) + 2] = 10 ** 6

    for i in range(n * m + n + 1, (n + 1) * m + n + 2):
        K0[3 * (i - 1)][3 * (i - 1)] = 10 ** 6
        K0[3 * (i - 1) + 2][3 * (i - 1) + 2] = 10 ** 6

    i = 0
    while (i < m * n + n + 1):
        K0[3 * (i - 1) + 3][3 * (i - 1) + 3] = 10 ** 6
        K0[3 * (i - 1) + 4][3 * (i - 1) + 4] = 10 ** 6
        i += m + 1

    i = m + 1
    while (i < (n + 1) * m + n + 2):
        K0[3 * (i - 1)][3 * (i - 1)] = 10 ** 6
        K0[3 * (i - 1) + 1][3 * (i - 1) + 1] = 10 ** 6
        i += m + 1

    return K0

def MatK0UTW2zaw(K0,m,n):

    for i in range(m + 1):
        K0[3 * (i)][3 * (i)] = 10 ** 6
        K0[3 * (i) + 2][3 * (i) + 2] = 10 ** 6

    for i in range(n * m + n + 1, (n + 1) * m + n + 2):
        K0[3 * (i - 1)][3 * (i - 1)] = 10 ** 6
        K0[3 * (i - 1) + 2][3 * (i - 1) + 2] = 10 ** 6


    return K0

def MatK0UTW3SW(K0,m,n):
   K0[0][0]=10**50
   x=(m+1)*(n+1)-1

   K0[m*3][m*3]=10**50
   K0[x*3][x*3]=10**50

   return K0


def Obciazenie(A,B,le,lw,mp):
    F0exy = ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]
    F0e = ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]
    a = Symbol("a")
    b = Symbol("b")
    x = Symbol("x")
    y = Symbol("y")
    q = Symbol("q")

    nxy = []
    for i in range(12):
        nxy.append("0")

    nxy[0] = (1 - x) * (1 - y) * (1 + x + y - 2 * x ** 2 - 2 * y ** 2)
    nxy[1] = b * (1 - x) * (1 - y) ** 2 * y
    nxy[2] = -a * (1 - y) * (1 - x) ** 2 * x
    nxy[3] = (1 - y) * (3 * x + y - 2 * x ** 2 - 2 * y ** 2) * x
    nxy[4] = b * (1 - y) ** 2 * x * y
    nxy[5] = a * (1 - x) * (1 - y) * x ** 2
    nxy[6] = (-1 + 3 * x + 3 * y - 2 * x ** 2 - 2 * y ** 2) * x * y
    nxy[7] = -b * (1 - y) * x * y ** 2
    nxy[8] = a * (1 - x) * x ** 2 * y
    nxy[9] = (1 - x) * (x + 3 * y - 2 * x ** 2 - 2 * y ** 2) * y
    nxy[10] = -b * (1 - x) * (1 - y) * y ** 2
    nxy[11] = -a * (1 - x) ** 2 * x * y

    for i in range(12):
        F0exy[i] = nxy[i] * q

    for i in range(12):
        C = nxy[i]
        CX = integrate(C, (x, 0, 1))
        CY = integrate(CX, (y, 0, 1))
        Temp1 = a * b * CY
        F0e[i] = Temp1.subs([(a, A), (b, B)])

    F0 = []
    for i in range(lw * 3):
        F0.append(0)
    i = 0

    for k in range(le):
        for i in range(4):
            i1 = mp[k][i]
            for r in range(3):
                F0[(i1 - 1) * 3 + r] = F0[(i1 - 1) * 3 + r] + F0e[(i) * 3 + r]

    return F0e,F0
def ObcPInit(lw):
    F0 = []
    for i in range(lw * 3):
        F0.append(0)
    i = 0
    return F0

def ObciazenieP(F0,P,nrwez):
    F0[nrwez*3]+=P
    return F0

def Q0(lw,K0,F0):

    Z = Matrix(F0)
    Y = np.matrix(K0)
    Y2=np.linalg.cholesky(Y)
    Y3=Matrix(Y2)
    X=Y3*Z
    return X


def Q0Chol(K0,F0):
    T0 = time.clock()
    L=np.linalg.cholesky(K0)
    T1 = time.clock()
    #print("Cholewski1 :",T1-T0)

    T0 = time.clock()
    y=np.linalg.solve(L,F0)
    T1 = time.clock()
    #print("Rozwiązanie Macierzy 1 :", T1 - T0)

    T0 = time.clock()
    LT=np.matrix.getT(L)
    T1 = time.clock()
    #print("Transponowanie Macierzy 1 :", T1 - T0)

    T0 = time.clock()
    X=np.linalg.solve(LT,y)
    T1 = time.clock()
    return X

def Ugiecie(Q0):
    U=[]
    I=len(Q0)/3
    for i in range(int(I)):
        U.append(Q0[3*i])
    return U

def TablicaUgiec(wsp,U):
    Uwsp=[]
    print(np.shape(wsp))
    i=0
    for i in wsp:

        d=wsp[i,0]

        print(type(U),np.shape(U),U[0],type(0),type(i))
        Uwsp.append([d[0],d[1],U[i]])

    return Uwsp

def UgieciaINT(U,m,n,pm,pn):
    x= m +1 + m*pm
    y= n+1 + n*pn

    Uk=np.zeros((y,x))
    d=0
    X,Y=0,0

    for i in range(n+1):
        for j in range(m+1):
            #print(x*i," ",y*j," ",i*j," ",U[i*j])
            Uk[(pn+1)*i][(pm+1)*j]=U[d]
            d=d+1
    return Uk

def DobudowanieUgiec(Uk,Q0,pm,pn,a0,b0,mp,nrel,m,n):
    T1=time.clock()
    nxy = []
    W1=mp[0]
    W2=mp[1]
    W3=mp[2]
    W4=mp[3]

    DeltaA=1/(pm+1)
    DeltaB=1/(pn+1)

    WymiarMacierzyX=(m*pm)+1
    WymiarMacierzyY = (n * pn) + 1


    w1, fix1, fiy1 = Q0[3 * (W1 - 1)], Q0[3 * (W1 - 1) + 1], Q0[3 * (W1 - 1) + 2]
    w2, fix2, fiy2 = Q0[3 * (W2 - 1)], Q0[3 * (W2 - 1) + 1], Q0[3 * (W2 - 1) + 2]
    w3, fix3, fiy3 = Q0[3 * (W3 - 1)], Q0[3 * (W3 - 1) + 1], Q0[3 * (W3 - 1) + 2]
    w4, fix4, fiy4 = Q0[3 * (W4 - 1)], Q0[3 * (W4 - 1) + 1], Q0[3 * (W4 - 1) + 2]
    T2=time.clock()
    #print(W1," ",W2," ",W3," ",W4)

    for i in range(12):
        nxy.append("0")


# MAMY TO - DOBRZE NAN WRZUCA _ TERAZ DORZUCIC TYLKO DYNAMICZNY ROZMIAR TABLICY I OGOLNIE UK
    #NIECH OBRACA A NIE TEMP


    Xprim = 0
    Yprim = 0

    DoXStep = (nrel % m) * (pm + 1)
    DoYStep = (np.int_(nrel / n)) * (pn + 1)
    DoX = DoXStep
    DoY = DoYStep

    for i in range(pn+2):
        for j in range(pm+2):
            nxy1 = AFK.n0(a0,b0,Xprim,Yprim)
            nxy2 =AFK.n1(a0,b0,Xprim,Yprim)
            nxy3 =AFK.n2(a0,b0,Xprim,Yprim)

            nxy4 =AFK.n3(a0,b0,Xprim,Yprim)
            nxy5 =AFK.n4(a0,b0,Xprim,Yprim)
            nxy6 =AFK.n5(a0,b0,Xprim,Yprim)

            nxy7 =AFK.n6(a0,b0,Xprim,Yprim)
            nxy8 =AFK.n7(a0,b0,Xprim,Yprim)
            nxy9 =AFK.n8(a0,b0,Xprim,Yprim)

            nxy10 =AFK.n9(a0,b0,Xprim,Yprim)
            nxy11 =AFK.n10(a0,b0,Xprim,Yprim)
            nxy12 =AFK.n11(a0,b0,Xprim,Yprim)

            Ug1 = w1 * nxy1 + fix1 * nxy2 + fiy1 * nxy3
            Ug2 = w2 * nxy4 + fix2 * nxy5 + fiy2 * nxy6
            Ug3 = w3 * nxy7 + fix3 * nxy8 + fiy3 * nxy9
            Ug4 = w4 * nxy10 + fix4 * nxy11 + fiy4 * nxy12
            Ugiecie=Ug1 + Ug2 + Ug3 + Ug4


            Uk[DoY][DoX]=Ugiecie
            DoX = DoX + 1
            Xprim = Xprim + DeltaA

        Yprim = Yprim + DeltaB
        Xprim=0
        DoY = DoY + 1
        DoX = DoXStep

    #print(TEMP)
    #print(w1," ",w2," ",w3," ",w4)
    #print(nrel,' ',T2-T1," ",T3-T2,' ', T4-T3)

    return Uk

def LiniaWpływu(Uk):
    I,J=Uk.shape
    UL=[]
    for i in range(I):
        for j in range(J):
            UL.append(Uk[i][j])
    return UL

def MacierzWIJ(lm,ln,E0,v,a,b):
    Wij = []
    M = lm
    N = ln
    I = 0
    J = 0
    MAXJ, MAXI = -N, -M
    Fki = CAL.TablicaFki(a,b,lm,ln)

    while (J > MAXJ):
        while (I > MAXI):
            # FORMULA
            ZEM2 = []
            i, j, m, n = I, J, M, N
            while (j < n):
                while (i < m):
                    ZEM2 = ZE.OblWij(Fki, ZEM2, abs(i), abs(j))
                    i = i + 1
                i = I
                j = j + 1
            Wij.append(ZEM2)
            # KONIEC
            I = I - 1
            M = M - 1
        I = 0
        M = -MAXI
        J = J - 1
        N = N - 1

    # ListaList[i]-lista osiadań wszystkich węzłów dla obciążonego i tego węzła
    wym=len(Wij)
    Wij2=np.zeros((wym,wym))
    Pivot=0
    x=0
    y=0
    for k in Wij:
        for l in k:
            Wij2[y][x]=l
            Wij2[y][x]=l*((1-v**2)/(3.141592653589793*E0))
            Pivot=Pivot+1
            if(x==wym-1):
                x=0
                y=y+1
            else:
                x=x+1


    return Wij2


def MacierzB(Yij,pm,pn,a,b):
    N,_=Yij.shape
    B=np.zeros((N+3,N+3))

    # Zapełnienie Macierzy B Yij
    for i in range(N):
        for j in range(N):
            B[i][j]=Yij[i][j]

    #Sily SigmaP
    for i in range(N):
        B[N][i]=1

    #Momenty MY
    X=a/2
    StepX=a
    MaxStepX=0
    for i in range(N):
        B[N+1][i]=X
        X=X+a
        MaxStepX=MaxStepX+1
        if(MaxStepX==pm):
            X=a/2
            MaxStepX=0
    #Koniec

    #Momenty MX
    Y = b/2
    StepY = b
    MaxStepY = 0
    for i in range(N):
        B[N+2][i]=Y
        MaxStepY=MaxStepY+1
        if(MaxStepY==pm):
            Y=Y+b
            MaxStepY=0
    #B[N + 1][i] = 2
    #B[N + 2][i] = 3
    return B

### Do przetestowania
def CzysteP(F0uzyt):
    P0uzyt=[]
    for i in range(len(F0uzyt)):
        if(i%3==0):
            P0uzyt.append(F0uzyt[i])
    return P0uzyt
### Do przetestowania
def P0INT(P0,m,n,pm,pn):
    x= m +1 + m*pm
    y= n+1 + n*pn

    PRO=np.zeros((y,x))
    d=0
    for i in range(n+1):
        for j in range(m+1):
            #print(x*i," ",y*j," ",i*j," ",U[i*j])
            l,p=(pn+1)*i,(pm+1)*j
            PRO[(pn+1)*i][(pm+1)*j]=P0[d]
            d=d+1

    Pk=[]

    for i in range(y):
        for j in range(x):
            Pk.append(PRO[i][j])

    Pk.append(0)
    Pk.append(0)
    Pk.append(0)

    # Ale czy PK ma odpowiedni kształt - czy to jest wektor ???  - jeśli nie - to transponuj
    return Pk

### Do przetestowania
def BPD(B,P,D):
    b=np.matrix(B)
    p=np.matrix(P)
    d=np.matrix(D)

    return (b*p-d)

def MacierzA(Wij,Yij,pm,pn,La,Lb):
    N,_=Yij.shape
    A=np.zeros((N+3,N+3))
    # Zapełnienie Macierzy B Yij
    for i in range(N):
        for j in range(N):
            A[i][j] = Yij[i][j]+Wij[i][j] ### JEDNOSTKI- JAKIE JEDNOSTKI??????
    # Sily SigmaP
    for i in range(N):
        A[N][i] = 1
    # Momenty MY
    X = 0
    StepX = 0.25
    MaxStepX = 0
    for i in range(N):
        A[N + 1][i] = X
        X = X + 0.25
        MaxStepX = MaxStepX + 1
        if (MaxStepX == pm):
            X = 0
            MaxStepX = 0
    # Koniec
    # Momenty MX
    Y = 0
    StepY = 0.25
    MaxStepY = 0
    for i in range(N):
        A[N + 2][i] = Y
        MaxStepY = MaxStepY + 1
        if (MaxStepY == pm):
            Y = Y + 0.25
            MaxStepY = 0

    #NO TO TERAZ ZAPELNIAMY PIOWNOWO

    X = 0
    StepX = 0.25
    MaxStepX = 0
    for i in range(N):
        A[i][N]= (X/La)
        A[i][N+1] = -X/La
        X = X + 0.25
        MaxStepX = MaxStepX + 1
        if (MaxStepX == pm):
            X = 0
            MaxStepX = 0
       ### PRZERWA - 3 kolumna
    Y = 0
    StepY = 0.25
    MaxStepY = 0
    for i in range(N):
        A[i][N]=A[i][N]+(Y/Lb)-1
        A[i][N+2] = -Y/Lb
        MaxStepY = MaxStepY + 1
        if (MaxStepY == pm):
            Y = Y + 0.25
            MaxStepY = 0


    return A


# Do Dopisania
def SzkodyGórnicze(R,kat,wezel):
    D=[]
    return D

def ZEM_SIN(A,B):


    T0 = time.clock()
    X = np.linalg.solve(A, B)
    T1 = time.clock()
    return X




