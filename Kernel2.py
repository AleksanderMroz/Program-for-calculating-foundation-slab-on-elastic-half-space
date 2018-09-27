import time as time

import numpy as np
from sympy import *

init_printing(use_unicode=True)
import Zemoczkin as ZE
import AproxFunkKsztal as AF
import Dx2 as AFK
import Dy2 as AFK2
import Dxy as AFK3
import CalkowanieWIj as CAL


def WyznaczenieMy(Uk, Q0, pm, pn, a0, b0, mp, nrel, m, n, v):
    nxy = []
    for i in range(12):
        nxy.append("0")
    W1, W2, W3, W4 = mp[nrel]
    w1, fix1, fiy1 = Q0[3 * (W1 - 1)], Q0[3 * (W1 - 1) + 1], Q0[3 * (W1 - 1) + 2]
    w2, fix2, fiy2 = Q0[3 * (W2 - 1)], Q0[3 * (W2 - 1) + 1], Q0[3 * (W2 - 1) + 2]
    w3, fix3, fiy3 = Q0[3 * (W3 - 1)], Q0[3 * (W3 - 1) + 1], Q0[3 * (W3 - 1) + 2]
    w4, fix4, fiy4 = Q0[3 * (W4 - 1)], Q0[3 * (W4 - 1) + 1], Q0[3 * (W4 - 1) + 2]

    Xprim = 0.5
    Yprim = 0.5

    DoXStep = (nrel % m) * 1
    DoYStep = (np.int_(nrel / m)) * 1
    DoX = DoXStep
    DoY = DoYStep

    nxy1 = AFK.n0(a0, b0, Xprim, Yprim) + v * AFK2.n0(a0, b0, Xprim, Yprim)
    nxy2 = AFK.n1(a0, b0, Xprim, Yprim) + v * AFK2.n1(a0, b0, Xprim, Yprim)
    nxy3 = AFK.n2(a0, b0, Xprim, Yprim) + v * AFK2.n2(a0, b0, Xprim, Yprim)

    nxy4 = AFK.n3(a0, b0, Xprim, Yprim) + v * AFK2.n3(a0, b0, Xprim, Yprim)
    nxy5 = AFK.n4(a0, b0, Xprim, Yprim) + v * AFK2.n4(a0, b0, Xprim, Yprim)
    nxy6 = AFK.n5(a0, b0, Xprim, Yprim) + v * AFK2.n5(a0, b0, Xprim, Yprim)

    nxy7 = AFK.n6(a0, b0, Xprim, Yprim) + v * AFK2.n6(a0, b0, Xprim, Yprim)
    nxy8 = AFK.n7(a0, b0, Xprim, Yprim) + v * AFK2.n7(a0, b0, Xprim, Yprim)
    nxy9 = AFK.n8(a0, b0, Xprim, Yprim) + v * AFK2.n8(a0, b0, Xprim, Yprim)

    nxy10 = AFK.n9(a0, b0, Xprim, Yprim) + v * AFK2.n9(a0, b0, Xprim, Yprim)
    nxy11 = AFK.n10(a0, b0, Xprim, Yprim) + v * AFK2.n10(a0, b0, Xprim, Yprim)
    nxy12 = AFK.n11(a0, b0, Xprim, Yprim) + v * AFK2.n11(a0, b0, Xprim, Yprim)

    Ug1 = w1 * nxy1 + fix1 * nxy2 + fiy1 * nxy3
    Ug2 = w2 * nxy4 + fix2 * nxy5 + fiy2 * nxy6
    Ug3 = w3 * nxy7 + fix3 * nxy8 + fiy3 * nxy9
    Ug4 = w4 * nxy10 + fix4 * nxy11 + fiy4 * nxy12
    Ugiecie = Ug1 + Ug2 + Ug3 + Ug4

    Uk[DoY][DoX] = Ugiecie
    return Uk


def WyznaczenieMx(Uk, Q0, pm, pn, a0, b0, mp, nrel, m, n, v):
    nxy = []
    for i in range(12):
        nxy.append("0")
    W1, W2, W3, W4 = mp[nrel]
    w1, fix1, fiy1 = Q0[3 * (W1 - 1)], Q0[3 * (W1 - 1) + 1], Q0[3 * (W1 - 1) + 2]
    w2, fix2, fiy2 = Q0[3 * (W2 - 1)], Q0[3 * (W2 - 1) + 1], Q0[3 * (W2 - 1) + 2]
    w3, fix3, fiy3 = Q0[3 * (W3 - 1)], Q0[3 * (W3 - 1) + 1], Q0[3 * (W3 - 1) + 2]
    w4, fix4, fiy4 = Q0[3 * (W4 - 1)], Q0[3 * (W4 - 1) + 1], Q0[3 * (W4 - 1) + 2]

    Xprim = 0.5
    Yprim = 0.5

    DoXStep = (nrel % m) * 1
    DoYStep = (np.int_(nrel / m)) * 1
    DoX = DoXStep
    DoY = DoYStep

    nxy1 = AFK2.n0(a0, b0, Xprim, Yprim) + v * AFK.n0(a0, b0, Xprim, Yprim)
    nxy2 = AFK2.n1(a0, b0, Xprim, Yprim) + v * AFK.n1(a0, b0, Xprim, Yprim)
    nxy3 = AFK2.n2(a0, b0, Xprim, Yprim) + v * AFK2.n2(a0, b0, Xprim, Yprim)

    nxy4 = AFK2.n3(a0, b0, Xprim, Yprim) + v * AFK.n3(a0, b0, Xprim, Yprim)
    nxy5 = AFK2.n4(a0, b0, Xprim, Yprim) + v * AFK.n4(a0, b0, Xprim, Yprim)
    nxy6 = AFK2.n5(a0, b0, Xprim, Yprim) + v * AFK.n5(a0, b0, Xprim, Yprim)

    nxy7 = AFK2.n6(a0, b0, Xprim, Yprim) + v * AFK.n6(a0, b0, Xprim, Yprim)
    nxy8 = AFK2.n7(a0, b0, Xprim, Yprim) + v * AFK.n7(a0, b0, Xprim, Yprim)
    nxy9 = AFK2.n8(a0, b0, Xprim, Yprim) + v * AFK.n8(a0, b0, Xprim, Yprim)

    nxy10 = AFK2.n9(a0, b0, Xprim, Yprim) + v * AFK.n9(a0, b0, Xprim, Yprim)
    nxy11 = AFK2.n10(a0, b0, Xprim, Yprim) + v * AFK.n10(a0, b0, Xprim, Yprim)
    nxy12 = AFK2.n11(a0, b0, Xprim, Yprim) + v * AFK.n11(a0, b0, Xprim, Yprim)

    Ug1 = w1 * nxy1 + fix1 * nxy2 + fiy1 * nxy3
    Ug2 = w2 * nxy4 + fix2 * nxy5 + fiy2 * nxy6
    Ug3 = w3 * nxy7 + fix3 * nxy8 + fiy3 * nxy9
    Ug4 = w4 * nxy10 + fix4 * nxy11 + fiy4 * nxy12
    Ugiecie = Ug1 + Ug2 + Ug3 + Ug4

    Uk[DoY][DoX] = Ugiecie
    return Uk


def WyznaczenieMxy(Uk, Q0, pm, pn, a0, b0, mp, nrel, m, n, v):
    nxy = []
    for i in range(12):
        nxy.append("0")
    W1, W2, W3, W4 = mp[nrel]

    w1, fix1, fiy1 = Q0[3 * (W1 - 1)], Q0[3 * (W1 - 1) + 1], Q0[3 * (W1 - 1) + 2]
    w2, fix2, fiy2 = Q0[3 * (W2 - 1)], Q0[3 * (W2 - 1) + 1], Q0[3 * (W2 - 1) + 2]
    w3, fix3, fiy3 = Q0[3 * (W3 - 1)], Q0[3 * (W3 - 1) + 1], Q0[3 * (W3 - 1) + 2]
    w4, fix4, fiy4 = Q0[3 * (W4 - 1)], Q0[3 * (W4 - 1) + 1], Q0[3 * (W4 - 1) + 2]

    Xprim = 0.5
    Yprim = 0.5

    DoXStep = (nrel % m) * 1
    DoYStep = (np.int_(nrel / m)) * 1
    DoX = DoXStep
    DoY = DoYStep

    nxy1 = AFK3.n0(a0, b0, Xprim, Yprim)
    nxy2 = AFK3.n1(a0, b0, Xprim, Yprim)
    nxy3 = AFK3.n2(a0, b0, Xprim, Yprim)

    nxy4 = AFK3.n3(a0, b0, Xprim, Yprim)
    nxy5 = AFK3.n4(a0, b0, Xprim, Yprim)
    nxy6 = AFK3.n5(a0, b0, Xprim, Yprim)

    nxy7 = AFK3.n6(a0, b0, Xprim, Yprim)
    nxy8 = AFK3.n7(a0, b0, Xprim, Yprim)
    nxy9 = AFK3.n8(a0, b0, Xprim, Yprim)

    nxy10 = AFK3.n9(a0, b0, Xprim, Yprim)
    nxy11 = AFK3.n10(a0, b0, Xprim, Yprim)
    nxy12 = AFK3.n11(a0, b0, Xprim, Yprim)

    Ug1 = w1 * nxy1 + fix1 * nxy2 + fiy1 * nxy3
    Ug2 = w2 * nxy4 + fix2 * nxy5 + fiy2 * nxy6
    Ug3 = w3 * nxy7 + fix3 * nxy8 + fiy3 * nxy9
    Ug4 = w4 * nxy10 + fix4 * nxy11 + fiy4 * nxy12
    Ugiecie = Ug1 + Ug2 + Ug3 + Ug4

    Uk[DoY][DoX] = Ugiecie
    return Uk


def PrzebudowanieUgięc(Uk, Q0, a0, b0, mp, nrel, m, n, X, Y):
    T1 = time.clock()
    nxy = []
    # print(nrel,mp[nrel])
    # print(type(mp[0]))
    W1 = mp[0]
    W2 = mp[1]
    W3 = mp[2]
    W4 = mp[3]

    w1, fix1, fiy1 = Q0[3 * (W1 - 1)], Q0[3 * (W1 - 1) + 1], Q0[3 * (W1 - 1) + 2]
    w2, fix2, fiy2 = Q0[3 * (W2 - 1)], Q0[3 * (W2 - 1) + 1], Q0[3 * (W2 - 1) + 2]
    w3, fix3, fiy3 = Q0[3 * (W3 - 1)], Q0[3 * (W3 - 1) + 1], Q0[3 * (W3 - 1) + 2]
    w4, fix4, fiy4 = Q0[3 * (W4 - 1)], Q0[3 * (W4 - 1) + 1], Q0[3 * (W4 - 1) + 2]

    nxy1 = AF.n0(a0, b0, 0.5, 0.5)
    nxy2 = AF.n1(a0, b0, 0.5, 0.5)
    nxy3 = AF.n2(a0, b0, 0.5, 0.5)

    nxy4 = AF.n3(a0, b0, 0.5, 0.5)
    nxy5 = AF.n4(a0, b0, 0.5, 0.5)
    nxy6 = AF.n5(a0, b0, 0.5, 0.5)

    nxy7 = AF.n6(a0, b0, 0.5, 0.5)
    nxy8 = AF.n7(a0, b0, 0.5, 0.5)
    nxy9 = AF.n8(a0, b0, 0.5, 0.5)

    nxy10 = AF.n9(a0, b0, 0.5, 0.5)
    nxy11 = AF.n10(a0, b0, 0.5, 0.5)
    nxy12 = AF.n11(a0, b0, 0.5, 0.5)

    Ug1 = w1 * nxy1 + fix1 * nxy2 + fiy1 * nxy3
    Ug2 = w2 * nxy4 + fix2 * nxy5 + fiy2 * nxy6
    Ug3 = w3 * nxy7 + fix3 * nxy8 + fiy3 * nxy9
    Ug4 = w4 * nxy10 + fix4 * nxy11 + fiy4 * nxy12
    Ugiecie = Ug1 + Ug2 + Ug3 + Ug4
    # print(Ugiecie, nxy1,nxy2,nxy3,nxy4,nxy5,nxy6,nxy7,nxy8,nxy9,nxy10,nxy11,nxy12)

    # print(Y,X)
    # print(nrel, Ugiecie)
    Uk[Y][X] = Ugiecie

    return Uk


def MacierzA(Wij, Yij, pm, pn, La, Lb, a, b):
    N, _ = Yij.shape
    A = np.zeros((N + 3, N + 3))
    # Zapełnienie Macierzy B Yij
    for i in range(N):
        for j in range(N):
            A[i][j] = Yij[i][j] + Wij[i][j]  ### JEDNOSTKI- JAKIE JEDNOSTKI??????
    # Sily SigmaP
    for i in range(N):
        A[N][i] = 1
    # Momenty MY
    X = La / 2
    StepX = La
    MaxStepX = 0
    for i in range(N):
        A[N + 1][i] = X
        X = X + La
        MaxStepX = MaxStepX + 1
        if (MaxStepX == pm):
            X = La / 2
            MaxStepX = 0
    # Koniec
    # Momenty MX
    Y = Lb / 2
    StepY = Lb
    MaxStepY = 0
    for i in range(N):
        A[N + 2][i] = Y
        MaxStepY = MaxStepY + 1
        if (MaxStepY == pm):
            Y = Y + Lb
            MaxStepY = 0

    # NO TO TERAZ ZAPELNIAMY PIOWNOWO

    X = La / 2
    StepX = 0.25
    MaxStepX = 0
    for i in range(N):
        A[i][N] = (X / a)
        A[i][N + 1] = -X / a
        X = X + La
        MaxStepX = MaxStepX + 1
        if (MaxStepX == pm):
            X = La / 2
            MaxStepX = 0
    ### PRZERWA - 3 kolumna
    Y = Lb / 2
    StepY = Lb
    MaxStepY = 0
    for i in range(N):
        A[i][N] = A[i][N] + (Y / b) - 1
        A[i][N + 2] = -Y / b
        MaxStepY = MaxStepY + 1
        if (MaxStepY == pm):
            Y = Y + Lb
            MaxStepY = 0

    return A


def ObcPInit(le):
    F0 = []
    for i in range(le * 3):
        F0.append(0)
    i = 0
    return F0


def ObciazenieP(F0, P, nrwez):
    F0[nrwez * 3] = P
    return F0


def CzysteP(F0uzyt):
    P0uzyt = []
    for i in range(len(F0uzyt)):
        if (i % 3 == 0):
            P0uzyt.append(F0uzyt[i])
    return P0uzyt


### Do przetestowania


def P0INT(P0, m, n, pm, pn):
    x = m
    y = n

    PRO = np.zeros((y, x))
    d = 0
    for i in range(n):
        for j in range(m):
            # print(x*i," ",y*j," ",i*j," ",U[i*j])

            PRO[(1) * i][(1) * j] = P0[d]
            d = d + 1

    Pk = []

    for i in range(y):
        for j in range(x):
            Pk.append(PRO[i][j])

    Pk.append(0)
    Pk.append(0)
    Pk.append(0)

    # Ale czy PK ma odpowiedni kształt - czy to jest wektor ???  - jeśli nie - to transponuj
    return Pk


def MacierzOsiadań(R, lm, ln, E0, v, a, b):
    Fki = CAL.TablicaFki(a, b, lm, ln)

    U = np.zeros((ln, lm))

    wprawo = 0
    wgore = 0
    for i in range(ln):
        for j in range(lm):
            U[i][j] = ((1 - v ** 2) / (3.141592653589793 * E0)) * ZE.UgiecieEl(Fki, R, wprawo, wgore, lm, ln)
            wprawo += 1
        wprawo = 0
        wgore += 1

    return U


def MacierzOsiadańMOD(Rbasic, lm, ln, E0, v, a, b):
    Fki = CAL.TablicaFki(a, b, lm + 2, ln + 2)
    R = np.zeros((ln + 2, lm + 2))
    for i in range(ln):
        for j in range(lm):
            R[i + 1][j + 1] = Rbasic[i][j]

    U = np.zeros((ln + 2, lm + 2))

    wprawo = 0
    wgore = 0
    for i in range(ln + 2):
        for j in range(lm + 2):
            U[i][j] = ((1 - v ** 2) / (3.141592653589793 * E0)) * ZE.UgiecieEl(Fki, R, wprawo, wgore, lm + 2, ln + 2)
            wprawo += 1
        wprawo = 0
        wgore += 1

    return U


def Dx(W, a):
    n, m = np.shape(W)

    Out = np.zeros((n, m - 1))

    for i in range(n):
        for j in range(m - 1):
            Out[i][j] = (W[i][j + 1] - W[i][j]) / a
    return Out


def Dxx(W, a):
    Out = Dx(W, a)
    Out2 = Dx(Out, a)
    n, m = np.shape(Out2)
    Out3 = np.zeros((n - 2, m))
    for i in range(n - 2):
        for j in range(m):
            Out3[i][j] = Out2[i + 1][j]
    return Out3


def Dxy(W, a, b):
    Out = Dx(W, a)
    Out2 = Dy(Out, b)
    n, m = np.shape(Out2)
    Out3 = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            Out3[i][j] = Out2[i][j]
    return Out3


def Dy(W, b):
    n, m = np.shape(W)

    Out = np.zeros((n - 1, m))

    for i in range(n - 1):
        for j in range(m):
            Out[i][j] = (W[i + 1][j] - W[i][j]) / b
    return Out


def Dyy(W, b):
    Out = Dy(W, b)
    Out2 = Dy(Out, b)

    n, m = np.shape(Out2)
    Out3 = np.zeros((n, m - 2))
    for i in range(n):
        for j in range(m - 2):
            Out3[i][j] = Out2[i][j + 1]
    return Out3


def MYY(dxx, dyy, E, h, v):
    n, m = np.shape(dxx)
    MXX = np.zeros((n, m))
    D = E * h ** 3 / (12 * (1 - v ** 2))
    for i in range(n):
        for j in range(m):
            MXX[i][j] = -D * (dxx[i][j] + v * dyy[i][j])
    return MXX


def MXX(dxx, dyy, E, h, v):
    n, m = np.shape(dxx)
    MYY = np.zeros((n, m))
    D = E * h ** 3 / (12 * (1 - v ** 2))
    for i in range(n):
        for j in range(m):
            MYY[i][j] = -D * (v * dxx[i][j] + dyy[i][j])
    return MYY


def MXY(dxy, E, h, v):
    n, m = np.shape(dxy)
    MXY = np.zeros((n, m))
    D = E * h ** 3 / (12 * (1 - v ** 2))
    for i in range(n):
        for j in range(m):
            MXY[i][j] = -(1 - v) * (D) * dxy[i][j]
    return MXY


def drukuj(R, w, MXX, MYY, MXY, a0, b0, h0, m, n, E0, v, E0G, vG, txt):
    Raport = "****************************** Plyta Fundamentowa ******************************\n"
    Raport += "Polprzestrzen sprezysta: ni: "
    Raport += str(vG)
    Raport += "    E0: "
    Raport += str(E0G)
    Raport += " MPa\n"
    Raport += "********************************************************************************\n"
    Raport += "Dane:\n"
    Raport += "Plyta:\n"
    Raport += " a = " + str(a0) + " m \n"
    Raport += " b = " + str(b0) + " m \n"
    Raport += " h = " + str(h0) + " m \n"
    Raport += " Podzial boku a na m = " + str(m) + " elementow \n"
    Raport += " Podzial boku b na n = " + str(n) + " elementow \n"
    Raport += "********************************************************************************\n\n\n"
    Raport += "*************************************Reakcje************************************\n"
    Raport += "********************(wielkosci sa wyrazone w [m] oraz [MN])*********************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:
            TEMP = ""
            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                Raport += "%7.4f" % R[i][j - 1] + "  "
        Raport += "\n"
    Raport += "********************************************************************************\n\n\n"

    Raport += "*************************************Naprezenia************************************\n"
    Raport += "********************(wielkosci sa wyrazone w [m] oraz [MPa])*********************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:
            TEMP = ""
            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                X = R[i][j - 1] / (SX * SY)
                Raport += "%7.4f" % X + "  "
        Raport += "\n"
    Raport += "********************************************************************************\n\n\n"

    Raport += "************************************Osiadania***********************************\n"
    Raport += "***********************(wielkosci sa wyrazone w [m] )***************************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:

            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                Raport += "%7.4f" % w[i][j - 1] + "  "
        Raport += "\n"
    Raport += "********************************************************************************\n\n\n"

    Raport += "**********************************Momenty MXX***********************************\n"
    Raport += "***********************(wielkosci sa wyrazone w [MNm] )***************************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:

            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                Raport += "%7.4f" % MXX[i][j - 1] + "  "
        Raport += "\n"
    Raport += "********************************************************************************\n\n\n"

    Raport += "**********************************Momenty MYY***********************************\n"
    Raport += "***********************(wielkosci sa wyrazone w [MNm] )***************************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:

            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                Raport += "%7.4f" % MYY[i][j - 1] + "  "
        Raport += "\n"
    Raport += "********************************************************************************\n\n\n"

    Raport += "**********************************Momenty MXY***********************************\n"
    Raport += "***********************(wielkosci sa wyrazone w [MNm] )***************************\n"
    SX = a0 / m
    SY = b0 / n
    SXX = a0 / (2 * m)
    SYY = b0 / (2 * n)
    for i in range(m + 1):
        if (i == 0):
            print(" b \ a ", " ", end="")
            TEMP = " b \ a " + " "
            Raport += TEMP
        else:

            TEMP = "%7.4f " % SXX + "|"
            SXX += SX
            print(TEMP, " ", end="")
            Raport += TEMP
    Raport += "\n"
    for i in range(n):
        for j in range(m + 1):
            if (j == 0):
                TEMP = "%7.4f" % SYY + "|"
                SYY += SY
                Raport += TEMP

            else:
                Raport += "%7.4f" % MXY[i][j - 1] + "  "
        Raport += "\n"

    Raport += "********************************************************************************\n\n\n\n"
    Raport += "Autor pracy magisterskiej: inz. Aleksander Mroz\n"
    Raport += "Promotor: dr hab. inz. Wlodzimierz Brzakala \n"
    Raport += "Semestr letni \n"
    Raport += "Rok akademicki: 2017/2018 \n"

    print(Raport)

    #plik = open(txt, "w")
    #plik.write(Raport)
   # plik.close()


def mapa(wsp, mp, R, alpha, X, Y):
    # print(np.sin(alpha))

    Ry = np.zeros((3, 3))
    Ry[0][0] = np.cos(alpha)
    Ry[0][2] = np.sin(alpha)
    Ry[1][1] = 1
    Ry[2][0] = -np.sin(alpha)
    Ry[2][2] = np.cos(alpha)

    wsps = np.zeros((len(mp), 2))
    for i in range(len(mp)):
        w1, w2, w3, w4 = mp[i] - 1
        x1, y1 = wsp[w1]
        x2, y2 = wsp[w2]
        x3, y3 = wsp[w3]
        x4, y4 = wsp[w4]
        SX = X + (x1 + x2 + x3 + x4) / 4
        SY = Y + (y1 + y2 + y3 + y4) / 4
        wsps[i][0] = SX
        wsps[i][1] = SY

    wsps2 = np.zeros((len(mp), 2))
    for i in range(len(mp)):
        Wektor = np.zeros((3, 1))

        Wektor[0][0] = wsps[i][0]
        Wektor[1][0] = 0
        Wektor[2][0] = wsps[i][1]
        # print(Wektor)
        Wynik = np.dot(Ry, Wektor)
        wsps2[i][0] = Wynik[0][0]
        wsps2[i][1] = Wynik[2][0]
        # print(Wynik,"\n\n\n\n")

    Max = -999999999999999999999

    PR = []
    for i in range(len(wsps2)):
        x1, z1 = wsps2[i]
        Tem = np.arccos(x1 / R)
        y1 = R * np.sin(Tem)
        if (y1 > Max):
            Max = y1
        PR.append(y1)

    for i in range(len(PR)):
        PR[i] = Max - PR[i]

    return PR


def Walec(R, x1, y1, x2, y2):
    print("NOPE")
    pass
