import numpy as np
from sympy import *
import Kernel as KE
import Zemoczkin as ZE
import time
import Kernel2 as KE2


def Oblicz(*args):
    a0, b0, h0, v, E0, E0G, vG, m, n, TabelaF, txt, TabelaR = args
    a0 = float(a0)
    b0 = float(b0)
    h0 = float(h0)
    v = float(v)
    E0 = float(E0)
    E0G = float(E0G)
    vG = float(vG)
    m = int(m)
    n = int(n)

    a = a0 / m
    b = b0 / n

    T0 = time.clock()
    K1 = KE.fK0()
    K0e = KE.SymIntegrate(K1, a, b, v, E0, h0)

    le = m * n
    lw = (m + 1) * (n + 1)

    mp, wsp = KE.MatWspAndMp(le, lw, n, m, a, b)

    K0 = KE.MatK0(lw, le, mp, K0e)
    K0 = KE.MatK0UTW3SW(K0, m, n)
    T1 = time.clock()

    print("Koniec wstępnych czynności - teraz są pętle", T1 - T0)

    wym = le
    Yij = np.zeros((wym, wym))
    for i in range(le):
        w1, w2, w3, w4 = mp[i]


        T0 = time.clock()
        F0 = KE.ObcPInit(lw)
        F0 = KE.ObciazenieP(F0, 0.25, w1 - 1)
        F0 = KE.ObciazenieP(F0, 0.25, w2 - 1)
        F0 = KE.ObciazenieP(F0, 0.25, w3 - 1)
        F0 = KE.ObciazenieP(F0, 0.25, w4 - 1)
        Q0 = KE.Q0Chol(K0, F0)

        X = 0
        Y = 0
        Uk = np.zeros((n, m))
        for j in range(le):
            Uk = KE2.PrzebudowanieUgięc(Uk, Q0, a, b, mp[j][:], j, m, n, X, Y)
            X = X + 1
            if (X == m):
                X = 0
                Y = Y + 1


        UL = KE.LiniaWpływu(Uk)

        Yij[i][:] = UL
        T1 = time.clock()
        print("Pętla: ", i + 1, " trwała :", T1 - T0)


    ####DOBRA TERAZ LECIMY Z JAKIMŚ PRZYKLADOWYM OBCIAZENIEM
    F0 = KE.ObcPInit(le)


    def Przylozdo1(X, Y, mp, wsp, Ele):
        w1, w2, w3, w4 = mp
        W1X, W1Y = wsp[w1 - 1]
        W3X, W3Y = wsp[w3 - 1]
        if (X >= W1X and Y >= W1Y and X < W3X and Y < W3Y):
            return True
        return false

    def TablicaF(X, Y, mp, wsp, F0, WarP):
        for Ele in range(le):
            Out = Przylozdo1(X, Y, mp[Ele], wsp, Ele)
            if (Out == true):
                print("Przyłożono do Ele nr: ", Ele, " o War", WarP)
                F0 = KE.ObciazenieP(F0, WarP, Ele)
                return F0


    F = TabelaF

    max, _ = np.shape(F)
    for i in range(max):
        TablicaF(F[i][0], F[i][1], mp, wsp, F0, F[i][2])

    P0uzyt = KE.CzysteP(F0)
    # print(len(P0uzyt), P0uzyt)
    P0uzyt = KE2.P0INT(P0uzyt, m, n, 0, 0)
    print("P0", len(P0uzyt), P0uzyt)

    # PODZIAL WEZLOW po ile 25 cm bloczków jest
    pm = m
    pn = n

    Binit = KE.MacierzB(Yij, pm, pn, a, b)

    B = np.matrix(Binit)

    P = np.matrix(P0uzyt)
    P = np.transpose(P)

    X2 = B * P
    B = np.zeros((P.shape))
    X, l = X2.shape

    for i in range(X):
        B[i][0] = X2[i][0]

    for i in range(len(TabelaR)):
        XXX = float(TabelaR[i][1])
        YYY = float(TabelaR[i][2])
        alpha = float(TabelaR[i][3])
        R = float(TabelaR[i][0])
        Gornictwo = KE2.mapa(wsp, mp, R, alpha, XXX, YYY)

        for j in range(len(Gornictwo)):
            B[j][0] -= Gornictwo[j]



    Wij = KE.MacierzWIJ(m, n, E0G, vG, a, b)

    A = KE2.MacierzA(Wij, Yij, m, n, a, b, a0, b0)

    BRZA_MROZ = np.linalg.lstsq(A, B,0)

    WYNIK = BRZA_MROZ[0]

    TEM = 0
    for i in range(n + 1):
        print("%7.4f " % (i - 1), " ", end="")
    print('yhm')

    Poz = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (WYNIK[Poz]), " ", end="")
            TEM += WYNIK[Poz]
            Poz = Poz + 1
        print()
    Tem = 0
    Wyn = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            Wyn[i][j] = WYNIK[Tem]
            Tem += 1

    print("U1,U2,U3")
    print(WYNIK[-3:])
    print(TEM)

    #### Do ZWROTU ####
    Ugiecia = WYNIK[-3:]

    ###############################################################################################

    Wosiadania = KE2.MacierzOsiadań(Wyn, m, n, E0G, vG, a, b)

    min = 100000000000
    max = -10000000000
    for i in range(n):
        for j in range(m):
            if (Wosiadania[i][j] <= min):
                min = Wosiadania[i][j]
            if (Wosiadania[i][j] >= max):
                max = Wosiadania[i][j]

    print("WOSIADANIA")
    for i in range(n):
        print()
        print(" ", end="")

        for j in range(m):
            print("%7.4f " % (Wosiadania[i][j]), " ", end="")

    ##################################################################################################

    # Policzenie płyty MESEM, ale z siłami zadanymi

    F0 = KE.ObcPInit(lw)

    #### Siły zadane ####
    def Przylozdo1(X, Y, mp, wsp, Ele):
        w1, w2, w3, w4 = mp
        W1X, W1Y = wsp[w1 - 1]
        W3X, W3Y = wsp[w3 - 1]
        if (X >= W1X and Y >= W1Y and X < W3X and Y < W3Y):
            return True
        return false

    def TablicaF2(X, Y, mp, wsp, F0, WarP):
        for Ele in range(le):
            Out = Przylozdo1(X, Y, mp[Ele], wsp, Ele)
            if (Out == true):
                # print("Przyłożono do Ele nr: ", Ele)
                w1, w2, w3, w4 = mp[Ele]
                F0 = KE.ObciazenieP(F0, WarP / 4, w1 - 1)
                F0 = KE.ObciazenieP(F0, WarP / 4, w2 - 1)
                F0 = KE.ObciazenieP(F0, WarP / 4, w3 - 1)
                F0 = KE.ObciazenieP(F0, WarP / 4, w4 - 1)
                return F0

    print("No Dobra22222")
    F = TabelaF
    max, _ = np.shape(F)
    for i in range(max):
        TablicaF2(F[i][0], F[i][1], mp, wsp, F0, F[i][2])

    print("TAK WYGLADA WEKTOR PO PRZYLOZENIU F ZEW")
    print(F0)

    ##### Siły odporu #####

    for i in range(m * n):
        w1, w2, w3, w4 = mp[i]
        ODPOR = WYNIK[i][0] / 4
        F0 = KE.ObciazenieP(F0, -ODPOR, w1 - 1)
        F0 = KE.ObciazenieP(F0, -ODPOR, w2 - 1)
        F0 = KE.ObciazenieP(F0, -ODPOR, w3 - 1)
        F0 = KE.ObciazenieP(F0, -ODPOR, w4 - 1)
        # print("Przyłożono odpror w ele ",i, " o sile", ODPOR *4 )

    print("TAK WYGLADA WEKTOR PO PRZYLOZENIU Odporow")
    print(F0)
    Q0 = KE.Q0Chol(K0, F0)
    UFIN = np.zeros((n + 1, m + 1))
    Mark = 0
    for i in range(n + 1):
        for j in range(m + 1):
            UFIN[i][j] = Q0[Mark]
            Mark += 3

    print(len(Q0))

    ######################################################################################
    # Osiadania

    print("Deformacja")
    OSI = np.zeros((n, m))
    X = 0
    Y = 0

    for j in range(le):
        OSI = KE2.PrzebudowanieUgięc(OSI, Q0, a, b, mp[j][:], j, m, n, X, Y)
        X = X + 1
        if (X == m):
            X = 0
            Y = Y + 1

    for i in range(n + 1):
        print("%7.4f " % (i - 1), " ", end="")
    print()

    Tem = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (OSI[i][j]), " ", end="")
            Tem += 1
        print()

    print("Różnica")
    Tem = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (OSI[i][j] - Wosiadania[i][j]), " ", end="")
            Tem += 1
        print()

    print("MXX")
    MXX = np.zeros((n, m))
    for i in range(le):
        MXX = KE2.WyznaczenieMx(MXX, Q0, 0, 0, a, b, mp, i, m, n, v)

    for i in range(n):
        for j in range(m):
            MXX[i][j] = -MXX[i][j] * E0 * h0 ** 3 / (12 * (1 - v ** 2))

    TEM = 0
    for i in range(n + 1):
        print("%7.4f " % (i - 1), " ", end="")
    print()

    Tem = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (MXX[i][j]), " ", end="")
            Tem += 1
        print()

    #############################################################################3
    print("MYY")
    MYY = np.zeros((n, m))
    for i in range(le):
        MYY = KE2.WyznaczenieMy(MYY, Q0, 0, 0, a, b, mp, i, m, n, v)

    for i in range(n):
        for j in range(m):
            MYY[i][j] = -MYY[i][j] * E0 * h0 ** 3 / (12 * (1 - v ** 2))

    TEM = 0
    for i in range(n + 1):
        print("%7.4f " % (i - 1), " ", end="")
    print()

    Tem = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (MYY[i][j]), " ", end="")
            Tem += 1
        print()

        #############################################################################3

    MXY = np.zeros((n, m))
    for i in range(le):
        MXY = KE2.WyznaczenieMxy(MXY, Q0, 0, 0, a, b, mp, i, m, n, v)

    for i in range(n):
        for j in range(m):
            MXY[i][j] = -MXY[i][j] * E0 * (1 - v) * h0 ** 3 / (12 * (1 - v ** 2))

    TEM = 0
    for i in range(n + 1):
        print("%7.4f " % (i - 1), " ", end="")
    print()

    Tem = 0
    for i in range(n):
        print("%7.4f " % (i), " ", end="")
        for j in range(m):
            print("%7.4f " % (MXY[i][j]), " ", end="")
            Tem += 1
        print()

    KE2.drukuj(Wyn, Wosiadania, MXX, MYY, MXY, a0, b0, h0, m, n, E0, v, E0G, vG, txt)

    return Wyn, Wosiadania, MXX, MXY, MYY
