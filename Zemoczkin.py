import numpy as np
from sympy import *
import time as time
import CalkowanieWIj as CAL

# Siatka Zwraca tablicę MxN kwadracików o wymiarach 0,25x0,25 w formie listy
def TabWsp():

    Fki=np.zeros((8,8))
    Fki[0][0] = 2.974
    Fki[0][1] = 1.111
    Fki[0][2] = 0.511
    Fki[0][3] = 0.336
    Fki[0][4] = 0.251
    Fki[0][5] = 0.201
    Fki[0][6] = 0.167
    Fki[0][7] = 0.143

    Fki[1][0] = 1.111
    Fki[1][1] = 0.749
    Fki[1][2] = 0.455
    Fki[1][3] = 0.319
    Fki[1][4] = 0.244
    Fki[1][5] = 0.197
    Fki[1][6] = 0.165

    Fki[2][0] = 0.511
    Fki[2][1] = 0.455
    Fki[2][2] = 0.357
    Fki[2][3] = 0.279
    Fki[2][4] = 0.225
    Fki[2][5] = 0.186

    Fki[3][0] = 0.336
    Fki[3][1] = 0.319
    Fki[3][2] = 0.279
    Fki[3][3] = 0.237
    Fki[3][4] = 0.201

    Fki[4][0] = 0.251
    Fki[4][1] = 0.244
    Fki[4][2] = 0.225
    Fki[4][3] = 0.201


    Fki[5][0] = 0.201
    Fki[5][1] = 0.197
    Fki[5][2] = 0.186

    Fki[6][0] = 0.167
    Fki[6][1] = 0.165

    Fki[7][0] = 0.143

    return Fki

def TabWsp2(m,n,a,b):

    Fki=np.zeros((m,n))
    Fki[0][0]=CAL.SamNaSiebie(a,b)
    return Fki

#Funkcja przyjmuje tablcię Fki oraz ile jednostek w prawo oraz w góre od punktu przyłożenia obciążenia jest dany sektor
#Rozumiem, że Tam, gdzie obciążenie w metodzie wpływu, tam jest osiadanie od samego siebie 2,974
def OblWij(Fki,ZEM,wprawo,wgore):
   ZEM.append(Fki[wgore][wprawo])
   return ZEM


def UgiecieEl(Fki,R,wprawo,wgóre,lm,ln):
    U=0
    X=-wprawo
    Y=-wgóre

    for i in range(ln):
        for j in range(lm):
            U+=R[i][j]*Fki[abs(Y)][abs(X)]
            X+=1
        Y+=1
        X=-wprawo
    return U




