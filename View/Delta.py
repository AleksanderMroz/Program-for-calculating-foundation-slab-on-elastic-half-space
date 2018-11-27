
import numpy as np
import sympy as sp



D = np.zeros((6,6))


MP= np.zeros((6,5))
MK= np.zeros((6,5))
MS= np.zeros((6,5))

LoverEI = np.zeros((1,5))
LoverEI[0][0]=5/(6*2)
LoverEI[0][1]=2/(6*2)
LoverEI[0][2]=3/(6*2)
LoverEI[0][3]=4/(6*1.5)
LoverEI[0][4]=4/(6*1.5)

## q1
MP[0][0]=-7
MK[0][0]=-2

MP[0][1]=-2
MK[0][1]=-0

## q2
MP[1][0]=-10
MK[1][0]=-5

MP[1][1]=-5
MK[1][1]=-3

MP[1][2]=-3
MK[1][2]=0


## q3
MP[2][0]=-14
MK[2][0]=-9

MP[2][1]=-9
MK[2][1]=-7

MP[2][2]=-7
MK[2][2]=-4

MP[2][3]=-4
MK[2][3]=0



## q4
MP[3][0]=-1
MK[3][0]=-1

MP[3][1]=-1
MK[3][1]=-1

MP[3][2]=-1
MK[3][2]=-1


##x5
MP[4][0]=5
MK[4][0]=0

##x6
MP[5][0]=4
MK[5][0]=4

MP[5][1]=4
MK[5][1]=4

MP[5][2]=4
MK[5][2]=4

MP[5][4]=4
MK[5][4]=0

for i in range(6):
    for j in range(5):
        MS[i][j]=MP[i][j]*0.5 + MK[i][j]*0.5


for i in range(6):
    for j in range(6):
        D[i][j]= sum(LoverEI[0]*MP[i]*MP[j]+4*LoverEI[0]*MS[i]*MS[j]+LoverEI[0]*MK[i]*MK[j])

print(MP)
print(MS)
print(MK)
print(D)

DQQ= np.zeros((4,4))
DXX=np.zeros((2,2))
DXQ=np.zeros((2,4))

for i in range(4):
    for j in range(4):
        DQQ[i][j]=D[i][j]

for i in range(2):
    for j in range(2):
        DXX[i][j]=D[i+4][j+4]


for i in range(2):
    for j in range(4):
        DXQ[i][j]=D[i+4][j]


DQX=np.transpose(DXQ)

print(DQQ)
print(DXX)
print(DXQ)
print(DQX)
print()
print()
DXXi=np.linalg.inv(DXX)


D1=DQQ-DQX.dot(DXXi.dot(DXQ))



K=np.zeros((4,4))
K[0][0]=0.648
K[0][1]=-0.506
K[0][2]=0
K[0][3]=0.952


K[1][1]=0.498
K[1][2]=-0.07
K[1][3]=-0.753

K[2][2]=0.07
K[2][3]=-0.281

K[3][3]=3.325

for i in range(4):
    for j in range(4):
        if(K[i][j]==0):
            K[i][j]=K[j][i]

print("Macierz D1 red")
print(D1)


print("DIJ")
for i in range (6):
    for j in range (6):
        print(i+1,j+1,D[i][j])
    print()