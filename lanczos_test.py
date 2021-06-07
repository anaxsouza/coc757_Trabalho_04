import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

n = 100
#rng = np.random.default_rng(12345)
# print(rng)
B = np.random.default_rng().uniform(0, 1, (n, n))
Q, R = np.linalg.qr(B)
D = np.diag(np.arange(1, n+1, 1))
QT = np.transpose(Q)
A = Q @ D @ QT

'''
print('Q')
print(Q)
print('R')
print(R)

print('D')
print(D)

print('QT')
print(QT)

print('A')
print(A)
'''

Alfa = np.zeros(n)
Beta = np.zeros(n-1)

q0 = 0  # initialize
b0 = 0  # initialize
# arbitrary nonzero starting vector
x0 = np.random.default_rng().uniform(0, 1, n)
qk = x0 / np.linalg.norm(x0)  # normalize
ritz = np.zeros((n,n))

for k in range(1, n+1):
    uk = A @ qk  # generate next vector
    ak = np.conj(qk) @ uk
    uk = uk - b0*q0 - ak*qk
    bk = np.linalg.norm(uk)
    b0 = bk
    Alfa[(k-1)] = ak
    if k != n:
        Beta[(k-1)] = bk
    q0 = qk
    qk = uk/bk

    T = np.zeros((k, k))
    cont1 = 0
    cont2 = 0
    cont3 = 0

    for i in range(k):
        for j in range(k):
            if i == j:
                T[i][j] = Alfa[cont1]
                cont1 = cont1 + 1
            if (i-j) == -1:
                T[i][j] = Beta[cont2]
                cont2 = cont2 + 1
            if (i-j) == 1:
                T[i][j] = Beta[cont3]
                cont3 = cont3 + 1


    w, _ = np.linalg.eig(T)
    #print(w)

    for m in range(w.shape[0]):
        ritz[(k-1)][m] = w[m]
    #ritz[(k-1)] = w

cont1 = 0
#print(ritz)
val = int((n*(1+n))/2)
x_ = np.zeros(val)
y_ = np.zeros(val)
for i in range(n):
    for j in range(n):
        if i >= j:
            x_[cont1] = ritz[i][j]
            y_[cont1] = i+1
            cont1 +=1
        
print(x_)
print(y_)

fig_1 = plt.figure(figsize=(16, 16))
ax1 = fig_1.subplots(1, 1)
#ax1.plot(ritz, 'k.')
ax1.plot(x_,y_,'k.')
ax1.set_xlabel('Ritz values', fontsize=12)
ax1.set_ylabel('iteration', fontsize=12)

#plt.savefig('figure_1.png')
plt.show()

'''
print('\n\nAlfa')
print(Alfa)

print('\n\nBeta')
print(Beta)

T = np.zeros((n, n))
cont1 = 0
cont2 = 0
cont3 = 0

for i in range(n):
    for j in range(n):
        if i == j:
            T[i][j] = Alfa[cont1]
            cont1 = cont1 + 1
        if (i-j) == -1:
            T[i][j] = Beta[cont2]
            cont2 = cont2 + 1
        if (i-j) == 1:
            T[i][j] = Beta[cont3]
            cont3 = cont3 + 1

print('\nT')
print(T)

w, v = np.linalg.eig(T)
print('E-value:', w)
print('E-vector', v)
'''