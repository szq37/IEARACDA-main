import numpy as np
import copy
def DC(D,mu,T0):
    U,S,V = np.linalg.svd(D)
    T1 = np.zeros(np.size(T0))
    for i in range(1,100):
        T1 = DCInner(S,mu,T0)
        err = np.sum(np.square(T1-T0))
        if err < 1e-6:
            break
        T0 = T1
    J_1 = np.dot(U, np.diag(T1))
    J = np.dot(J_1, V)
    return J,T1
def DCInner(S,mu,T_k):
    lamb = 1/mu
    grad = 1/(1+np.square(T_k))
    T_k1 = S-lamb*grad
    T_k1[T_k1<0]=0
    return T_k1
def errorsol(Y_1,A_1,X_1,lamb,mu,type):
    if type == 1:
        D=A_1-A_1@X_1+(Y_1 / mu)
        E=np.zeros(np.shape(D))
        eps = lamb/mu
        DD = np.abs(D)-eps
        DD2= DD*np.sign(D)
        ID = np.abs(D)>eps
        E[ID]=DD2[ID]
    else:
        alpha = lamb /mu
        G = A_1-A_1@X_1+(Y_1 / mu)
        G1=np.sqrt(np.sum(np.square(G),1))
        G1[G1 == 0] = alpha
        G2 = (G1 - alpha)/ G1
        E = np.dot(G ,np.diag((G1 > alpha)* G2))
    return E
#主函数
def ARA(A,mu):
    lamb = 1e-3
    rho = 1.1
    tol = 1e-6
    m, n = np.shape(A)  
    X = np.zeros((n, n))
    E = np.zeros((m, n))  
    Y1 = np.zeros((m, n))  
    Y2 = np.zeros((n, n))  
#开始迭代更新
    for i in range(0, 1000):
        D = X+(Y2 / mu)
        sig = np.zeros(min(m, n))  
        J, sig = DC(copy.deepcopy(D), mu, copy.deepcopy(sig))
        X=update_X(I,A,E,J,Y1,Y2,mu)
        E= errorsol(copy.deepcopy(Y1),copy.deepcopy(A),copy.deepcopy(X),lamb,mu,type)
        Y1= Y1+mu*(A-A@X-E)
        Y2=Y2+mu*(X-J)
        mu = mu*rho
        sigma = np.linalg.norm(A-A@X-E,'fro')
        RRE = sigma/np.linalg.norm(X-J,'fro')
        if RRE < tol:
            break