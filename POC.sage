def Gen_UOV(v,o,N):#-------------Generate a UOV layer 
    vp=v
    op=o
    n=v+o
    K=GF(q)
    P=PolynomialRing(K,'x',n)
    x=P.gens()
    x_vec=vector(x)
    VV=[0 for i in range(op)]
    VO=[0 for i in range(op)]
    Bb=[0 for i in range(op)]
    B=[[] for i in range(op)]
    A=[0 for i in range(op)]
    D=[0 for i in range(op)]
    for i in range(op):
        VO[i]=random_matrix(K,vp,op)
    for i in range(op):
        VV[i]=random_matrix(K,vp,vp)
    for i in range(op):
        A[i]=block_matrix([[VV[i],VO[i]],[matrix(K,op,vp),matrix(K,op,op)]])
    F=[i for i in range(o)]
    for i in range(o):
        F[i]=matrix(K,N,N)
        F[i].set_block(0,0,A[i])
    return F

def Gen_MB(v,d,op,N):#-------------------------Generate a MB UOV layer.
    o=d*op
    n=v+o
    K.<a>=GF(q)
    P=PolynomialRing(K,'x',n)
    x=P.gens()
    x_vec=vector(x)
    A00=[0 for i in range(o)]
    A01=[0 for i in range(o)]
    B0=[0 for i in range(o)]
    B=[[] for i in range(o)]
    A=[0 for i in range(o)]
    D=[0 for i in range(o)]
    D[0]=random_matrix(K,v,v)
    Ar=[0 for i in range(op)]
    for r in range(op):
        A01[r]=[]
        for i in range(v*op):
            A01[r].append(K.random_element())
        Ar[r]=matrix(K,op,v,A01[r]).transpose()
        for i in range(v*(o-op)):
            A01[r].append(K(0))
        A01[r]=matrix(K,o,v,A01[r]).transpose()

    rotate_matrix=[]
    for i in range(op):
        for j in range(o):
            if(j==(o-op+i)):
                rotate_matrix.append(K(1))
            else:
                rotate_matrix.append(K(0))
    for i in range(o-op):
        for j in range(o):
            if(j==i):
                rotate_matrix.append(K(1))
            else:
                rotate_matrix.append(K(0))
    rotate_matrix=matrix(K,o,o,rotate_matrix).transpose()
    for i in range(o):
        A01[i]=A01[i%op]*rotate_matrix^(i//op)

    for i in range(o):
        A00[i]=D[i]
    for i in range(o):
        A[i]=matrix.block([[A00[i],A01[i]],[matrix(K,o,v),matrix(K,o,o)]])
    F=[0 for i in range(o)]
    for i in range(o):
        F[i]=A[i]#-------------------central quadratic matrix
    return F
v1=3
o1=4
d=2
op=4
o2=op*d
N=v1+o1+o2
m=o1+o2
v2=v1+o1
q=31
K=GF(q)
F1=Gen_UOV(v1,o1,N)
F2=Gen_MB(v1+o1,d,op,N)

S=identity_matrix(K,N)
B1=random_matrix(K,v1,o1)
B2=random_matrix(K,v1,o2)
B3=random_matrix(K,o1,o2)
S.set_block(0,v1,B1)
S.set_block(0,v2,B2)
S.set_block(o1,v2,B3)
Z=matrix(K,o1+o2,1)
Sn=S
Sn.set_block(0,N-1,Z)
T=identity_matrix(K,o1+o2)
T1=random_matrix(K,o1,o2)
T.set_block(0,o1,T1)
for i in range(o1,o1+o2,1):
    T[o1-1,i]=0

F=[]
for i in range(o1):
    F.append(S.transpose()*F1[i]*S)
for i in range(o2):
    F.append(S.transpose()*F2[i]*S)

P=[]
for i in range(o1+o2):
    tmp=0
    for j in range(o1+o2):
        tmp=tmp+T[i,j]*F[j]
    P.append(tmp+tmp.transpose())
for i in range(o1+o2):
    print i
    print P[i]
    print ""
