import numpy as np

def resize(sim, H_new, ctr):
    H_old_inv = np.array(sim.B);
    H_old = np.array(sim.H);
    F = np.transpose(H_old_inv @ H_new);
    sim.strain(F - np.identity(3));
    
    H_new_inv = np.array(sim.B);
    H_new = np.array(sim.H);
    F_inv = np.transpose(H_new_inv @ H_old);
    F = np.transpose(H_old_inv @ H_new);
    
    ctr_trnsf = np.dot(F, ctr)
    
    def retrieve(x):
        x[:]= np.dot(F_inv, x - ctr_trnsf) + ctr_trnsf;

    sim.do(retrieve);

def displace(sim,dx):
    def f(x):
        x+=dx;
    sim.do(f);

def normalize(x):
    return x/np.linalg.norm(x);
class Loop:
    def __init__(self,b,X_0,X_1,X_2):
        self.X=np.array([X_0,X_1,X_2]);
        self.t=np.array([normalize(X[2]-X[1]),\
        normalize(X[0]-X[2]),normalize(X[1]-X[0])]);
        self.bt=np.array([np.cross(b,t[0]),\
        np.cross(b,t[1]),np.cross(b,t[2])]);


def calc_nr(x,X,n,r):
    r[0]=np.linalg.norm(X[0]-x);
    n[0][:]=(X[0]-x)/r[0];
    r[1]=np.linalg.norm(X[1]-x);
    n[1][:]=(X[1]-x)/r[1];
    r[2]=np.linalg.norm(X[2]-x);
    n[2][:]=(X[2]-x)/r[2];
    

def calc_gs(n,b,gs):
    gs[:]=\
    np.dot(b,np.cross(n[1],n[2]))/(1.0+np.dot(n[1],n[2]))*0.5*(n[1]+n[2])+\
    np.dot(b,np.cross(n[2],n[0]))/(1.0+np.dot(n[2],n[0]))*0.5*(n[2]+n[0])+\
    np.dot(b,np.cross(n[0],n[1]))/(1.0+np.dot(n[0],n[1]))*0.5*(n[0]+n[1]);
def calc_fs(t,bt,n,r,fs):
    fs[:]=\
    bt[0]*np.log((r[2]*(1.0+np.dot(n[2],t[0])))/(r[1]*(1.0+np.dot(n[1],t[0]))))+\
    bt[1]*np.log((r[0]*(1.0+np.dot(n[0],t[1])))/(r[2]*(1.0+np.dot(n[2],t[1]))))+\
    bt[2]*np.log((r[1]*(1.0+np.dot(n[1],t[2])))/(r[0]*(1.0+np.dot(n[0],t[2]))));


def solid_angle(N,n):
    a=np.array([np.arccos(np.dot(n[1],n[2])),np.arccos(np.dot(n[2],n[0])),np.arccos(np.dot(n[0],n[1]))]);
    s=0.5*(a[0]+a[1]+a[2]);
    return -(4.0*np.arctan(np.sqrt(np.tan(0.5*s)*np.tan(0.5*(s-a[0]))*np.tan(0.5*(s-a[1]))*np.tan(0.5*(s-a[2])))))*np.sign(np.dot(n[0],N));

def get_u(nu,b,X,t,bt,x):
    n=np.full((3,3),0.0);
    r=np.full((3),0.0);
    fs=np.full((3),0.0);
    gs=np.full((3),0.0);
    calc_nr(x,X,n,r);
    calc_gs(n,b,gs);
    calc_fs(t,bt,n,r,fs);
    return(-solid_angle(N,n)*b/4.0*np.pi)\
    (gs-(1.0-2.0*nu)*fs)/(8.0*np.pi*(1.0-nu));

def res_shear0(b,s):
    tau=np.linalg.norm(b)*np.linalg.norm(s);
    Q=np.array([normalize(b),normalize(s),normalize(np.cross(b,s))]);
    QT=np.transpose(Q);
    F=np.identity(3);
    F[0][1]=tau;
    E=0.5*(np.dot(np.transpose(F),F)-np.identity(3))
    E=np.full((3,3),0.0);
    E[0][1]=E[1][0]=0.5*tau;
    #print(F)
    #print(E)
    return np.dot(np.dot(QT,E),Q)

def res_shear(b,s):
    return np.tensordot(s,b,0)+np.tensordot(b,s,0)




#__p=np.dot(o,s)
#__t=__p-np.floor(__p)
#o_lo=o-__t*s/np.dot(s,s)
#o_up=o+(1.0-__t)*s/np.dot(s,s)


def matrix_2_tensor(C_m):
    v = np.array([[0,0], [1,1], [2,2], [1,2], [2,0], [0,1]]);
    C = np.full((3,3,3,3),0.0);
    for i in range(0,3):
        C[i, i, i, i] = C_m[i, i];
        C[v[i+3, 0], v[i+3, 1], v[i+3, 0], v[i+3, 1]] = \
        C[v[i+3, 1], v[i+3, 0], v[i+3, 0], v[i+3, 1]] = \
        C[v[i+3, 0], v[i+3, 1], v[i+3, 1], v[i+3, 0]] = \
        C[v[i+3, 1], v[i+3, 0], v[i+3, 1], v[i+3, 0]] = \
        C_m[3+i, 3+i];
    for i in range(0,3):
        for j in range(i+1,3):
            C[i, i, j, j] = C[j, j, i, i] = C_m[i, j];
            C[v[i+3, 0], v[i+3, 1], v[j+3, 0], v[j+3, 1]] = \
            C[v[i+3, 1], v[i+3, 0], v[j+3, 0], v[j+3, 1]] = \
            C[v[i+3, 0], v[i+3, 1], v[j+3, 1], v[j+3, 0]] = \
            C[v[i+3, 1], v[i+3, 0], v[j+3, 1], v[j+3, 0]] = \
            C[v[j+3, 0], v[j+3, 1], v[i+3, 0], v[i+3, 1]] = \
            C[v[j+3, 0], v[j+3, 1], v[i+3, 1], v[i+3, 0]] = \
            C[v[j+3, 1], v[j+3, 0], v[i+3, 0], v[i+3, 1]] = \
            C[v[j+3, 1], v[j+3, 0], v[i+3, 1], v[i+3, 0]] = \
            C_m[3+i, 3+j];
    for i in range(0,3):
        for j in range(0,3):
            C[i, i, v[j+3, 0], v[j+3, 1]] = \
            C[i, i, v[j+3, 1], v[j+3, 0]] = \
            C[v[j+3, 0], v[j+3, 1], i, i] = \
            C[v[j+3, 1], v[j+3, 0], i, i] = \
            C_m[i, 3+j];

    return C;

def tensor_2_matrix(C):
    v = [[0,0], [1,1], [2,2], [1,2], [2,0], [0,1]];
    C_m = np.full((6,6), 0.0);
    for i in range(0,3):
        C_m[i, i] = C[i, i, i, i];
        C_m[3+i, 3+i] = C[v[i+3, 0], v[i+3, 1], v[i+3, 0], v[i+3, 1]];
    for i in range(0,3):
        for j in range(i+1,3):
            C_m[i, j] = C_m[j, i]=C[i, i, j, j];
            C_m[3+i, 3+j] = C_m[3+j, 3+i] = C[v[i+3, 0], v[i+3, 1], v[j+3, 0], v[j+3, 1]];
    for i in range(0,3):
        for j in range(0,3):
            C_m[i, 3+j] = C_m[3+j, i] = C[i, i, v[j+3, 0], v[j+3, 1]];
    return C_m;

def rot(T, Q, depth = 1):
    if len(T.shape) > depth:
        return np.tensordot(rot(T, Q, depth+1), Q, axes=([0],[1]))
    else: 
        return np.tensordot(T, Q, axes=([0],[1]))

def get_str(C,sh):
    sts=np.full((6),0.0);
    sts[5]=sh;
    
    str=np.linalg.solve(C,sts)
    FT_F=[
          [str[0],str[5]*0.5,str[4]*0.5],
          [str[5]*0.5,str[1],str[3]*0.5],
          [str[4]*0.5,str[3]*0.5,str[2]]
          ]+np.identity(3);



    E=np.full((3,3),0.0);
    E[2][2]=np.sqrt(FT_F[2][2])
    E[2][1]=FT_F[2][1]/E[2][2];
    E[2][0]=FT_F[2][0]/E[2][2];
    E[1][1]=np.sqrt(FT_F[1][1]-E[2][1]*E[2][1])
    E[1][0]=(FT_F[1][0]-E[2][0]*E[2][1])/E[1][1];
    E[0][0]=np.sqrt(FT_F[0][0]-E[1][0]*E[1][0]-E[2][0]*E[2][0]);
    E-=np.identity(3);
    return E;


def poly_modulus(C):
    S=np.linalg.inv(C);
    G_0=15.0/(\
              4.0*(S[0][0]+S[1][1]+S[2][2])-\
              4.0*(S[1][2]+S[2][0]+S[0][1])+\
              3.0*(S[3][3]+S[4][4]+S[5][5])\
              );
    B_0=1.0/(
             S[0][0]+S[1][1]+S[2][2]+\
             2.0*(S[1][2]+S[2][0]+S[0][1])\
             );
    G_1=(
         C[0][0]+C[1][1]+C[2][2]-\
         (C[1][2]+C[2][0]+C[0][1])+\
         3.0*(C[3][3]+C[4][4]+C[5][5])
         )/15.0;
    B_1=(C[0][0]+C[1][1]+C[2][2]+2.0*(C[1][2]+C[2][0]+C[0][1]))/9.0;
    G=0.5*(G_0+G_1);
    B=0.5*(B_0+B_1);
    Y=9.0*B*G/(3.0*B+G);
    nu=(1.5*B-G)/(3.0*B+G);
    return [nu,Y,B,G];


def cubic(c11,c12,c44):
    C = np.full((6,6),0.0);
    C[0, 0] = C[1, 1] = C[2][2] = c11
    C[1, 2] = C[2, 0] = C[0][1] = \
    C[2, 1] = C[0, 2] = C[1][0] = c12
    C[3, 3] = C[4, 4] = C[5, 5] = c44;
    return matrix_2_tensor(C);


def sth(C, b, threshold=1.0e-9):
    B = np.stack([C[:, 0, :, 0], C[:, 0, :, 1] + C[:, 1, :, 0], C[:, 1, :, 1]])
    def _(i, j, k):
        return np.linalg.det(np.stack([B[i][0], B[j][1], B[k][2]]))
    
    coefs = np.array([
        _(2, 2, 2),
        
        _(1, 2, 2) + _(2, 1, 2) + _(2, 2, 1),
        
        _(0, 2, 2) + _(2, 0, 2) + _(2, 2, 0) + \
        _(2, 1, 1) + _(1, 2, 1) + _(1, 1, 2),
        
        _(0, 1, 2) + _(2, 0, 1) + _(1, 2, 0) + \
        _(0, 2, 1) + _(1, 0, 2) + _(2, 1, 0) + \
        _(1, 1, 1), 
        
        _(0, 1, 1) + _(1, 0, 1) + _(1, 1, 0) + \
        _(2, 0, 0) + _(0, 2, 0) + _(0, 0, 2),
        
        _(1, 0, 0) + _(0, 1, 0) + _(0, 0, 1), 
        
        _(0, 0, 0)
    ])
    
    
    ps=np.roots(coefs)
    
    for i in range(0,6,2):
        for j in range(5,i,-1):
            if np.linalg.norm(ps[i]-np.conj(ps[j]))<1.0e-9:
                t=ps[i+1];
                ps[i+1]=ps[j];
                ps[j]=t;
        if np.imag(ps[i])<0.0:
            t=ps[i];
            ps[i]=ps[i+1];
            ps[i+1]=t;

    ps=np.array([ps[0],ps[2],ps[4]]);
    _ = np.stack([np.full((3), 1.0),ps, ps*ps])
    a = np.tensordot(_, B, axes=([0],[0])) 
    
    A= np.full((3,3), 0.0, dtype=np.complex128);
    for i in range(0,3):
        A[0, i]= a[i, 0, 1] * a[i, 1, 2] - a[i, 1, 1] * a[i, 0, 2];
        A[1, i]= a[i, 0, 2] * a[i, 1, 0] - a[i, 1, 2] * a[i, 0, 0];
        A[2, i]= a[i, 0, 0] * a[i, 1, 1] - a[i, 1, 0] * a[i, 0, 1];
        for j in range(0,3):
            if np.linalg.norm(A[j, i])<threshold:
                A[j, i]=0.0;
        if A[2, i]!=0.0:
            A[1, i]/=A[2, i];
            A[0, i]/=A[2, i];
            A[2, i]=1.0;
        elif A[1, i]!=0.0:
            A[0, i]/=A[1, i];
            A[1, i]=1.0;
        elif A[0, i]!=0.0:
            A[0, i]=1.0;
        else:
            A[0, i]=1.0;
            A[1, i]=-a[i, 0, 0]/a[i, 0, 1];

    B=np.full((3,3),0.0,dtype=np.complex128);
    for i in range(0,3):
        for n in range(0,3):
            B[i, n]=\
            (C[i, 1, 0, 0]+ps[n]*C[i, 1, 0, 1])*A[0, n]+\
            (C[i, 1, 1, 0]+ps[n]*C[i, 1, 1, 1])*A[1, n]+\
            (C[i, 1, 2, 0]+ps[n]*C[i, 1, 2, 1])*A[2, n];

    Krn = np.full((6, 6), 0.0);
    bb = np.full((6),0.0);
    for i in range(0,3):
        bb[i] = b[i] / np.linalg.norm(b);
        for j in range(0,3):
            Krn[i, j] = np.real(A[i, j]);
            Krn[i, j+3] = -np.imag(A[i, j])
            Krn[i+3, j] = np.real(B[i, j]);
            Krn[i+3, j+3] = -np.imag(B[i, j])

    DD = np.linalg.solve(Krn, bb);
    T=np.full((3,3), 0.0, dtype=np.complex128);
    for i in range(0,3):
        for j in range(0,3):
            T[i, j]=A[i, j] * (DD[j] + DD[j+3] * 1j);            
    return ps,T    

class HirthEdge:
    def __init__(self, C, b):
        self.C = C;
        self.b = b;
        self.ps, self.T = sth(C, b)
    
    def ave_disp_helper(self, z):
        return np.pi*np.imag(z)-np.log(2.0)-np.pi*np.sign(np.imag(z))*np.real(z)*1j;

    def disp_helper(self, z):
        n=1.0+np.floor(np.real(z)-0.5)
        return np.log(np.sin((z-n)*np.pi))-np.sign(np.imag(z))*n*np.pi*1j;
    
    def disp(self, x):
        return np.real(np.dot(self.T,self.disp_helper(x[0]*np.full((3),1.0)+x[1]*self.ps))/(-2.0*np.pi*1j))
    def ave_disp(self, x):
        return np.array([
            (0.5*x[0]-0.25)*np.sign(x[1]),
            np.real(np.dot(self.T,self.ave_disp_helper(np.full((3),-0.5)+np.abs(x[1])*self.ps))/(-2.0*np.pi*1j))[1],
            0.0])
        
class HirthScrew:
    def __init__(self, C, b):
        self.C = C;
        self.b = b;
        self.ps, self.T = sth(C, b)
    
    def ave_disp_helper(self, z):
        return np.pi*np.imag(z)-np.log(2.0)-np.pi*np.sign(np.imag(z))*np.real(z)*1j;

    def disp_helper(self, z):
        n=1.0+np.floor(np.real(z)-0.5)
        return np.log(np.sin((z-n)*np.pi))-np.sign(np.imag(z))*n*np.pi*1j;
    
    def disp(self, x):
        return np.real(np.dot(self.T,self.disp_helper(x[0]*np.full((3),1.0)+x[1]*self.ps))/(-2.0*np.pi*1j)) \
            - np.array([0.0,0.0,0.5*x[0]])
    def ave_disp(self, x):
        return np.array([0.0, 
        self.disp(np.array([0.25,x[1],0.0]))[1],
        -0.25 if x[1] > 0.0 else 0.25 - x[0]])
        
    
def crack(C):
    B = np.linalg.inv(
        np.array([
        [C[0, 0, 0, 0], C[0, 0, 1, 1], C[0, 0, 0, 1]],
        [C[0, 0, 1, 1], C[1, 1, 1, 1], C[1, 1, 0, 1]],
        [C[0, 0, 0, 1], C[1, 1, 0, 1], C[0, 1, 0, 1]]
        ]
    ))

    _ = np.roots([B[0, 0], -2.0*B[0, 2],2.0*B[0, 1]+B[2, 2], -2.0*B[1, 2], B[1, 1]])

    mu = np.array([_[0],0.0]);

    if np.absolute(np.conjugate(mu[0]) - _[1]) > 1.0e-12:
        mu[1] = _[1];
    else:
        mu[1] = _[2]

    alpha = np.real(mu);
    beta = np.imag(mu);

    p = B[0,0] * mu**2 - B[0,2] * mu + B[0, 1]
    q = B[0,1] * mu - B[0, 2] + B[1, 1]/ mu

    _ = np.stack([p, q]) * np.array(mu[1], mu[0]) /(mu[1] - mu[0])
    K_r = np.real(_)
    K_i = np.imag(_)
 
    Tr = np.stack([
        np.array(np.array([[1.0, alpha[0]], [0.0, beta[0]]])), 
        np.array([[1.0, alpha[1]], [0.0, beta[1]]])
    ], axis=1)

    def u_f0(x): return np.sqrt(np.sqrt(x[0] * x[0] + x[1] * x[1]) + x[0])
    def u_f1(x): return np.sqrt(np.sqrt(x[0] * x[0] + x[1] * x[1]) - x[0]) * np.sign(x[1]) 

    def disp(x): 
        _ = Tr @ x[:2]
        return K_r @ u_f0(_) + K_i @ u_f1(_)
    
    return disp