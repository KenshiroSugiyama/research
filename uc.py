# coding: UTF-8
import math
import numpy as np
import warnings
warnings.simplefilter('ignore', category=RuntimeWarning)  # RuntimeWarningを無視扱いに設定
import re
import csv
import time
import matplotlib.pyplot as plt

def TDMAsolver(a, b, c, d):
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def main():
    
    n=1000

    #配列作成
    z = np.empty(n+1) 
    zm = np.empty(n+1)
    nu_t = np.empty(n+1) 
    nu_tm = np.empty(n+1) 
    #---------------------
    nu_tcm = np.empty(n+1) 
    cof = np.empty(n+1)
    #---------------------
    u = np.empty(n+1) 
    uo = np.empty(n+1) 
    um = np.empty(n+1) 
    D_u = np.empty(n+1) 
    #---------------------    
    dudzm = np.empty(n+1)
    dudz = np.empty(n+1)
    #---------------------   
    a_u = np.empty(n+1)
    b_u = np.empty(n+1)   
    c_u = np.empty(n+1)
    r_u = np.empty(n+1)
    #---------------------
    c = np.empty(n+1)
    co = np.empty(n+1)
    cm = np.empty(n+1)
    zeta = np.empty(n+1)
    Xm = np.empty(n+1)
    fum = np.empty(n+1)
    #---------------------
    dcdzm = np.empty(n+1)
    dcdz = np.empty(n+1)
    err_u = np.empty(n+1)
    err_c = np.empty(n+1)
    #---------------------
    Ri = np.empty(n+1)
    #---------------------
    kappa = 0.4
    sgm_c = 1.3
    
    z = np.empty(n+1)  
    u = np.empty(n+1)  

    # print(z)
    vs=1.00

    S=0.1
    ds=1.0/(1.0*10**2)
    ks=2.5*ds
    zb=1.0*ds/12
    ub=1.0/kappa*np.log(30.0*zb/ks)
    Rit=1.0/S

    dt=1.0*10**(-11)
    dz=(1.0-zb)/n
    for i in range(0,n+1):
        zm[i]=zb+dz*i
    z[1:n+1]=zm[1:n+1]-0.5*dz

    #-------- Initial velocity	-----------------------
    
    u = 1.0/kappa*np.log(30.0*z/ks)
    
    c=1.0-z-zb
    
    maxerr=1.0
    cb = 0
    start_time = time.perf_counter()
    while maxerr>1*10**(-6):
        #------------  gradient	-------------

        uo[1:n+1]=u[1:n+1]
        um[0]=ub
        um[1:n]=(uo[1:n]+uo[2:n+1])/2
        um[n]=uo[n]

        co[1:n+1]=c[1:n+1]
        cm[0]=cb
        cm[1:n]=(co[1:n]+co[2:n+1])/2
        cm[n]=0

        dudzm[0]=2*(uo[1]-ub)/dz
        dudzm[1:n]=(uo[2:n+1]-uo[1:n])/dz
        dudzm[n]=2*(um[n]-uo[n])/dz
        dudz[1:n+1]=(um[1:n+1]-um[0:n])/dz

        dcdzm[0]=2*(co[1]-cm[0])/dz
        dcdzm[1:n]=(co[2:n+1]-co[1:n])/dz
        dcdzm[n]=2*(cm[n]-co[n])/dz
        dcdz[1:n+1]=(cm[1:n+1]-cm[0:n])/dz
    
        #------eddy viscosity ----------------
        nu_tm[0:n+1]=kappa**2*zm[0:n+1]**2*(1-zm[0:n+1])*abs(dudzm[0:n+1])
        
        nu_t[1:n+1]=(nu_tm[0:n]+nu_tm[1:n+1])/2

        #------- 	velocity	----------------------

        D_u[0:n+1]=(nu_tm[0:n+1])/dz

        a_u[1]=0
        a_u[2:n+1]=-D_u[1:n]
        b_u[1]=dz/dt+2*D_u[0]+D_u[1]
        b_u[2:n]=dz/dt+D_u[1:n-1]+D_u[2:n]
        b_u[n]=dz/dt+D_u[n-1]
        c_u[1:n]=-D_u[1:n]
        c_u[n]=0
        r_u[1]=(uo[1]/dt+c[1])*dz+2*D_u[0]*ub
        r_u[2:n+1]=(uo[2:n+1]/dt+c[2:n+1])*dz

        # print(u)
        # u = tridag(a_u,b_u,c_u,r_u,u,n) #?サブルーチン呼び出し
        u[1:n+1] = TDMAsolver(a_u[2:n+1],b_u[1:n+1],c_u[1:n],r_u[1:n+1])

        #------------- concentration --------------------------

        zeta[0]=0
        for i in range(1,n+1):
        	zeta[i]=zeta[i-1]+vs/nu_t[i]*dz
            
        fum[0:n+1]=np.exp(-zeta[0:n+1])

        Xm[0]=0.0
        for i in range(1,n+1):
        	Xm[i]=Xm[i-1]+(fum[i-1]+fum[i])*dz/2.0

        cb=1/Xm[n]
        cm[0:n+1]=cb*fum[0:n+1]
        c[1:n+1]=(cm[0:n]+cm[1:n+1])/2.0


        #---------  error	-----------------

        err_u[1:n+1]=abs(uo[1:n+1]-u[1:n+1])
        err_c[1:n+1]=abs(co[1:n+1]-c[1:n+1])
        maxerr=max(np.max(err_u[1:n+1]),np.max(err_c[1:n+1]))
        print(maxerr,cm[0])

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(elapsed_time)

    #---------ファイルに出力---------------------------------------------------------------
    # outfile = open('output.csv','w', newline='')
    # writer = csv.writer(outfile)
    # writer.writerow(['um[i]', 'zm[i]'])

    # for i in range(0,n+1):
    #     writer.writerow([um[i], zm[i]])

    # outfile.close()

    
    #----------グラフ表示------------------------------------------------------------------
    plt.title('um-zm')
    plt.xlabel('um[cm]')
    plt.ylabel('zm[cm]')
    plt.plot(um,zm,'b:')
    plt.show()

    # plt.title('cm-zm')
    # plt.xlabel('cm[cm]')
    # plt.ylabel('zm[cm]')
    # plt.xlim(0,0.05)
    # plt.plot(cm,zm,'k:')
    # plt.show()
    # 
    


if __name__ == "__main__":
    main()