# coding: UTF-8
import math
import numpy as np
import warnings
warnings.simplefilter('ignore', category=RuntimeWarning)  # RuntimeWarningを無視扱いに設定
import re
import csv
import time
import matplotlib.pyplot as plt


def main():
    #配列作成
    n=100000
    zm=np.empty(n+1) 
    z=np.empty(n+1) 
    u=np.empty(n+1) 
    c=np.empty(n+1) 
    uo=np.empty(n+1) 
    co=np.empty(n+1)  
    dudzm=np.empty(n+1) 
    dudz=np.empty(n+1) 
    dcdzm=np.empty(n+1) 
    dcdz=np.empty(n+1) 
    cf=np.empty(n+1)
    nu_tm=np.empty(n+1)
    nu_t=np.empty(n+1)
    nu_tcm=np.empty(n+1)   
    D_u=np.empty(n+1)
    a_u=np.empty(n+1)
    b_u=np.empty(n+1)
    c_u=np.empty(n+1)
    r_u=np.empty(n+1)
    um=np.empty(n+1)
    zeta=np.empty(n+1)
    fum=np.empty(n+1)
    Xm=np.empty(n+1)
    X=np.empty(n+1)
    cm=np.empty(n+1)
    err_u=np.empty(n+1)
    err_c=np.empty(n+1)
    root1=np.empty(n+1)
    root2=np.empty(n+1)

    kappa = 0.4
    sgm_c = 1.3

    Rp=1
    S=0.1
    ds=1.0*10**(-4)
    ks=2.5*ds
    Rit=1.0/S
    zb=1.0*ds/12.0

    Fr=S**2 #???
    

    a1=2.891394478
    a2=0.95296
    a3=0.056834671
    a4=0.0002892046
    a5=0.000244647

    Rf=np.exp(-a1+a2*np.log(Rp)-a3*(np.log(Rp))**2.0-a4*(np.log(Rp))**3.0+a5*(np.log(Rp))**4.0) #Dietrich falling velocity
    
    dt=1.0*10**(-9)
    dz=(1.0-zb)/n
    for i in range(0,n+1):
        zm[i]=zb+dz*i
	
    z[1:n+1]=zm[1:n+1]-0.5*dz
    # print(z[n])    

#!-------- Initial velocity	-----------------------
    for k in range(0,91):
        vs=0.01+0.001*k	

        za=zm[100] #value of H=0.001+ds/12
        ips=za-zb #thichness of bottom layer
        h=za/zb 

        u[1:n+1]=1.0/kappa*np.log(30.0*z[1:n+1]/ks)
        # print(u[n])
        c[1:n+1]=1.0-z[1:n+1]-zb
        X[1:n+1]=z[1:n+1]
        
        E=(5.73/10**3)*((1.0/vs)**1.31)*(Fr**1.59)*(Rp**(-0.86))
        # print(E)
        alpha=-(vs/kappa)
        ca=1.6
        cb=ca*((h)**(-alpha))
        ub=1.0/kappa*np.log(30.0*zm[0]/ks)
        Xm[1:n]=(X[1:n]+X[2:n+1])/2
        Xm[n]=1
        Xm[0]=0

        maxerr=1.0
        while maxerr>1*10**(-6):
            #------------  gradient	-------------

            uo[1:n+1]=u[1:n+1]
            # print(uo[n])
            um[0]=ub
            um[1:n]=(uo[1:n]+uo[2:n+1])/2
            um[n]=uo[n]

            co[1:n+1]=c[1:n+1]
            cm[0]=cb
            cm[1:n]=(co[1:n]+co[2:n+1])/2
            cm[n]=0


            dudzm[0]=2*(uo[1]-um[0])/dz
            dudzm[1:n]=(uo[2:n+1]-uo[1:n])/dz
            dudzm[n]=2*(um[n]-uo[n])/dz
            # print(dudzm[n-3:n+1])
            dudz[1:n+1]=(um[1:n+1]-um[0:n])/dz


            dcdzm[0]=2*(co[1]-cm[0])/dz
            dcdzm[1:n]=(co[2:n+1]-co[1:n])/dz
            dcdzm[n]=2*(cm[n]-co[n])/dz
            dcdz[1:n+1]=(cm[1:n+1]-cm[0:n])/dz
        
            #-------	eddy viscosity	----------------

            cof=(dudzm**2-1.35*Rit*dcdzm)/(dudzm**2-14.85*Rit*dcdzm)
            nu_tm=cof*(kappa**2)*((zm)**2)*(1-zm)*abs(dudzm)   
            #nu_tm(0:n)=kappa**2*zm(0:n)**2*(1d0-zm(0:n))*dabs(dudzm(0:n))

            nu_t[1:n+1]=(nu_tm[0:n]+nu_tm[1:n+1])/2

            #-------    u     ----------------------
        
            #bottom layer equation

            for j in range(0,101):
                um[j]=(1.0/kappa)*np.log(zm[j]/zb)+ips*((zm[j]-zb)/(2.0*kappa)+zb*(np.log(zm[j]/zb))*(cb-kappa+vs)*(1.0/(kappa-vs))-ips*kappa*cb*zb*((zm[j]/zb)**(1.0+alpha)-1.0)*(0.5/(kappa-vs)**2.0))
            # print(um[100])
            #upper layer equation

            #coefficient of dudz equation
            ak=(kappa**2)*((zm)**2)*(1-zm)
            bk=-1.35*Rit*(kappa**2)*((zm)**2)*(1-zm)*dcdzm-(1-Xm)
            ck=14.85*Rit*dcdzm*(1-Xm)

            root1=(-bk+(bk**2.0-4*ak*ck)**0.5)/(2*ak)

            root2=root1**0.5  #solve Quartic equation
            
            dudz[101:n+1]=(dudzm[100:n]+dudzm[101:n+1])/2
            # print(dudz[n-3:n+1])
            for i in range(101,n+1):
                um[i]=um[i-1]+dudz[i]*dz

            u[1:101]=(um[0:100]+um[1:101])/2.0
            # print(u[n-2:n+1])
            #-------------- C,X --------------------------
            for j in range(0,101):
                cm[j]=cb*(zm[j]/zb)**alpha+ips*vs*cb/(2.0*kappa)*(3*(zm[j]-zb)+zb*(kappa*cb/(kappa-vs)-1.0)*np.log(zm[j]/zb)-kappa**2*cb*zb/(kappa-vs)**2*((zm[j]/zb)**(1.0+alpha)-1.0))*(zm[j]/zb)**alpha
            

            Xm[0]=0.0
            for j in range(1,101):
                Xm[j]=Xm[j-1]+(cm[j]+cm[j-1])/2.0*dz
            # print(Xm[100])
            zeta[100]=0
            for i in range(101,n+1):
                zeta[i]=zeta[i-1]+vs/nu_t[i]*dz
        # !			write(*,*) i, zeta(i), nu_t(i)
            # print(zeta[n-3:n+1])
            fum[100:n+1]=np.exp(-zeta[100:n+1])
            # print(fum[n-3:n+1])
            cm[100:n+1]=cm[100]*fum[100:n+1]
            # print(cm[n-3:n+1])
            for i in range(101,n+1):
                Xm[i]=Xm[i-1]+(cm[i-1]+cm[i])*dz/2.0
        # !			write(*,*) i, fum(i)
            


            cb=cb/Xm[n]
            for j in range(0,101):
                cm[j]=cb*(zm[j]/zb)**alpha+ips*vs*cb/(2.0*kappa)*(3*(zm[j]-zb)+zb*(kappa*cb/(kappa-vs)-1.0)*np.log(zm[j]/zb)-kappa**2*cb*zb/(kappa-vs)**2*((zm[j]/zb)**(1.0+alpha)-1.0))*(zm[j]/zb)**alpha
            
            cm[100:n+1]=cm[100]*fum[100:n+1]
            c[1:n+1]=(cm[0:n]+cm[1:n+1])/2.0
            Xm[0]=0.0
            for i in range(1,n+1):
                Xm[i]=Xm[i-1]+(cm[i-1]+cm[i])*0.5*dz

            
            X[1:n+1]=(Xm[0:n]+Xm[1:n+1])/2
            # print(X[n-2:n+1])
            #!---------  error	-----------------
            err_u=abs(uo-u)
            # print(co[0:3], end=' ')
            # print(c[0:3])
            # print(err_u[0:2])
            err_c=abs(co-c)
            # print(err_c[0:2])
            maxerr=max(np.max(err_u),np.max(err_c))

            # write(*,*)'maxerr=', maxerr,'k=', k
            # print(rf'$maxerr={{{maxerr}}}$')
            # print(rf'$k={{{k}}}$')
        # print(E)
        Hds=cm[0]*Rf**2/(E*S*vs**2)
        # print("H==")
        print(Hds, end=' ')
        print(Xm[n])
        # write(*,*)'H=',HDs,'X_all=',Xm(n)
        # !write(50,*)vs, cm(0)
        # write(60,*)vs, Hds
    #----------グラフ表示------------------------------------------------------------------
    plt.title('E-vs')
    plt.xlabel('vs')
    plt.ylabel('E')
    plt.xlim(0,1)
    # color=["b:","g:","r:","c","y:","k:"]
    plt.plot(cm,zm,"b:")
    # plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
