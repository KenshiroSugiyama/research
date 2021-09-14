# coding: UTF-8
import math
import numpy as np

# def tridag(a,b,c,r,u,n):
#     NMAX=50000
#     if (b(1) == 0):
#         exit()
#         print('tridag: rewrite equations')
#     bet=b(1)
#     u(1)=r(1)/bet
#     for j in range(2, n):
#         gam(j)=c(j-1)/bet
#         bet=b(j)-a(j)*gam(j)
#         if (bet == 0): 
#             exit() 
#             print('tridag failed')
#         u(j)=(r(j)-a(j)*u(j-1))/bet

#     for j in range(n-1, 1, -1):
#         u(j)=u(j)-gam(j+1)*u(j+1)
    
def main():
    # a=[0]*4
    # b=[1]*4
    # c=[0]*4
    # c[0:3]=a[0:3]-3
    # print(c)
    kappa=0.4
    sgm_c=1.3

    n=5000

    Ri=0.01
    ds=1.0/(1.0*10**3)
    ks=2.5*ds
    zb=1.0*ds
    ub=1.0/kappa*math.log(30.0*zb/ks)

    dt=1.0*10**(-7)
    dz=(1.0-zb)/n
    zm=np.empty(n+1)
    # print(zm)
    for i in range(0,n):
        zm[i]=zb+dz*i
    
    z=np.empty(n+1)
    for j in range(1,n):
        z[j]=zm[j]-0.5*dz
    # print(zm)

    #--------- Initial velocity	-----------------------
    u=np.empty(n+1)
    for j in range(1,n):
        u[j]=1.0/kappa*math.log(30.0*z[j]/ks)
    
    c=np.empty(n+1)
    for j in range(1,n):
        c[j]=1.0-z[j]

    uo=np.empty(n+1)
    co=np.empty(n+1)
    for j in range(0,50):
        vs=0.001*500**(0.02*(50-j))
        maxerr=1.0
        # while maxerr>1*10**(-6):
        #     uo[j]=u[j]
        #     co[j]=c[j]

        #------------  gradient	-------------
        dudzm=np.empty(n+1)
        dudzm[0]=2.0*(uo[1]-ub)/dz
        for j in range(1,n-1):
            dudzm[j]=(uo[j+1]-uo[j])/dz
        dudzm[n]=0.0

        dudz=np.empty(n+1)
        for j in range(1,n):
            dudz[j]=(dudzm[j-1]+dudzm[j])/2.0

        dcdzm=np.empty(n+1)
        dcdzm[0]=2.0*(co[1]-cb)/dz
        for j in range(1,n-1):
            dcdzm[j]=(co[j+1]-co[j])/dz
        dcdzm[n]=2.0*(0-co[n])/dz
        
        dcdz=np.empty(n+1)
        for j in range(1,n):
            dcdz[j]=(dcdzm[j-1]+dcdzm[j])/2.0

        print('終了')

		# #-------	eddy viscosity	----------------

		# cf(0:n)=(dudzm(0:n)**2-1.350*Ri*dcdzm(0:n))&
		# 			 /(dudzm(0:n)**2-14.85*Ri*dcdzm(0:n))
		# nu_tm(0:n)=kappa**2*zm(0:n)**2*dudzm(0:n)*cf(0:n)
		# nu_t(1:n)=(nu_tm(0:n-1)+nu_tm(1:n))/2.0
		# nu_tcm(0:n)=nu_tm(0:n)/sgm_c

		# #------- 	velocity	----------------------

		# D_u(0:n)=nu_tm(0:n)/dz

		# a_u(1)=0.0
		# a_u(2:n)=-D_u(1:n-1)
		# b_u(1)=dz/dt+2.0*D_u(0)+D_u(1)
		# b_u(2:n-1)=dz/dt+D_u(1:n-2)+D_u(2:n-1)
		# b_u(n)=dz/dt+D_u(n-1)
		# c_u(1:n-1)=-D_u(1:n-1)
		# c_u(n)=0.0

		# r_u(1)=(uo(1)/dt+c(1))*dz+2.0*D_u(0)*ub
		# r_u(2:n)=(uo(2:n)/dt+c(2:n))*dz

		# call tridag(a_u,b_u,c_u,r_u,u,n)

		# #-------------- concentration --------------------------

		# zeta(0)=0
		# for i in range(1,n):
		# 	zeta[i]=zeta[i-1]+vs/nu_t[i]*dz
		# 	write(*,*) i, zeta(i), nu_t(i)

		# fum[0:n]=math.exp(-zeta[0:n])

		# Xm(0)=0.0
		# for i in range(1,n):
		# 	Xm[i]=Xm(i-1)+(fum(i-1)+fum(i))*dz/2.0
		# 	write(*,*) i, fum(i)
		

		# cb=1/Xm(n)
		# cm(0:n)=cb*fum(0:n)
		# c(1:n)=(cm(0:n-1)+cm(1:n))/2.0

		# #---------  error	-----------------

		# err_u(1:n)=dabs(uo(1:n)-u(1:n))
		# err_c(1:n)=dabs(co(1:n)-c(1:n))
		# maxerr=max(maxval(err_u(1:n)),maxval(err_c(1:n)))



if __name__ == '__main__':
    main()

    
