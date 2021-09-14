# coding: UTF-8
import math

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
    zm=[0]*n
    for i in range(0,n):
        zm[i]=zb+dz*i
        zm.append(zm[i])

    print(zm[n])

    

if __name__ == '__main__':
    main()

    
