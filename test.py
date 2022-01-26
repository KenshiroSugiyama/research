# coding: UTF-8
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

def main():
    S=0.1
    Fr=np.sqrt(S)
    # Fr=S**2
    Rp_s=[1,10,100]
    for index, Rp in enumerate(Rp_s):
        print(Rp)
        E_all=[]
        vs_all=[]
        for k in range(1,101):
            vs=0.01*k	
            # print(vs)
            vs_all.append(vs)
            E=(5.73*10**(-3))*((1.0/vs)**1.31)*(Fr**1.59)*(Rp**(-0.86))
            E_all.append(E)
            # print(E)
        
        
        plt.xlabel('vs')
        plt.ylabel('E')
        # plt.xlim(0.01,1)
        # plt.ylim(0.0001,1)
        plt.xscale('log')
        plt.yscale('log')
        color=["b:","g:","r:","c","y:","k:"]
        label=["Rp=1","Rp=10","Rp=100"]
        plt.plot(vs_all,E_all,color[index], label=label[index])
        plt.legend()
    plt.show()

    

if __name__ == "__main__":
    main()