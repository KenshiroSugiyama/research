# coding: UTF-8
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

def main():
    Rps = []
    Rfs = []
    for i in range(1,101):
        Rp = i
        a1=2.891394478
        a2=0.95296
        a3=0.056834671
        a4=0.0002892046
        a5=0.000244647

        Rf=np.exp(-a1+a2*np.log(Rp)-a3*(np.log(Rp))**2.0+a4*(np.log(Rp))**3.0-a5*(np.log(Rp))**4.0) #Dietrich falling velocity
        
        # print(Rp, end=' ')
        # print(Rf)
        Rps.append(Rp)
        print(Rps)
        Rfs.append(Rf)
        #グラフ出力
        # plt.title('E-vs')
        
        # color=["b:","g:","r:","c","y:","k:"]
    
    
    plt.xlabel('Rp')
    plt.ylabel('Rf')
    plt.xscale("log")
    # plt.yscale("log")
    plt.xlim(0,100)
    # pos = [1, 10, 100] 
    # ticks = ['1', '10']
    # plt.set_xticks(pos)
    # plt.set_xticklabels(ticks)
    plt.ylim(0,1.5)
    plt.plot(Rps,Rfs,"k:", linestyle = 'solid')
    # plt.legend()
    
    plt.show() # プロットを表示
    

if __name__ == "__main__":
    main()