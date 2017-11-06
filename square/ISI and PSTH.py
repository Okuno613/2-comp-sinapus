# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 17:47:44 2017

@author: tarol
"""
import math, numpy, scipy.optimize
import matplotlib.pyplot as plt



def main():

    count=0
    n=0
    m=0
    i=0
    l=0
    j=0
    k=0
    tmps=0
    size=60000
    time=3600
    cutsize=20

    numbers1 =numpy.zeros(size)
    numbers2 =numpy.zeros(size)
    numbers3 =numpy.zeros(size)
    isi =numpy.zeros(64)
    x_isi =numpy.arange(64)
    psth =numpy.zeros(time)
    x_psth =numpy.arange(time)
    x_psth_plot =numpy.arange(time/int(cutsize))
    psth_sum = numpy.zeros(time/int(cutsize))
    
    for line in open(OpenfileName+'.txt', 'r'):
        items = line.split()
        if( cutsize<count and count<size+cutsize+1 ):
            numbers1[n] = int(items[0])
            numbers2[n] = int(items[1])
            numbers3[n] = int(items[2])
            n=n+1
        count=count+1

    argmax1=numpy.argmax(numbers1)
    argmax2=numpy.argmax(numbers2)
    argmax3=numpy.argmax(numbers3)
    print argmax1
    print argmax2
    print argmax3


    for i in range(time):
        for j in range(size):
            if numbers2[j]==i and 5<i:
                psth[i]=psth[i]+1
                tmps=tmps+psth[i]
        if i%cutsize==0:
            psth_sum[m]=tmps
            m=m+1
            tmps=0
            
                
    for i in range(time):
        for j in range(64):
            if numbers3[i]==j and 0<j:
                isi[j]=isi[j]+1                
        
    for i in x_psth_plot:
        x_psth_plot[i]=x_psth_plot[i]*cutsize
    
    print psth_sum
        
                
    #plt.bar(x_psth,psth)  # plot true curve
    plt.bar(x_psth_plot,psth_sum,width=cutsize)
    plt.title('Least-squares fit to noisy data')
    plt.rcParams['font.family'] = 'IPAGothic'
    plt.title("PSTH" )
    plt.xlabel('Time[ms]')
    plt.ylabel('#Spikes')
    #plt.ylim(-0.01, 0.04)
    plt.savefig('PSTH '+str(OpenfileName)+'.png', dpi=150)
    plt.show()
    
    plt.bar(x_isi,isi)  # plot true curve
    plt.title('Least-squares fit to noisy data')
    plt.rcParams['font.family'] = 'IPAGothic'
    plt.title("ISI" )
    plt.xlabel('Time[ms]')
    plt.ylabel('Interspike interval[ms]')
    #plt.ylim(-0.01, 0.04)
    plt.savefig('ISI '+str(OpenfileName)+'.png', dpi=150)
    plt.show()    



OpenfileName='V_moved'
main()
OpenfileName='Vs_moved'
main()
