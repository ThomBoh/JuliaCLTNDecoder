import math
import numpy as np
import scipy.linalg as la
import quimb as qu
import quimb.tensor as qtn
import time
import matplotlib.pyplot as plt
import random
from decimal import *

def get_binall(varray):
    """Convert an array of uint16 values to binary"""
    return np.unpackbits(varray.reshape(varray.shape[0], 1).view(np.uint8), axis=1)

def get_bin(x, n=0):
    y=format(x, 'b').zfill(n)
    l=len(y)
    r=[]
    for i in range(0,l):
        r.append(int(ord(y[i])-48))
    return r

def get_int(x):
    l=len(x)
    y=0
    for i in range(0,l):
        b=x[i]
        y=y+b*2**(l-i-1)
    return int(y)

#helper functions that build probability distributions into local tensors

#for single qubit location

def sqgz(z,psq):
    zm2=z%2
    out = psq[zm2]
    return out

#for two-qubit location

def tqgz(zc,zt,ptq):
    ret=0
    zcm2=zc%2
    ztm2=zt%2
    if (zcm2==0) & (ztm2==0):
        ret=ptq[0]
    elif (zcm2==1) & (ztm2==0):
        ret=ptq[1]
    elif (zcm2==0) & (ztm2==1):
        ret=ptq[2]
    else:
        ret=ptq[3]
    return ret


#create local tensor modeling CNOT location

def RepCNOTTensorHW(ptq):
    TOut=np.zeros((2,2,2,2))
    for itera in range(0,2**4):
        [zcl,zcr,ztl,ztr]=get_bin(itera,4)
        TOut[ztl][zcl][ztr][zcr]=tqgz(zcl+zcr+ztl,ztl+ztr,ptq)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling a simple wait location

def IdTensor(psq):
    TOut=np.zeros((2,2))
    for itera in range(0,2**2):
        [zl,zr]=get_bin(itera,2)
        TOut[zl][zr]=sqgz(zl+zr,psq)
    TOut=TOut.astype('float128')
    return TOut

#local tensor modeling terminal data qubit wait location. Two extra indices, one for applying logical operator, and one for applying a pure error

def DatTermTensor(psq):
    TOut=np.zeros((2,2,2))
    for itera in range(0,2**3):
        [zl,zlog,zpur]=get_bin(itera,3)
        TOut[zl][zlog][zpur]=sqgz(zl+zlog+zpur,psq)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling a first data qubit wait location

def DTermWait(psq):
    TOut=np.zeros(2)
    for itera in range(0,2):
        TOut[itera]=sqgz(itera,psq)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling an ancilla qubit measurement location

def AMeasTensor1(pm):
    TOut=np.zeros(2)
    for itera in range(0,2):
        TOut[itera]=sqgz(itera,pm)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling an ancilla qubit measurement location, with extra leg for injecting a syndrome bit flip

def AMeasTensor2(pm):
    TOut=np.zeros((2,2))
    for itera in range(0,2**2):
        [zl,zbf]=get_bin(itera,2)
        TOut[zl][zbf]=sqgz(zl+zbf,pm)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling ancilla qubit preparation (which has extra leg to inject a syndrome bit flip)

def APrepTensor1(pp):
    TOut= np.zeros((2,2))
    for itera in range(0,2**2):
        [zl,zr]=get_bin(itera,2)
        TOut[zl][zr]=sqgz(zl+zr,pp)
    TOut=TOut.astype('float128')
    return TOut

#create local tensor modeling ancilla qubit preparation

def APrepTensor2(pp):
    TOut= np.zeros(2)
    for itera in range(0,2):
        TOut[itera]=sqgz(itera,pp)
    TOut=TOut.astype('float128')
    return TOut

#The following function produces a circuit level tensor network modeling a measurement circuit for the repetition code of size/distance nq, where the syndrome measurement is repeated nr times. p is the noise parameter k1/k2.

def CLTNRepHWFull(nq,nr,p):
    
    #First initialize all HW noise parameters. Note that p=k1/k2
    
    pxm=0.004 #X readout failure probability
    kappa2=10**7
    alpha2=8 #photon number
    p1=10*p #fail prob for data qubit wait location during very first ancilla initialization
    p2=0.299*(p**0.5) #fail prob for data qubit wait location during execution of a CNOT on other qubits
    p3=p*kappa2*alpha2*(350*(10**(-9)) +(10/(kappa2*alpha2))) #fail prob for data qubit wait location during readout of ancilla
    pz1=0.845*(p**0.5) #prob of z1 otimes I failure after CNOT (qubit 1 is control, 2 is target)
    pz2=0.133*(p**0.5) #prob of I otimes z2 failure after CNOT (qubit 1 is control, 2 is target)
    pz1z2=0.133*(p**0.5) #prob of z1 otimes z2 failure after CNOT
    
    #create probability distributions to be passed to local tensor constructors which model different locations
    pi1=[1-p1,p1] #data qubit wait during first ancilla init
    pi2=[1-p2,p2] #data qubit wait during CNOT
    pi3=[1-p3,p3] #data qubit wait during ancilla measurement and re-init
    pm=[1-pxm,pxm] #X measurement location
    pp=[1-15*p/2,15*p/2] #ancilla initialization
    perf=[1,0] #perfect wait location
    ptq=[1-pz1-pz2-pz1z2,pz1,pz2,pz1z2] #CNOT location
    
    #Create a GHZ tensor that will be used to actuate application of logical operators
    gh0=format(0, 'b').zfill(nq+1)
    gh1=format(2**(nq+1)-1, 'b').zfill(nq+1)
    ghz = (qtn.MPS_computational_state(gh0) + qtn.MPS_computational_state(gh1))
    ghz=ghz.contract(all)
    ghzd=ghz.data
    
    #create the tensor network
    for r in range(0,nr):
        
        #create initalization/first wait locations
        if r==0:
            tagsD=('D{:d}T0'.format(nq-1))
            indsD=('d{:d}g0'.format(nq-1),)
            TN=qtn.Tensor(DTermWait(pi1),indsD,tagsD)
            for dq in range(0,nq-1):
                tagsD=('D{:d}T0'.format(dq))
                indsD=('d{:d}g0'.format(dq),)
                TN=TN&qtn.Tensor(DTermWait(pi1),indsD,tagsD)
        else:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r))
                indsD=('d{:d}g{:d}'.format(dq,4*r-1),'d{:d}g{:d}'.format(dq,4*r))
                TN=TN&qtn.Tensor(IdTensor(perf),indsD,tagsD)
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r))
            indsA=('sf{:d}r{:d}'.format(aq,r),'a{:d}g{:d}'.format(aq,4*r))
            TN=TN&qtn.Tensor(APrepTensor1(pp),indsA,tagsA)
        
        #first round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq,aq,4*r+1))
            indsD=('d{:d}g{:d}'.format(aq,4*r),'a{:d}g{:d}'.format(aq,4*r),'d{:d}g{:d}'.format(aq,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(nq-1,4*r+1))
        indsD=('d{:d}g{:d}'.format(nq-1,4*r),'d{:d}g{:d}'.format(nq-1,4*r+1))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        
        #second round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq+1,aq,4*r+2))
            indsD=('d{:d}g{:d}'.format(aq+1,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1),'d{:d}g{:d}'.format(aq+1,4*r+2),'a{:d}g{:d}'.format(aq,4*r+2))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(0,4*r+2))
        indsD=('d{:d}g{:d}'.format(0,4*r+1),'d{:d}g{:d}'.format(0,4*r+2))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        
        #measurement/wait locations
        #first do the data wait locations
        if r==nr-1:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}lz'.format((nq-1)//2))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq-1))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
        else:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}g{:d}'.format((nq-1)//2,4*r+3))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
        
        #now do the ancilla measurement failure locations
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r+3))
            indsA=('a{:d}g{:d}'.format(aq,4*r+2),)
            TN=TN&qtn.Tensor(AMeasTensor1(pm),indsA,tagsA)
    
    #finally, put in the logical control tensor
    tagsL=('LZ')
    indsL=('lz',)
    for dq in range(0,nq):
        indsL=indsL+('d{:d}lz'.format(dq),)
    TN=TN&qtn.Tensor(ghzd,indsL,tagsL)
    
    return TN

#The following function produces a circuit level tensor network as above, but traces over the pure error degrees of freedom

def CLTNRepHWPETrace(nq,nr,p):
    px=0.004
    kappa2=10**7
    alpha2=8
    pi1=[1-10*p,10*p]
    pi2=[1-0.299*(p**0.5),0.299*(p**0.5)]
    p3=p*kappa2*alpha2*(350*(10**(-9)) +(10/(kappa2*alpha2)))
    pi3=[1-p3,p3]
    pz1=0.845*(p**0.5)
    pz2=0.133*(p**0.5)
    pz1z2=0.133*(p**0.5)
    pm=[1-px,px]
    pp=[1-15*p/2,15*p/2]
    perf=[1,0]
    ptq=[1-pz1-pz2-pz1z2,pz1,pz2]
    #print(psq)
    #print(ptq)
    #print(pm)
    #print(pp)
    gh0=format(0, 'b').zfill(nq+1)
    gh1=format(2**(nq+1)-1, 'b').zfill(nq+1)
    ghz = (qtn.MPS_computational_state(gh0) + qtn.MPS_computational_state(gh1))
    ghz=ghz.contract(all)
    ghzd=ghz.data
    for r in range(0,nr):
        if r==0:
            tagsD=('D{:d}T0'.format(nq-1))
            indsD=('d{:d}g0'.format(nq-1),)
            TN=qtn.Tensor(DTermWait(pi1),indsD,tagsD)
            for dq in range(0,nq-1):
                tagsD=('D{:d}T0'.format(dq))
                indsD=('d{:d}g0'.format(dq),)
                TN=TN&qtn.Tensor(DTermWait(pi1),indsD,tagsD)
        else:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r))
                indsD=('d{:d}g{:d}'.format(dq,4*r-1),'d{:d}g{:d}'.format(dq,4*r))
                TN=TN&qtn.Tensor(IdTensor(perf),indsD,tagsD)
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r))
            indsA=('sf{:d}r{:d}'.format(aq,r),'a{:d}g{:d}'.format(aq,4*r))
            TN=TN&qtn.Tensor(APrepTensor(pp),indsA,tagsA)
        #first round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq,aq,4*r+1))
            indsD=('d{:d}g{:d}'.format(aq,4*r),'a{:d}g{:d}'.format(aq,4*r),'d{:d}g{:d}'.format(aq,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(nq-1,4*r+1))
        indsD=('d{:d}g{:d}'.format(nq-1,4*r),'d{:d}g{:d}'.format(nq-1,4*r+1))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #second round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq+1,aq,4*r+2))
            indsD=('d{:d}g{:d}'.format(aq+1,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1),'d{:d}g{:d}'.format(aq+1,4*r+2),'a{:d}g{:d}'.format(aq,4*r+2))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(0,4*r+2))
        indsD=('d{:d}g{:d}'.format(0,4*r+1),'d{:d}g{:d}'.format(0,4*r+2))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #measurement round/wait locations
        #first do the terminal data wait locations
        if r==nr-1:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
                tagsT=('D{:d}P'.format(dq))
                indsT=('d{:d}p'.format(dq),)
                TN=TN&qtn.Tensor(np.ones(2),indsT,tagsT)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}lz'.format((nq-1)//2))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq-1))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
                tagsT=('D{:d}P'.format(dq-1))
                indsT=('d{:d}p'.format(dq-1),)
                TN=TN&qtn.Tensor(np.ones(2),indsT,tagsT)
        else:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}g{:d}'.format((nq-1)//2,4*r+3))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
        #now do the ancilla measurement failure locations
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r+3))
            indsA=('a{:d}g{:d}'.format(aq,4*r+2),)
            TN=TN&qtn.Tensor(AMeasTensor1(pm),indsA,tagsA)
    #finally, put in the logical control tensor
    tagsL=('LZ')
    indsL=('lz',)
    for dq in range(0,nq):
        indsL=indsL+('d{:d}lz'.format(dq),)
    TN=TN&qtn.Tensor(ghzd,indsL,tagsL)
    return TN

def CLTNRepHWTraceMeas(nq,nr,p):
        #First initialize all HW noise parameters. Note that p=k1/k2
    
    pxm=0.004 #X readout failure probability
    kappa2=10**7
    alpha2=8 #photon number
    p1=10*p #fail prob for data qubit wait location during very first ancilla initialization
    p2=0.299*(p**0.5) #fail prob for data qubit wait location during execution of a CNOT on other qubits
    p3=p*kappa2*alpha2*(350*(10**(-9)) +(10/(kappa2*alpha2))) #fail prob for data qubit wait location during readout of ancilla
    pz1=0.845*(p**0.5) #prob of z1 otimes I failure after CNOT (qubit 1 is control, 2 is target)
    pz2=0.133*(p**0.5) #prob of I otimes z2 failure after CNOT (qubit 1 is control, 2 is target)
    pz1z2=0.133*(p**0.5) #prob of z1 otimes z2 failure after CNOT
    
    #create probability distributions to be passed to local tensor constructors which model different locations
    pi1=[1-p1,p1] #data qubit wait during first ancilla init
    pi2=[1-p2,p2] #data qubit wait during CNOT
    pi3=[1-p3,p3] #data qubit wait during ancilla measurement and re-init
    pm=[1-pxm,pxm] #X measurement location
    pp=[1-15*p/2,15*p/2] #ancilla initialization
    perf=[1,0] #perfect wait location
    ptq=[1-pz1-pz2-pz1z2,pz1,pz2,pz1z2] #CNOT location
    
    gh0=format(0, 'b').zfill(nq+1)
    gh1=format(2**(nq+1)-1, 'b').zfill(nq+1)
    ghz = (qtn.MPS_computational_state(gh0) + qtn.MPS_computational_state(gh1))
    ghz=ghz.contract(all)
    ghzd=ghz.data
    for r in range(0,nr):
        if r==0:
            tagsD=('D{:d}T0'.format(nq-1))
            indsD=('d{:d}g0'.format(nq-1),)
            TN=qtn.Tensor(DTermWait(pi1),indsD,tagsD)
            for dq in range(0,nq-1):
                tagsD=('D{:d}T0'.format(dq))
                indsD=('d{:d}g0'.format(dq),)
                TN=TN&qtn.Tensor(DTermWait(pi1),indsD,tagsD)
        else:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r))
                indsD=('d{:d}g{:d}'.format(dq,4*r-1),'d{:d}g{:d}'.format(dq,4*r))
                TN=TN&qtn.Tensor(IdTensor(perf),indsD,tagsD)
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r))
            indsA=('sf{:d}r{:d}'.format(aq,r),'a{:d}g{:d}'.format(aq,4*r))
            TN=TN&qtn.Tensor(APrepTensor1(pp),indsA,tagsA)
            tagsST=('A{:d}T{:d}TR'.format(aq,4*r))
            indsST=('sf{:d}r{:d}'.format(aq,r),)
            TN=TN&qtn.Tensor(np.ones(2),indsST,tagsST)
        #first round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq,aq,4*r+1))
            indsD=('d{:d}g{:d}'.format(aq,4*r),'a{:d}g{:d}'.format(aq,4*r),'d{:d}g{:d}'.format(aq,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(nq-1,4*r+1))
        indsD=('d{:d}g{:d}'.format(nq-1,4*r),'d{:d}g{:d}'.format(nq-1,4*r+1))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #second round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq+1,aq,4*r+2))
            indsD=('d{:d}g{:d}'.format(aq+1,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1),'d{:d}g{:d}'.format(aq+1,4*r+2),'a{:d}g{:d}'.format(aq,4*r+2))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(0,4*r+2))
        indsD=('d{:d}g{:d}'.format(0,4*r+1),'d{:d}g{:d}'.format(0,4*r+2))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #measurement round/wait locations
        #first do the terminal data wait locations
        if r==nr-1:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}lz'.format((nq-1)//2))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq-1))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
        else:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}g{:d}'.format((nq-1)//2,4*r+3))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
        #now do the ancilla measurement failure locations
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r+3))
            indsA=('a{:d}g{:d}'.format(aq,4*r+2),)
            TN=TN&qtn.Tensor(AMeasTensor1(pm),indsA,tagsA)
    #finally, put in the logical control tensor
    tagsL=('LZ')
    indsL=('lz',)
    for dq in range(0,nq):
        indsL=indsL+('d{:d}lz'.format(dq),)
    TN=TN&qtn.Tensor(ghzd,indsL,tagsL)
    return TN

#The following function produces a circuit level tensor network as above, but traces over all degrees of freedom. Contracting it should produce 1, showing that the tensor network gives a genuine probability distribution.

def CLTNRepHWTraceAll(nq,nr,p):
    px=0.004
    kappa2=10**7
    alpha2=8
    pi1=[1-10*p,10*p]
    pi2=[1-0.299*(p**0.5),0.299*(p**0.5)]
    p3=p*kappa2*alpha2*(350*(10**(-9)) +(10/(kappa2*alpha2)))
    pi3=[1-p3,p3]
    pz1=0.845*(p**0.5)
    pz2=0.133*(p**0.5)
    pz1z2=0.133*(p**0.5)
    pm=[1-px,px]
    pp=[1-15*p/2,15*p/2]
    perf=[1,0]
    ptq=[1-pz1-pz2-pz1z2,pz1,pz2]
    #print(psq)
    #print(ptq)
    #print(pm)
    #print(pp)
    gh0=format(0, 'b').zfill(nq+1)
    gh1=format(2**(nq+1)-1, 'b').zfill(nq+1)
    ghz = (qtn.MPS_computational_state(gh0) + qtn.MPS_computational_state(gh1))
    ghz=ghz.contract(all)
    ghzd=ghz.data
    for r in range(0,nr):
        if r==0:
            tagsD=('D{:d}T0'.format(nq-1))
            indsD=('d{:d}g0'.format(nq-1),)
            TN=qtn.Tensor(DTermWait(pi1),indsD,tagsD)
            for dq in range(0,nq-1):
                tagsD=('D{:d}T0'.format(dq))
                indsD=('d{:d}g0'.format(dq),)
                TN=TN&qtn.Tensor(DTermWait(pi1),indsD,tagsD)
        else:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r))
                indsD=('d{:d}g{:d}'.format(dq,4*r-1),'d{:d}g{:d}'.format(dq,4*r))
                TN=TN&qtn.Tensor(IdTensor(perf),indsD,tagsD)
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r))
            indsA=('sf{:d}r{:d}'.format(aq,r),'a{:d}g{:d}'.format(aq,4*r))
            TN=TN&qtn.Tensor(APrepTensor(pp),indsA,tagsA)
            tagsST=('A{:d}T{:d}TR'.format(aq,4*r))
            indsST=('sf{:d}r{:d}'.format(aq,r),)
            TN=TN&qtn.Tensor(np.ones(2),indsST,tagsST)
        #first round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq,aq,4*r+1))
            indsD=('d{:d}g{:d}'.format(aq,4*r),'a{:d}g{:d}'.format(aq,4*r),'d{:d}g{:d}'.format(aq,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(nq-1,4*r+1))
        indsD=('d{:d}g{:d}'.format(nq-1,4*r),'d{:d}g{:d}'.format(nq-1,4*r+1))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #second round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq+1,aq,4*r+2))
            indsD=('d{:d}g{:d}'.format(aq+1,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1),'d{:d}g{:d}'.format(aq+1,4*r+2),'a{:d}g{:d}'.format(aq,4*r+2))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(0,4*r+2))
        indsD=('d{:d}g{:d}'.format(0,4*r+1),'d{:d}g{:d}'.format(0,4*r+2))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #measurement round/wait locations
        #first do the terminal data wait locations
        if r==nr-1:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
                tagsT=('D{:d}P'.format(dq))
                indsT=('d{:d}p'.format(dq),)
                TN=TN&qtn.Tensor(np.ones(2),indsT,tagsT)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}lz'.format((nq-1)//2))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}lz'.format(dq),'d{:d}p'.format(dq-1))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
                tagsT=('D{:d}P'.format(dq-1))
                indsT=('d{:d}p'.format(dq-1),)
                TN=TN&qtn.Tensor(np.ones(2),indsT,tagsT)
        else:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}g{:d}'.format((nq-1)//2,4*r+3))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
        #now do the ancilla measurement failure locations
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r+3))
            indsA=('a{:d}g{:d}'.format(aq,4*r+2),)
            TN=TN&qtn.Tensor(AMeasTensor1(pm),indsA,tagsA)
    #finally, put in the logical control tensor
    tagsL=('LZ')
    indsL=('lz',)
    for dq in range(0,nq):
        indsL=indsL+('d{:d}lz'.format(dq),)
    TN=TN&qtn.Tensor(ghzd,indsL,tagsL)
    tagsLT=('LZT')
    indsLT=('lz',)
    TN=TN&qtn.Tensor(np.ones(2),indsLT,tagsLT)
    return TN

#The following function produces a circuit level tensor network with extra measurement failure locations on the data qubits to simulate additional faulty readout before ideal decoding.

def CLTNRepHWQC(nq,nr,p):
    px=0.004
    kappa2=10**7
    alpha2=8
    pi1=[1-10*p,10*p]
    pi2=[1-0.299*(p**0.5),0.299*(p**0.5)]
    p3=p*kappa2*alpha2*(350*(10**(-9)) +(10/(kappa2*alpha2)))
    pi3=[1-p3,p3]
    pz1=0.845*(p**0.5)
    pz2=0.133*(p**0.5)
    pz1z2=0.133*(p**0.5)
    pm=[1-px,px]
    pp=[1-15*p/2,15*p/2]
    perf=[1,0]
    ptq=[1-pz1-pz2-pz1z2,pz1,pz2,pz1z2]
    #print(psq)
    #print(ptq)
    #print(pm)
    #print(pp)
    gh0=format(0, 'b').zfill(nq+1)
    gh1=format(2**(nq+1)-1, 'b').zfill(nq+1)
    ghz = (qtn.MPS_computational_state(gh0) + qtn.MPS_computational_state(gh1))
    ghz=ghz.contract(all)
    ghzd=ghz.data
    for r in range(0,nr):
        if r==0:
            tagsD=('D{:d}T0'.format(nq-1))
            indsD=('d{:d}g0'.format(nq-1),)
            TN=qtn.Tensor(DTermWait(pi1),indsD,tagsD)
            for dq in range(0,nq-1):
                tagsD=('D{:d}T0'.format(dq))
                indsD=('d{:d}g0'.format(dq),)
                TN=TN&qtn.Tensor(DTermWait(pi1),indsD,tagsD)
        else:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r))
                indsD=('d{:d}g{:d}'.format(dq,4*r-1),'d{:d}g{:d}'.format(dq,4*r))
                TN=TN&qtn.Tensor(IdTensor(perf),indsD,tagsD)
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r))
            indsA=('sf{:d}r{:d}'.format(aq,r),'a{:d}g{:d}'.format(aq,4*r))
            TN=TN&qtn.Tensor(APrepTensor1(pp),indsA,tagsA)
        #first round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq,aq,4*r+1))
            indsD=('d{:d}g{:d}'.format(aq,4*r),'a{:d}g{:d}'.format(aq,4*r),'d{:d}g{:d}'.format(aq,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(nq-1,4*r+1))
        indsD=('d{:d}g{:d}'.format(nq-1,4*r),'d{:d}g{:d}'.format(nq-1,4*r+1))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #second round of CNOTs
        for aq in range(0,nq-1):
            tagsD=('D{:d}A{:d}T{:d}'.format(aq+1,aq,4*r+2))
            indsD=('d{:d}g{:d}'.format(aq+1,4*r+1),'a{:d}g{:d}'.format(aq,4*r+1),'d{:d}g{:d}'.format(aq+1,4*r+2),'a{:d}g{:d}'.format(aq,4*r+2))
            TN=TN&qtn.Tensor(RepCNOTTensorHW(ptq),indsD,tagsD)
        tagsD=('D{:d}T{:d}'.format(0,4*r+2))
        indsD=('d{:d}g{:d}'.format(0,4*r+1),'d{:d}g{:d}'.format(0,4*r+2))
        TN=TN&qtn.Tensor(IdTensor(pi2),indsD,tagsD)
        #measurement round/wait locations
        #first do the terminal data wait locations
        if r==nr-1:
            for dq in range(0,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pm),indsD,tagsD)
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+4))
                indsD=('d{:d}g{:d}'.format(dq,4*r+3),'d{:d}lz'.format(dq),'d{:d}p'.format(dq))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+4))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+3),'d{:d}lz'.format((nq-1)//2))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+4))
                indsD=('d{:d}g{:d}'.format(dq,4*r+3),'d{:d}lz'.format(dq),'d{:d}p'.format(dq-1))
                TN=TN&qtn.Tensor(DatTermTensor(pi3),indsD,tagsD)
        else:
            for dq in range(0,(nq-1)//2):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            tagsD=('D{:d}T{:d}'.format((nq-1)//2,4*r+3))
            indsD=('d{:d}g{:d}'.format((nq-1)//2,4*r+2),'d{:d}g{:d}'.format((nq-1)//2,4*r+3))
            TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
            for dq in range(((nq-1)//2)+1,nq):
                tagsD=('D{:d}T{:d}'.format(dq,4*r+3))
                indsD=('d{:d}g{:d}'.format(dq,4*r+2),'d{:d}g{:d}'.format(dq,4*r+3))
                TN=TN&qtn.Tensor(IdTensor(pi3),indsD,tagsD)
        #now do the ancilla measurement failure locations
        for aq in range(0,nq-1):
            tagsA=('A{:d}T{:d}'.format(aq,4*r+3))
            indsA=('a{:d}g{:d}'.format(aq,4*r+2),)
            TN=TN&qtn.Tensor(AMeasTensor1(pm),indsA,tagsA)
    #finally, put in the logical control tensor
    tagsL=('LZ')
    indsL=('lz',)
    for dq in range(0,nq):
        indsL=indsL+('d{:d}lz'.format(dq),)
    TN=TN&qtn.Tensor(ghzd,indsL,tagsL)
    return TN


#The following function takes a TN pure error pattern and translates it to the corresponding syndrome, or vice-versa

def syntranslator(nq,syn):
    nsb=nq-1
    outsyn=np.zeros(nsb)
    reps=np.zeros((nsb,nsb))
    for i in range(0,nsb//2):
        reps[i][i]=1
        if i>0:
            reps[i][i-1]=1
    for i in range(nsb//2,nsb):
        reps[i][i]=1
        if i<nsb-1:
            reps[i][i+1]=1
    for i in range(0,nsb):
        if syn[i]==1:
            outsyn=np.mod(outsyn+reps[i],2)
    return outsyn

#The following function creates a table of translations between TN pure error patterns and their corresponding syndromes

def syntranslatortable(nq):
    syntrantable=np.zeros(2**(nq-1))
    nsb=nq-1
    for sa in range(0,2**(nq-1)):
        syn=bitmapta[sa]
        outsyn=np.zeros(nsb)
        reps=np.zeros((nsb,nsb))
        for i in range(0,nsb//2):
            reps[i][i]=1
            if i>0:
                reps[i][i-1]=1
        for i in range(nsb//2,nsb):
            reps[i][i]=1
            if i<nsb-1:
                reps[i][i+1]=1
        for i in range(0,nsb):
            if syn[i]==1:
                outsyn=np.mod(outsyn+reps[i],2)
        syntrantable[sa]=int(get_int(outsyn))
    return syntrantable

#the following two functions are helper functions that aid Pfail and PfailQC in doing their calculations
#in a way that avoids memory issues for large code sizes

def TNReducer(TN,sobits,nq,nr,shave):
    counter=0
    for rep in range(0,nr):
        for b in range(0,nq-1):
            if counter<shave:
                bit=sobits[rep*(nq-1)+b]
                indb=('sf{:d}r{:d}'.format(b,rep),)
                tagb=('SF{:d}r{:d}BIT{:d}'.format(b,rep,bit))
                if bit==0:
                    TN=TN&qtn.Tensor(np.array([1,0]),indb,tagb)
                elif bit==1:
                    TN=TN&qtn.Tensor(np.array([0,1]),indb,tagb)
                else:
                    print('BAD')
            counter=counter+1
    TNC=TN.contract(all,optimize='auto-hq')
    out=TNC.data
    return out

def PfailSlice(TNSlice,nq,nr,shave,pfail):
    TNdat=np.asarray(TNSlice.data)
    TNdat=TNdat.reshape((2**(nr*(nq-1)+nq-shave-1),2))
    for soa in range(0,2**(nr*(nq-1)+nq-shave-1)):
        idp=TNdat[soa][0]
        lzp=TNdat[soa][1]
        if idp<lzp:
            pfail=pfail+idp
        else:
            pfail=pfail+lzp
    return pfail

#The following function calculates the exact logical failure rate for the repetition code on nq qubits with
#nr syndrome measurement repetitions using exact maximum likelihood decoding (via tensor network). Noise strength
# is p=k1/k2

def Pfail(nq,nr,p):
    TN=CLTNRepHWFull(nq,nr,p)
    tbits=nr*(nq-1)+nq
    pfail=0
    if tbits>24:
        nshave=tbits-14-nq
        shavebits=get_binall(np.arange(2**(nshave)).byteswap())[:,(64-nshave):]
        for shave in range(0,int(2**nshave)):
            bits=shavebits[shave]
            TNR=TNReducer(TN,bits,nq,nr,nshave)
            pfail=PfailSlice(TNR,nq,nr,nshave,pfail)
    else:
        TNC=TN.contract(all,optimize='auto-hq')
        TNdat=TNC.data
        TNdat=TNdat.reshape((2**(nr*(nq-1)),2**(nq-1),2))
        for so in range(0,2**(nr*(nq-1))):
            for sa in range(0,2**(nq-1)):
                idp=TNdat[so][sa][0]
                lzp=TNdat[so][sa][1]
                if idp<lzp:
                    pfail=pfail+idp
                else:
                    pfail=pfail+lzp
    return pfail

#The following function does the same thing as Pfail, but it also allows the final "perfect round" of measurement
#on the data qubits to experience measurement noise.

def PfailQC(nq,nr,p):
    TN=CLTNRepHWQC(nq,nr,p)
    tbits=nr*(nq-1)+nq
    pfail=0
    if tbits>24:
        nshave=tbits-14-nq
        shavebits=get_binall(np.arange(2**(nshave)).byteswap())[:,(64-nshave):]
        for shave in range(0,2**nshave):
            bits=shavebits[shave]
            TNR=TNReducer(TN,bits,nq,nr,nshave)
            pfail=PfailSlice(TNR,nq,nr,nshave,pfail)
    else:
        TNC=TN.contract(all,optimize='auto-hq')
        TNdat=TNC.data
        TNdat=TNdat.reshape((2**(nr*(nq-1)),2**(nq-1),2))
        for so in range(0,2**(nr*(nq-1))):
            for sa in range(0,2**(nq-1)):
                idp=TNdat[so][sa][0]
                lzp=TNdat[so][sa][1]
                if idp<lzp:
                    pfail=pfail+idp
                else:
                    pfail=pfail+lzp
    return pfail

