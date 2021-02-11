# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:35:21 2020

@author: Matteo Leonardi
"""

import gwinc
#import opts
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import collections

from gwinc       import plot
from gwinc_kagra import gwinc_kagra
from gwinc_kagra import precomp
from gwinc_kagra import kagra_sensitivity

ifom = 'ifo/bKAGRA.yaml'

freq = np.logspace(0, 4, 1000)
ifo  = gwinc.load_ifo(ifom)
ifo.Optics.Type = 'SignalRecycled'
#print(ifo.Optics.PhotoDetectorEfficiency)

#ifo.Squeezer.FilterCavity.Lrt = 10e-6
#precomp.temp_to_power(ifo)
#print(ifo)


#f = open('300m.dat', 'r')
#prova = f.read()


cvf = 'KAGRAcomp' # FIStest, FISopt, AFCtest, AFCopt, AFCSQZopt, FCopt, KAGRAcomp, KAGRASRM, KAGRAFDScomp

#conf = 'DRSE'
#noises = gwinc_kagra(freq, ifo, conf, True)
#ax.loglog(freq, np.sqrt(noises.get('Quantum Vacuum')), color='red',  label=conf, linestyle='--')
#%% comparison between high QEPD and current one (QE = 0.9)
if cvf == 'FIStest':
    fig = plt.figure(dpi=100, facecolor='white')
    ax  = fig.add_subplot(2, 1, 1)

    conf = 'BRSE'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises_BRSE = gwinc_kagra(freq, ifo, conf, True)
    rng=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4)
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Quantum Vacuum')), color='red',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Total')), color='green',  label='Total '+conf, linestyle='--')

    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises = gwinc_kagra(freq, ifo, conf, True,9)
    rng2=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
    
    ax.loglog(freq, np.sqrt(noises.get('Quantum Vacuum')), color='blue',  label='Quantum w/ high QE PD', linestyle='-')
    ax.loglog(freq, np.sqrt(noises.get('Total')), color='black',  label='Total w/ high QE PD', linestyle='--')
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_ylim([1e-24, 5e-21])
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    #ax.set_xticklabels([])
    ax.legend(ncol=1, fontsize='small',loc='upper right')
    
    
    fig.text(0.51, 0.99, 'BRSE with FIS', ha='center', va='top', size=15)
    fig.show()
    
    print('Original BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4))
    print('Improved BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4))
    print('Original BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 30))
    print('Improved BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30))
    
    ax  = fig.add_subplot(2, 1, 2)
    ax.semilogx(freq, -20*np.log10(np.sqrt(noises_BRSE.get('Quantum Vacuum'))/np.sqrt(noises.get('Quantum Vacuum'))), color='red',  label='ratio', linestyle='-')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("Relative difference [dB]")
    
    
    # png = 'figs/highQEBRSEcomparison.png'
    # plt.savefig(png, dpi = 200)
    #print('Figure is saved as '+png)
    
#%% FIS optimization
if cvf == 'FISopt':
    SQZlevel = np.linspace(0,15,16)
    rangeBNS = np.full(len(SQZlevel),0.)
    rangeBBHl = np.full(len(SQZlevel),0.)
    rangeBBHh = np.full(len(SQZlevel),0.)
    
    for i in range(0,16):
        # print(SQZlevel[i])
        sqzinj = SQZlevel[i]
        conf = 'BRSEFIS'
        ifo.Optics.PhotoDetectorEfficiency = 0.99
        noises = gwinc_kagra(freq, ifo, conf, False,sqzinj)
        rngBNS=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
        rngBBHl=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
        rngBBHh=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
        if i==0:
            rngBNSorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
            rngBBHlorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
            rngBBHhorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
        
        rangeBNS[i]=rngBNS/rngBNSorig
        rangeBBHl[i]=rngBBHl/rngBBHlorig
        rangeBBHh[i]=rngBBHh/rngBBHhorig
        
    fig = plt.figure(dpi=100, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    ax.plot(SQZlevel, rangeBNS, label='BNS')
    ax.plot(SQZlevel, rangeBBHl, label='light BBH')
    ax.plot(SQZlevel, rangeBBHh, label='heavy BBH')
    
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([0, 15])
    ax.set_ylim([0.3, 1.1])
    ax.set_xlabel("Injected FIS [dB]")
    ax.set_ylabel("Range improvement")
    #ax.set_xticklabels([])
    ax.legend(ncol=1, fontsize='medium',loc='center right')
    fig.text(0.51, 0.99, 'Range improvement of BRSE with FIS', ha='center', va='top', size=15)
    fig.text(0.6, 0.3, 'BNS range w/o FIS:   %.1f MPc' %rngBNSorig, ha='right', va='top', size=10)
    fig.text(0.6, 0.25, 'light BBH range w/o FIS:   %.1f MPc' %rngBBHlorig, ha='right', va='top', size=10)
    fig.text(0.6, 0.2, 'heavy BBH range w/o FIS:   %.1f MPc' %rngBBHhorig, ha='right', va='top', size=10)
    # fig.show()
    png = 'figs/rangeBRSE_FIS.png'
    plt.savefig(png, dpi = 200)
    print('Figure is saved as '+png)
  
#%% AFC test
if cvf =='AFCtest':
    fig = plt.figure(dpi=100, facecolor='white')
    ax  = fig.add_subplot(2, 1, 1)

    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises_BRSE = gwinc_kagra(freq, ifo, conf, False, 5)
    rng=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4)
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Quantum Vacuum')), color='red',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Total')), color='green',  label='Total '+conf, linestyle='--')

    conf = 'BRSEAFC'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises = gwinc_kagra(freq, ifo, conf, False, 5, 10000e-6)
    rng2=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
    ax.loglog(freq, np.sqrt(noises.get('Quantum Vacuum')), color='blue',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises.get('Total')), color='black',  label='Total '+conf, linestyle='--')
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_ylim([1e-24, 5e-21])
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    #ax.set_xticklabels([])
    ax.legend(ncol=1, fontsize='small',loc='upper right')
    
    
    fig.text(0.51, 0.99, 'BRSE with FIS', ha='center', va='top', size=15)
    # fig.show()
    
    print('Original BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4))
    print('Improved BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4))
    print('Original BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 30))
    print('Improved BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30))
    
    ax  = fig.add_subplot(2, 1, 2)
    ax.semilogx(freq, -20*np.log10(np.sqrt(noises_BRSE.get('Quantum Vacuum'))/np.sqrt(noises.get('Quantum Vacuum'))), color='red',  label='ratio', linestyle='-')
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("Relative difference [dB]")
    
    
    # png = 'figs/highQEBRSEcomparison.png'
    # plt.savefig(png, dpi = 200)
    #print('Figure is saved as '+png)

#%% AFC optimization
if cvf == 'AFCopt':
    losses = np.linspace(0,1200e-6,121)
    rangeBNS = np.full(len(losses),0.)
    rangeBBHl = np.full(len(losses),0.)
    rangeBBHh = np.full(len(losses),0.)
    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises = gwinc_kagra(freq, ifo, conf, False,0)
    rngBNSns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
    rngBBHlns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
    rngBBHhns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)

    
    for i in range(0,121):
        # print(SQZlevel[i])
        lossparam = losses[i]
        conf = 'BRSEAFC'
        ifo.Optics.PhotoDetectorEfficiency = 0.99
        noises = gwinc_kagra(freq, ifo, conf, False,5,lossparam)
        rngBNS=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
        rngBBHl=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
        rngBBHh=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
        if i==0:
            # conf = 'BRSEFIS'
            # ifo.Optics.PhotoDetectorEfficiency = 0.99
            # noises = gwinc_kagra(freq, ifo, conf, False,0)
            rngBNSorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
            rngBBHlorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
            rngBBHhorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
        
        rangeBNS[i]=rngBNS/rngBNSorig
        rangeBBHl[i]=rngBBHl/rngBBHlorig
        rangeBBHh[i]=rngBBHh/rngBBHhorig
        
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(2, 1, 1)
    ax.plot(losses/60*1e6, rangeBNS, label='BNS')
    ax.plot(losses/60*1e6, rangeBBHl, label='light BBH')
    ax.plot(losses/60*1e6, rangeBBHh, label='heavy BBH')

    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([0, 20])
    # ax.set_ylim([0.95, 1.35])
    ax.set_ylabel("Range improvement\n w.r.t. BRSE w/ FIS")
    # ax.set_xlabel("Input transmission = Loss [ppm/m]")
    ax.legend(ncol=1, fontsize='medium',loc='center right')
    fig.text(0.5, 0.97, 'Range improvement of BRSE with AFC (5dB inj)', ha='center', va='top', size=15)

    ax  = fig.add_subplot(2, 1, 2)
    ax.plot(losses/60*1e6, rangeBNS/rngBNSns*rngBNSorig, label='BNS')
    ax.plot(losses/60*1e6, rangeBBHl/rngBBHlns*rngBBHlorig, label='light BBH')
    ax.plot(losses/60*1e6, rangeBBHh/rngBBHhns*rngBBHhorig, label='heavy BBH')
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([0, 20])
    # ax.set_ylim([0.95, 1.35])
    ax.set_ylabel("Range improvement\n w.r.t BRSE")
    ax.set_xlabel("Input transmission = Loss [ppm/m]")
    ax.legend(ncol=1, fontsize='medium',loc='best')
    
    png = 'figs/rangeBRSE_AFC.png'
    plt.savefig(png, dpi = 200)
    print('Figure is saved as '+png)
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(2, 1, 1)
    AFCbw = 3e8/120*losses/(3.14*np.sqrt(1-losses))
    ax.plot(AFCbw, rangeBNS, label='BNS')
    ax.plot(AFCbw, rangeBBHl, label='light BBH')
    ax.plot(AFCbw, rangeBBHh, label='heavy BBH')

    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([0, 950])
    # ax.set_ylim([0.95, 1.35])
    ax.set_ylabel("Range improvement\n w.r.t. BRSE w/ FIS")
    # ax.set_xlabel("Input transmission = Loss [ppm/m]")
    ax.legend(ncol=1, fontsize='medium',loc='center right')
    fig.text(0.5, 0.97, 'Range improvement of BRSE with AFC (5dB inj)', ha='center', va='top', size=15)

    ax  = fig.add_subplot(2, 1, 2)
    ax.plot(AFCbw, rangeBNS/rngBNSns*rngBNSorig, label='BNS')
    ax.plot(AFCbw, rangeBBHl/rngBBHlns*rngBBHlorig, label='light BBH')
    ax.plot(AFCbw, rangeBBHh/rngBBHhns*rngBBHhorig, label='heavy BBH')
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([0, 950])
    # ax.set_ylim([0.95, 1.35])
    ax.set_ylabel("Range improvement\n w.r.t BRSE")
    ax.set_xlabel("AFC bandwidth [Hz]")
    ax.legend(ncol=1, fontsize='medium',loc='best')
    
    png = 'figs/rangeBRSE_AFCbw.png'
    plt.savefig(png, dpi = 200)
    print('Figure is saved as '+png)
    
#%% AFC+SQZlvl optimization
if cvf == 'AFCSQZopt':
    cyclei = 61
    losses = np.linspace(0,1200e-6,cyclei)
    cyclej = 61
    SQZlevel = np.linspace(0,15,cyclej)
    rangeBNS = np.full((len(losses),len(SQZlevel)),0.)
    rangeBBHl = np.full((len(losses),len(SQZlevel)),0.)
    rangeBBHh = np.full((len(losses),len(SQZlevel)),0.)
    # losses = np.linspace(0,1200e-6,121)
    # rangeBNS = np.full(len(losses),0.)
    # rangeBBHl = np.full(len(losses),0.)
    # rangeBBHh = np.full(len(losses),0.)
    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises = gwinc_kagra(freq, ifo, conf, False,0)
    rngBNSns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
    rngBBHlns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
    rngBBHhns=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
    

    
    for i in range(0,cyclei):
        # print(SQZlevel[i])
        lossparam = losses[i]
        for j in range(0,cyclej):
            sqzinj = SQZlevel[j]
            conf = 'BRSEAFC'
            ifo.Optics.PhotoDetectorEfficiency = 0.99
            noises = gwinc_kagra(freq, ifo, conf, False,sqzinj,lossparam)
            rngBNS=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
            rngBBHl=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
            rngBBHh=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
            # if i==0:
            #     # conf = 'BRSEFIS'
            #     # ifo.Optics.PhotoDetectorEfficiency = 0.99
            #     # noises = gwinc_kagra(freq, ifo, conf, False,0)
            #     rngBNSorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
            #     rngBBHlorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
            #     rngBBHhorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
            
            rangeBNS[i,j]=rngBNS
            rangeBBHl[i,j]=rngBBHl
            rangeBBHh[i,j]=rngBBHh
    
    
    
    print('Max BNS range:  %(v1).1f MPc, with loss = %(v2).1f ppm/m and %(v3).1f dB of injected FIS' %{'v1': np.max(rangeBNS), 'v2': losses[np.where(rangeBNS==np.max(rangeBNS))[0]]/60*1e6, 'v3': SQZlevel[np.where(rangeBNS==np.max(rangeBNS))[1]]})
    # print('Max light BBH range:  %(v1).1f MPc, with loss = %(v2).1f ppm/m and %(v3).1f dB of injected FIS' %{'v1': np.max(rangeBBHl), 'v2': losses[np.where(rangeBBHl==np.max(rangeBBHl))[0]]/60*1e6, 'v3': SQZlevel[np.where(rangeBBHl==np.max(rangeBBHl))[1]]})
    # print('Max heavy BBH max range:  %.1f MPc' %np.max(rangeBBHh))
    
    # print(losses[np.where(rangeBBHl==np.max(rangeBBHl))[0]]/60*1e6)
    # print(SQZlevel[np.where(rangeBBHl==np.max(rangeBBHl))[1]])
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses/60*1e6, rangeBNS/rngBNSns)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("Input transmission = Loss [ppm/m]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'BNS improvement with AFC', ha='center', va='top', size=15)
    png = 'figs/BNSrangeBRSE_AFCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses/60*1e6, rangeBBHl/rngBBHlns)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("Input transmission = Loss [ppm/m]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'Light BBH improvement with AFC', ha='center', va='top', size=15)
    png = 'figs/BBHlrangeBRSE_AFCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses/60*1e6, rangeBBHh/rngBBHhns)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("Input transmission = Loss [ppm/m]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'Heavy BBH improvement with AFC', ha='center', va='top', size=15)
    png = 'figs/BBHhrangeBRSE_AFCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
        

#%% FC optimization
if cvf == 'FCopt':
    cyclei = 21
    losses = np.linspace(0,50e-6,cyclei)
    cyclej = 31
    SQZlevel = np.linspace(0,15,cyclej)
    rangeBNS = np.full((len(losses),len(SQZlevel)),0.)
    rangeBBHl = np.full((len(losses),len(SQZlevel)),0.)
    rangeBBHh = np.full((len(losses),len(SQZlevel)),0.)
    
    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises = gwinc_kagra(freq, ifo, conf, False,0)
    rngBNSorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
    rngBBHlorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
    rngBBHhorig=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
    
    for i in range(0,cyclei):
        # print(SQZlevel[i])
        lossparam = losses[i]
        for j in range(0,cyclej):
            sqzinj = SQZlevel[j]
            conf = 'BRSEFC'
            ifo.Optics.PhotoDetectorEfficiency = 0.99
            noises = gwinc_kagra(freq, ifo, conf, False,sqzinj,lossparam)
            rngBNS=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 1.4)
            rngBBHl=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 30)
            rngBBHh=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises.get('Total')), 75)
        
            rangeBNS[i,j]=rngBNS
            rangeBBHl[i,j]=rngBBHl
            rangeBBHh[i,j]=rngBBHh
        
    
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses*1e6, rangeBNS/rngBNSorig)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("RTL [ppm]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'BNS improvement with FC (L=60m)', ha='center', va='top', size=15)
    png = 'figs/BNSrangeBRSE_FCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses*1e6, rangeBBHl/rngBBHlorig)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("RTL [ppm]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'Light BBH improvement with FC (L=60m)', ha='center', va='top', size=15)
    png = 'figs/BBHlrangeBRSE_FCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)
    im = ax.contourf(SQZlevel, losses*1e6, rangeBBHh/rngBBHhorig)
    cbar = fig.colorbar(im)
    ax.set_xlabel("FIS injected [dB]")
    ax.set_ylabel("RTL [ppm]")
    cbar.set_label('Range improvement w.r.t. BRSE')
    fig.text(0.5, 0.97, 'Heavy BBH improvement with FC (L=60m)', ha='center', va='top', size=15)
    png = 'figs/BBHhrangeBRSE_FCSQZ.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
    # ax.xaxis.grid(True)
    # ax.yaxis.grid(True)
    # # ax.set_xlim([0, 20])
    # # ax.set_ylim([0.95, 1.35])
    # ax.set_ylabel("Range improvement")
    # ax.set_xlabel("RTL [ppm/m]")
    # #ax.set_xticklabels([])
    # ax.legend(ncol=1, fontsize='medium',loc='best')
    # # fig.text(0.51, 0.99, 'Range improvement of BRSE with AFC', ha='center', va='top', size=15)
    # # fig.text(0.51, 0.99, 'Range improvement of BRSE with FIS', ha='center', va='top', size=15)
    # # fig.show()
    # png = 'figs/rangeBRSE_FClosses.png'
    # plt.savefig(png, dpi = 200)
    # print('Figure is saved as '+png)

#%% sensitivity comparison
if cvf =='KAGRAcomp':
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(2, 1, 1)

    conf = 'BRSE'
    ifo.Optics.PhotoDetectorEfficiency = 0.9
    noises_BRSE = gwinc_kagra(freq, ifo, conf, False)
    rng_BRSE=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4)
    # ax.loglog(freq, np.sqrt(noises_BRSE.get('Quantum Vacuum')), color='red',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Total')), color='black',  label='Total '+conf, linestyle='--')
    f = open('BRSE.dat', 'w')
    for i in np.arange(1000):
        f.write('{:.4E} {:.4E} {:.4E}\n'.format(noises_BRSE.get('Freq')[i], (noises_BRSE.get('Quantum Vacuum'))[i], (noises_BRSE.get('Total'))[i]))
    f.close()

    conf = 'BRSEFIS'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises_BRSEFIS = gwinc_kagra(freq, ifo, conf, False, 5)
    rng_BRSEFIS=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFIS.get('Total')), 1.4)
    # # ax.loglog(freq, np.sqrt(noises_BRSEFIS.get('Quantum Vacuum')), color='blue',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSEFIS.get('Total')), color='red',  label='Total '+conf, linestyle='-')
    f = open('BRSEFIS.dat', 'w')
    for i in np.arange(1000):
        f.write('{:.4E} {:.4E} {:.4E}\n'.format(noises_BRSEFIS.get('Freq')[i], (noises_BRSEFIS.get('Quantum Vacuum'))[i], (noises_BRSEFIS.get('Total'))[i]))
    f.close()
    
    
    conf = 'BRSEAFC'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    firsttime = False
    if firsttime:
        noises_BRSEAFC = gwinc_kagra(freq, ifo, conf, False, 5, 180e-6)
        f = open('BRSEAFC.dat', 'w')
        for i in np.arange(1000):
            f.write('{:.4E} {:.4E} {:.4E}\n'.format(noises_BRSEAFC.get('Freq')[i], (noises_BRSEAFC.get('Quantum Vacuum'))[i], (noises_BRSEAFC.get('Total'))[i]))
        f.close()
    else:
        f=open("BRSEAFC.dat", "r")
        if f.mode == 'r':
            # df = pd.read_table('BRSEAFC.dat', delim_whitespace=True, header=None)
            data = np.loadtxt('BRSEAFC.dat', dtype=float)
            noises_BRSEAFC=collections.OrderedDict(Freq=data[:,0], QuantumVacuum=data[:,1], Total=data[:,2])
            # print(noises_BRSEAFC)
        
    rng_BRSEAFC=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEAFC.get('Total')), 1.4)
    # ax.loglog(freq, np.sqrt(noises_BRSEAFC.get('Quantum Vacuum')), color='blue',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSEAFC.get('Total')), color='blue',  label='Total '+conf, linestyle='-')
    
    conf = 'BRSEFC'
    ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises_BRSEFC = gwinc_kagra(freq, ifo, conf, False, 13.5, 15e-6)
    rng_BRSEFC=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFC.get('Total')), 1.4)
    # ax.loglog(freq, np.sqrt(noises_BRSEFC.get('Quantum Vacuum')), color='blue',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSEFC.get('Total')), color='green',  label='Total '+conf, linestyle='-')
    f = open('BRSEFC.dat', 'w')
    for i in np.arange(1000):
        f.write('{:.4E} {:.4E} {:.4E}\n'.format(noises_BRSEFC.get('Freq')[i], (noises_BRSEFC.get('Quantum Vacuum'))[i], (noises_BRSEFC.get('Total'))[i]))
    f.close()
    
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_ylim([1e-24, 5e-21])
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    #ax.set_xticklabels([])
    ax.legend(ncol=1, fontsize='small',loc='upper right')
    
    
    fig.text(0.51, 0.99, 'BRSE, FIS, AFC and FC comparison', ha='center', va='top', size=15)
    fig.show()
    
    print('BRSE BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4))
    print('BRSEFIS BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFIS.get('Total')), 1.4))
    print('BRSEAFC BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEAFC.get('Total')), 1.4))
    print('BRSEFC BNS range:   %.1f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFC.get('Total')), 1.4))
    print('BRSE light BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 30))
    print('BRSEFIS light BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFIS.get('Total')), 30))
    print('BRSEAFC light BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEAFC.get('Total')), 30))
    print('BRSEFC light BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFC.get('Total')), 30))
    print('BRSE heavy BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 75))
    print('BRSEFIS heavy BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFIS.get('Total')), 75))
    print('BRSEAFC heavy BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEAFC.get('Total')), 75))
    print('BRSEFC heavy BBH range:   %.0f MPc' % kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEFC.get('Total')), 75))

    
    ax  = fig.add_subplot(2, 1, 2)
    ax.semilogx(freq, -20*np.log10(np.sqrt(noises_BRSE.get('Quantum Vacuum'))/np.sqrt(noises_BRSEFIS.get('Quantum Vacuum'))), color='red',  label='BRSEFIS/BRSE QN ratio', linestyle='-')
    ax.semilogx(freq, -20*np.log10(np.sqrt(noises_BRSE.get('Quantum Vacuum'))/np.sqrt(noises_BRSEAFC.get('QuantumVacuum'))), color='blue',  label='BRSEAFC/BRSE QN ratio', linestyle='-')
    ax.semilogx(freq, -20*np.log10(np.sqrt(noises_BRSE.get('Quantum Vacuum'))/np.sqrt(noises_BRSEFC.get('Quantum Vacuum'))), color='green',  label='BRSEFC/BRSE QN ratio', linestyle='-')

    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("QN ratio \nw.r.t. BRSE [dB]") 
    ax.legend(ncol=1, fontsize='small',loc='upper right')
    png = 'figs/KAGRASQZcomp.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
#%% sensitivity comparison
if cvf =='KAGRASRM':
    fig = plt.figure(dpi=300, facecolor='white')
    ax  = fig.add_subplot(1, 1, 1)

    conf = 'BRSE'
    #ifo.Optics.PhotoDetectorEfficiency = 0.99
    noises_BRSE = gwinc_kagra(freq, ifo, conf, False)
    rng_BRSE=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSE.get('Total')), 1.4)
    # ax.loglog(freq, np.sqrt(noises_BRSE.get('Quantum Vacuum')), color='red',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSE.get('Total')), color='black',  label='Total '+conf, linestyle='-')    
    
    conf = 'BRSE'
    ifo.Optics.SRM.Transmittance = 0.3
    noises_BRSEnSRM = gwinc_kagra(freq, ifo, conf, False)
    rng_BRSEnSRM=kagra_sensitivity.inspiralrange(freq, np.sqrt(noises_BRSEnSRM.get('Total')), 1.4)
    # ax.loglog(freq, np.sqrt(noises_BRSE.get('Quantum Vacuum')), color='red',  label='Quantum '+conf, linestyle='-')
    ax.loglog(freq, np.sqrt(noises_BRSEnSRM.get('Total')), color='red',  label='Total '+conf+' 30% T_SRM', linestyle='-')    
    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim([5, 5000])
    ax.set_ylim([1e-24, 5e-21])
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    #ax.set_xticklabels([])
    ax.legend(ncol=1, fontsize='small',loc='upper right')
    png = 'figs/KAGRASRMcomp.png'
    plt.savefig(png, dpi = 300)
    print('Figure is saved as '+png)
    
