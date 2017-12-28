from __future__ import print_function, division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font', family='serif')
from operator import itemgetter
import os
import math

from astropy.table import Table
from astropy.io import fits, ascii
from astropy import units as u

AP_PATH = os.environ['APOGEE_DATA']

import apogee_tools as ap

def plotRegion(**kwargs):

    specs  = kwargs.get('specs')
    params = kwargs.get('params')
    unc    = kwargs.get('unc')

    spc1, spc2 = specs[0], specs[1]
    par1, par2 = params[0], params[1]
    unc1, unc2 = unc[0], unc[1]

    wave1, flux1, name1 = spc1.wave, spc1.flux, spc1.name
    wave2, flux2, name2 = spc2.wave, spc2.flux, spc2.name

    save = kwargs.get('save', False)
    out  = kwargs.get('out', 'comparison.pdf')
    xlim = kwargs.get('xlim', [wave1[0], wave1[-1]])

    # x_range = xlim[1] - xlim[0]
    # scale = math.ceil(math.log10(x_range))

    fig = plt.figure(figsize=(10,4)) 
    ax  = fig.add_subplot(1,1,1)

    label1 = name1 + r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par1[0], unc1[0], par1[1], unc1[1], par1[2], unc1[2])
    label2 = name2 + r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par2[0], unc2[0], par2[1], unc2[1], par2[2], unc2[2])                                                               

    plt.plot(wave1, flux1, label=label1, alpha=.8, color='b')
    plt.plot(wave2, flux2, label=label2, alpha=.8, color='g')

    # major_ticks = np.arange(15000, 17000, 200)
    # minor_ticks = np.arange(15000, 17000, 50)
    # ax.set_xticks(major_ticks)                                                       
    # ax.set_xticks(minor_ticks, minor=True)   

    plt.xlim(xlim)  
    
    plt.legend(loc='upper right', fontsize=10)    
    plt.ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]')   
    plt.xlabel(r'$\lambda$ [$\mathring{A}$]')
    plt.tight_layout()

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()

    return fig


def plotMultiRegion(**kwargs):

    specs  = kwargs.get('specs')
    params = kwargs.get('params')
    unc    = kwargs.get('unc')
    nplots = kwargs.get('nplots', 2) # can only plot 2 right now
    regs   = kwargs.get('regs', [[15650,15780], [16150,16280]])
    save   = kwargs.get('save', False)
    out    = kwargs.get('out', 'Cannon_Regions.pdf')

    spc1, spc2 = specs[0], specs[1]
    par1, par2 = params[0], params[1]
    unc1, unc2 = unc[0], unc[1]

    wave1, flux1, name1 = spc1.wave, spc1.flux, spc1.name
    wave2, flux2, name2 = spc2.wave, spc2.flux, spc2.name

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16,4))

    label1 = name1 + r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par1[0], unc1[0], par1[1], unc1[1], par1[2], unc1[2])
    label2 = name2 + r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par2[0], unc2[0], par2[1], unc2[1], par2[2], unc2[2])   

    ax1.plot(wave1, flux1, label='M dwarf', alpha=.8, color='b')
    ax1.plot(wave2, flux2, label='M giant', alpha=.8, color='g')
    ax1.set_title('Cannon Region 1')
    ax1.set_xlim(regs[0])
    ax1.set_xlabel(r'$\lambda$ [$\mathring{A}$]')
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax1.legend(loc='upper left', fontsize=10) 

    ax2.plot(wave1, flux1, label=label1, alpha=.8, color='b')
    ax2.plot(wave2, flux2, label=label2, alpha=.8, color='g')
    ax2.set_title('Cannon Region 2')
    ax2.set_xlim(regs[1])
    ax2.set_xlabel(r'$\lambda$ [$\mathring{A}$]')  
    ax2.legend(loc='upper right', fontsize=10)  
         
    plt.tight_layout()

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()

    return fig


def plotMultiSpec(**kwargs):

    specs  = kwargs.get('specs')
    params = kwargs.get('params')
    unc    = kwargs.get('unc')
    save   = kwargs.get('save', False)
    out    = kwargs.get('out', 'multispec.pdf')
    regs   = kwargs.get('regs', [[15650,15780], [16150,16280]])
    save   = kwargs.get('save', False)
    out    = kwargs.get('out', 'Cannon_Regions.pdf')
    xdim, ydim = kwargs.get('dim', [1,2])

    if len(regs) != xdim:
        print('Warning, number of regions not equal to number of dimensions.')

    nplots = xdim*ydim
    axes   = [ax[i] for i in range(nplots)]

    for i in range(len(specs)):
        spc[i] = specs[i]
        par[i] = params[i]
        unc[i] = unc[i]

        wave[i], flux[i], name[i] = spc[i].wave, spc[i].flux, spc[i].name
    
    fig, (axes) = plt.subplots(3, sharex=True, sharey=True, figsize=(16, ydim*4))

    for i in range(nplots):
        ax[i] = fig.add_subplot(xdim,ydim,i)                                                                   
        ax[i].plot(wave1, flux1, alpha=.8, color='b')
        ax[i].plot(wave2, flux2, alpha=.8, color='g')
        ax[i].set_xlim(bands[j])
        ax[i].set_ylabel(r'$F_{\lambda}$')

    # # ax1 = fig.add_subplot(1,1,1)                                                                   
    # ax1.plot(wave1, flux1, alpha=.8, color='b')
    # ax1.plot(wave2, flux2, alpha=.8, color='g')
    # ax1.set_xlim(bands[0])
    # ax1.set_ylabel(r'$F_{\lambda}$')

    # # ax2 = fig.add_subplot(2,1,1)                                                               
    # ax2.plot(wave1, flux1, alpha=.8, color='b')
    # ax2.plot(wave2, flux2, alpha=.8, color='g')
    # ax2.set_xlim(bands[1])
    # ax2.set_ylabel(r'$F_{\lambda}$')

    # # ax3 = fig.add_subplot(3,1,1)                                                               
    # ax3.plot(wave1, flux1, alpha=.8, color='b')
    # ax3.plot(wave2, flux2, alpha=.8, color='g')
    # ax3.set_xlim(bands[2])
    # ax3.set_ylabel(r'$F_{\lambda}$')
    
    major_ticks = np.arange(15000, 17000, 200)
    minor_ticks = np.arange(15000, 17000, 50)
    ax3.set_xticks(major_ticks)                                                       
    ax3.set_xticks(minor_ticks, minor=True)   

    # plt.xlim(x_range)  
        
    # plt.ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]')   
    plt.xlabel(r'$\lambda$ [$\mathring{A}$]')

    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()


def plotFourSpecComparison(**kwargs):

    """
    Compare 4 spectra side by side, over Cannon regions
    """

    # Lists of spectra and their parameters
    specs  = kwargs.get('specs')
    params = kwargs.get('params')
    uncert = kwargs.get('unc')

    # optional
    save   = kwargs.get('save', True)
    out    = kwargs.get('out', 'Four_Spec_Comparison.pdf')
    regs   = kwargs.get('regs', [[15650,15780], [16150,16280]])

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(16,12))

    nspecs = len(specs)

    wave1, flux1 = specs[0].wave, specs[0].flux
    wave2, flux2 = specs[1].wave, specs[1].flux
    wave3, flux3 = specs[2].wave, specs[2].flux
    wave4, flux4 = specs[3].wave, specs[3].flux

    label1a = specs[0].name
    label2a = specs[1].name
    label3a = specs[2].name
    label4a = specs[3].name

    label1b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[0][0], params[0][1], params[0][2])
    label2b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[1][0], params[1][1], params[1][2])
    label3b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[2][0], params[2][1], params[2][2])
    label4b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[3][0], params[3][1], params[3][2])

    ax1.plot(wave1, flux1, label=label1a, alpha=.8, color='k')
    ax1.set_title('Cannon Region 1')
    ax1.set_xlim(regs[0])
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax1.legend(loc='upper left', fontsize=10) 

    ax2.plot(wave1, flux1, label=label1b, alpha=.8, color='k')
    ax2.set_title('Cannon Region 2')
    ax2.set_xlim(regs[1])
    ax2.legend(loc='upper right', fontsize=10) 

    ax3.plot(wave2, flux2, label=label2a, alpha=.8, color='k')
    ax3.set_xlim(regs[0])
    ax3.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax3.legend(loc='upper left', fontsize=10) 

    ax4.plot(wave2, flux2, label=label2b, alpha=.8, color='k')
    ax4.set_xlim(regs[1])
    ax4.legend(loc='upper right', fontsize=10) 

    ax5.plot(wave3, flux3, label=label3a, alpha=.8, color='k')
    ax5.set_xlim(regs[0])
    ax5.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax5.legend(loc='upper left', fontsize=10) 

    ax6.plot(wave3, flux3, label=label3b, alpha=.8, color='k')
    ax6.set_xlim(regs[1]) 
    ax6.legend(loc='upper right', fontsize=10) 

    ax7.plot(wave4, flux4, label=label4a, alpha=.8, color='k')
    ax7.set_xlim(regs[0])
    ax7.set_xlabel(r'$\lambda$ [$\mathring{A}$]')
    ax7.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax7.legend(loc='upper left', fontsize=10) 

    ax8.plot(wave4, flux4, label=label4b, alpha=.8, color='k')
    ax8.set_xlim(regs[1])
    ax8.set_xlabel(r'$\lambda$ [$\mathring{A}$]')  
    ax8.legend(loc='upper right', fontsize=10) 

    fig.subplots_adjust(hspace=0, wspace=.01)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()


def plotFourSpecComparisonWithModels(**kwargs):

    """
    Compare 4 spectra side by side with Cannon models, over Cannon regions
    """

    # Lists of spectra and their parameters
    specs  = kwargs.get('specs')
    params = kwargs.get('params')
    # uncert = kwargs.get('unc')

    # Cannon model info
    models = kwargs.get('models')
    labels = kwargs.get('labels')

    # optional
    save   = kwargs.get('save', True)
    out    = kwargs.get('out', 'Four_Spec_Comparison.pdf')
    regs   = kwargs.get('regs', [[15650,15780], [16150,16280]])

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(16,12))

    nspecs = len(specs)

    # Read in spectra wave and flux
    wave1, flux1 = specs[0].wave, specs[0].flux
    wave2, flux2 = specs[1].wave, specs[1].flux
    wave3, flux3 = specs[2].wave, specs[2].flux
    wave4, flux4 = specs[3].wave, specs[3].flux

    label1a = specs[0].name
    label2a = specs[1].name
    label3a = specs[2].name
    label4a = specs[3].name

    # Read in model wave and flux
    wl = models[0].wave

    synth_1 = models[0].flux
    synth_2 = models[1].flux
    synth_3 = models[2].flux
    synth_4 = models[3].flux

    # Create labels for the plots
    label1b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[0][0], params[0][1], params[0][2])
    label2b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[1][0], params[1][1], params[1][2])
    label3b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[2][0], params[2][1], params[2][2])
    label4b = r'$T_{eff}=%s, logg=%s, [Fe/H]=%s$' %(params[3][0], params[3][1], params[3][2])

    label1c = r'$T_{eff}=%s, [Fe/H]=%s$' %(labels[0][0], labels[0][1])
    label2c = r'$T_{eff}=%s, [Fe/H]=%s$' %(labels[1][0], labels[1][1])
    label3c = r'$T_{eff}=%s, [Fe/H]=%s$' %(labels[2][0], labels[2][1])
    label4c = r'$T_{eff}=%s, [Fe/H]=%s$' %(labels[3][0], labels[3][1])

    ax1.plot(wave1, flux1, label=label1a, alpha=.8, color='k')
    ax1.plot(wl, synth_1, label='Cannon Model', alpha=.8, color='r')
    ax1.set_title('Cannon Region 1')
    ax1.set_xlim(regs[0])
    ax1.set_ylim([0.6,1.2])
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax1.legend(loc='upper left', fontsize=10) 

    ax2.plot(wave1, flux1, label=label1b, alpha=.8, color='k')
    ax2.plot(wl, synth_1, label=label1c, alpha=.8, color='r')
    ax2.set_title('Cannon Region 2')
    ax2.set_xlim(regs[1])
    ax2.legend(loc='upper right', fontsize=10) 

    ax3.plot(wave2, flux2, label=label2a, alpha=.8, color='k')
    ax3.plot(wl, synth_2, label='Cannon Model', alpha=.8, color='r')
    ax3.set_xlim(regs[0])
    ax1.set_ylim([0.6,1.2])
    ax3.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax3.legend(loc='upper left', fontsize=10) 

    ax4.plot(wave2, flux2, label=label2b, alpha=.8, color='k')
    ax4.plot(wl, synth_2, label=label2c, alpha=.8, color='r')
    ax4.set_xlim(regs[1])
    ax4.legend(loc='upper right', fontsize=10) 

    ax5.plot(wave3, flux3, label=label3a, alpha=.8, color='k')
    ax5.plot(wl, synth_3, label='Cannon Model', alpha=.8, color='r')
    ax5.set_xlim(regs[0])
    ax1.set_ylim([0.6,1.2])
    ax5.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax5.legend(loc='upper left', fontsize=10) 

    ax6.plot(wave3, flux3, label=label3b, alpha=.8, color='k')
    ax6.plot(wl, synth_3, label=label3c, alpha=.8, color='r')
    ax6.set_xlim(regs[1]) 
    ax6.legend(loc='upper right', fontsize=10) 

    ax7.plot(wave4, flux4, label=label4a, alpha=.8, color='k')
    ax7.plot(wl, synth_4, label='Cannon Model', alpha=.8, color='r')
    ax7.set_xlim(regs[0])
    ax1.set_ylim([0.6,1.2])
    ax7.set_xlabel(r'$\lambda$ [$\mathring{A}$]')
    ax7.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax7.legend(loc='upper left', fontsize=10) 

    ax8.plot(wave4, flux4, label=label4b, alpha=.8, color='k')
    ax8.plot(wl, synth_4, label=label4c, alpha=.8, color='r')
    ax8.set_xlim(regs[1])
    ax8.set_xlabel(r'$\lambda$ [$\mathring{A}$]')  
    ax8.legend(loc='upper right', fontsize=10) 

    fig.subplots_adjust(hspace=0, wspace=.01)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()


def plotChiComparison(**kwargs):

    """
    This function plots a spectrum against two different models; labeled w/ chi-squared values
    """

    spec   = kwargs.get('spec')
    sp_lbl = kwargs.get('slbl')
    models = kwargs.get('models')
    params = kwargs.get('params')
    uncert = kwargs.get('unc')
    chival = kwargs.get('chis')
    regs   = kwargs.get('regs', [[15650,15780], [16150,16280]])
    save   = kwargs.get('save', False)
    out    = kwargs.get('out', 'Chi_Squared_Comparison.pdf')

    mdl1, mdl2 = models[0], models[1]
    par1, par2 = params[0], params[1]
    unc1, unc2 = uncert[0], uncert[1]

    spc_wave, spc_flux, spc_name = spec.wave, spec.flux, spec.name

    mdl_wave1, mdl_flux1, mdl_name1 = mdl1.wave, mdl1.flux, mdl1.name
    mdl_wave2, mdl_flux2, mdl_name2 = mdl2.wave, mdl2.flux, mdl2.name

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharey=True, figsize=(16,6))

    label1 = r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par1[0], unc1[0], par1[1], unc1[1], par1[2], unc1[2])
    label2 = r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(par2[0], unc2[0], par2[1], unc2[1], par2[2], unc2[2])   

    if sp_lbl != None:
        spc_par, spc_unc = sp_lbl[0], sp_lbl[1]
        label3 = r'$, T_{eff}=%s \pm %s, logg=%s \pm %s, [Fe/H]=%s \pm %s$' %(spc_par[0], spc_unc[0], spc_par[1], spc_unc[1], spc_par[2], spc_unc[2])
    else:
        label3 = ''

    ax1.plot(mdl_wave1, mdl_flux1, label='M dwarf', alpha=.8, color='b')
    ax1.plot(spc_wave, spc_flux, label=r'\chi^{2} = '+str(chival[0]), alpha=.65, color='k')
    ax1.set_title('Cannon Region 1')
    ax1.set_xlim(regs[0])
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax1.legend(loc='upper left', fontsize=10) 

    ax2.plot(mdl_wave1, mdl_flux1, label=mdl_name1+label1, alpha=.8, color='b')
    ax2.plot(spc_wave, spc_flux, label=spec.name+label3, alpha=.65, color='k')
    ax2.set_title('Cannon Region 2')
    ax2.set_xlim(regs[1])
    ax2.legend(loc='upper right', fontsize=10) 

    ax3.plot(mdl_wave2, mdl_flux2, label='M giant', alpha=.8, color='g')
    ax3.plot(spc_wave, spc_flux, label=r'\chi^{2} = '+str(chival[1]), alpha=.65, color='k')
    ax3.set_xlim(regs[0])
    ax3.set_xlabel(r'$\lambda$ [$\mathring{A}$]')
    ax3.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]') 
    ax3.legend(loc='upper left', fontsize=10) 

    ax4.plot(mdl_wave2, mdl_flux2, label=mdl_name2+label2, alpha=.8, color='g')
    ax4.plot(spc_wave, spc_flux, label=spec.name+label3, alpha=.65, color='k')
    ax4.set_xlim(regs[1])
    ax4.set_xlabel(r'$\lambda$ [$\mathring{A}$]')  
    ax4.legend(loc='upper right', fontsize=10)  
         
    plt.tight_layout()

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()

    return fig


def plotChiComparisonScatter(**kwargs):

    """
    Plot scatter of chi values for model1 vs. model2
    Can have three different lists which will be plotted in different colors:
    1) dwarf list = blue, 2) giant list = green, 3) unknown list = black
    """
    rc('font', family='serif')

    dlist = kwargs.get('dlist') 
    glist = kwargs.get('glist')
    ulist = kwargs.get('ulist')

    # Get chi values for model1 and model2 respectively
    if dlist != None:
        dvals1 = list(map(itemgetter(0), dlist))
        dvals2 = list(map(itemgetter(1), dlist))

    if glist != None:
        gvals1 = list(map(itemgetter(0), glist))
        gvals2 = list(map(itemgetter(1), glist))

    if ulist != None:
        uvals1 = list(map(itemgetter(0), ulist))
        uvals2 = list(map(itemgetter(1), ulist))

    save = kwargs.get('save', False)
    out  = kwargs.get('out', 'Chi_Squared_Scatter.pdf')

    fig, ax = plt.subplots()

    if dlist != None:
        ax.scatter(dvals1, dvals2, label='Known M dwarfs', color='b', s=1)

    if glist != None:
        ax.scatter(gvals1, gvals2, label='Known M giants', color='g', s=1)

    if ulist != None:
        ax.scatter(uvals1, uvals2, label='Unknown M star', color='k', alpha=.8, s=1)

    lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), \
             np.max([ax.get_xlim(), ax.get_ylim()]) ]

    # now plot both limits against eachother
    ax.plot(lims, lims, '--', alpha=0.75, color='r', zorder=0, markersize=10)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel('Dwarf Model')
    ax.set_ylabel('Giant Model')
    ax.legend(loc='upper right', fontsize=10) 

    plt.savefig('plots/'+str(out))

    plt.tight_layout()
    plt.show()
    plt.close()


def plotChiComparisonScatterByParam(**kwargs):

    """
    Plot scatter of chi values for model1 vs. model2
    Color indicates a parameter such as logg or teff
    """

    # required
    dchi = kwargs.get('dchi')
    gchi = kwargs.get('gchi')
    color_param = kwargs.get('cparam')

    # optional
    save = kwargs.get('save', True)
    out  = kwargs.get('out', 'Chi_Squared_Scatter_By_Param.pdf')

    fig, ax = plt.subplots()

    scat = ax.scatter(dchi, gchi, alpha=.8, s=2, c=np.clip(color_param,0,500), cmap='viridis')

    lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), \
             np.max([ax.get_xlim(), ax.get_ylim()]) ]

    # now plot both limits against eachother
    ax.plot(lims, lims, '--', alpha=0.75, color='r', zorder=0, markersize=10)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel(r'$\chi^{2}$ Dwarf Model')
    ax.set_ylabel(r'$\chi^{2}$ Giant Model')
    ax.set_title('ASPCAP SNR')

    # plt.clim(-1,5)
    fig.colorbar(scat)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.tight_layout()
    plt.show()
    plt.close()


def plotParameterScatter(**kwargs):

    """
    Plot commparison of 3 parameters: 2 axes and 1 color
    """

    params = kwargs.get('params')

    plt.scatter(params[0], params[1], c=params[2])

    plt.tight_layout()
    plt.show()
    plt.close()


def plotCannonSynth(**kwargs):

    for n in range(len(fluxes)):
    
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16,4))
        
        ax1.plot(wl, pseudo_tr_flux[n], label=train_)
        ax1.plot(wl, fluxes[n], label=labels[n])
        ax1.set_ylim([0.6,1.2])
        ax1.set_xlim([15650,15780])
        ax1.legend(loc='upper left', fontsize=10) 

        ax2.plot(wl, pseudo_tr_flux[n], label=tr_ID[n])
        ax2.plot(wl, fluxes[n], label=r'Cannon Model, \chi^{2} = ')
        ax2.set_xlim([16150,16280])
        ax2.legend(loc='upper right', fontsize=10) 

        plt.show()


def plotCrossValidation(trn_labels, crv_labels, **kwargs):

    # required
    label_names = kwargs.get('label_names', ['TEFF', '[Fe/H]'])

    # optional
    save = kwargs.get('save', False)
    out  = kwargs.get('out', 'Cross_Validation_Scatter.pdf')

    if len(label_names) == 1:
        fig, ax1 = plt.subplots(1, 1)

        trn_teffs = list(map(itemgetter(0), trn_labels))

        crv_teffs = list(map(itemgetter(0), crv_labels))

        ax1.scatter(trn_teffs, crv_teffs, color='k')

        limits_1 = [ np.min([ax1.get_xlim(), ax1.get_ylim()]), \
                 np.max([ax1.get_xlim(), ax1.get_ylim()]) ]

        # now plot both limits against eachother
        ax1.plot(limits_1, limits_1, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax1.set_aspect('equal')
        ax1.set_xlim(limits_1)
        ax1.set_ylim(limits_1)

        ax1.set_title(label_names[0] + 'Labels')
        ax1.set_xlabel('Training Label')
        ax1.set_ylabel('Cross-Validated Label')

    elif len(label_names) == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2)

        trn_teffs = list(map(itemgetter(0), trn_labels))
        trn_metal = list(map(itemgetter(1), trn_labels))

        crv_teffs = list(map(itemgetter(0), crv_labels))
        crv_metal = list(map(itemgetter(1), crv_labels))

        ax1.scatter(trn_teffs, crv_teffs, color='k')

        limits_1 = [ np.min([ax1.get_xlim(), ax1.get_ylim()]), \
                 np.max([ax1.get_xlim(), ax1.get_ylim()]) ]

        # now plot both limits against eachother
        ax1.plot(limits_1, limits_1, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax1.set_aspect('equal')
        ax1.set_xlim(limits_1)
        ax1.set_ylim(limits_1)

        ax1.set_title(label_names[0] + 'Labels')
        ax1.set_xlabel('Training Label')
        ax1.set_ylabel('Cross-Validated Label')

        ax2.scatter(trn_metal, crv_metal, color='k')

        limits_2 = [ np.min([ax2.get_xlim(), ax2.get_ylim()]), \
                 np.max([ax2.get_xlim(), ax2.get_ylim()]) ]

        ax2.plot(limits_2, limits_2, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax2.set_aspect('equal')
        ax2.set_xlim(limits_2)
        ax2.set_ylim(limits_2)

        ax2.set_title(label_names[1] + 'Labels')
        ax2.set_xlabel('Training Label')
        ax2.set_ylabel('Cross-Validated Label')

    plt.savefig('plots/'+str(out))

    plt.tight_layout()
    plt.show()
    plt.close()


def plotOHBands(**kwargs):

    """
    Plot OH-bands which are sensitive to Teff.
    Unless specified, 1) black = spectrum, 2) red = cannon model, 3) blue = btsettl, 4) green = other
    """

    # required
    specs = kwargs.get('specs')

    # optional 
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])

    bands = [[15400,15450], [16350,16360], [16860,16890]]
    nbands = len(bands)

    fig, axs = plt.subplots(1, nbands, figsize=(12,3))

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j])
        ax.set_xlim(bands[i])

    plt.tight_layout()
    plt.show()
    plt.close()


def plotSpectralSequence(**kwargs):

    """
    Plot spectral sequence, with flux offset. Optional: overplot models.
    """

    # required
    specs  = kwargs.get('specs')
    models = kwargs.get('models')
    labels = kwargs.get('labels')

    # optional
    save   = kwargs.get('save', True)
    out    = kwargs.get('out', 'Spectral_Sequence.pdf')

    nspecs = len(specs)

    plotbands = [[15200,15800],[15870,16420],[16500,16930]] 

    highlights1 = [[15290,15350],[15530,15590],[15620,15650],[15715,15780]] 
    highlights2 = [[15960,15980],[16150,16220]]
    highlights3 = [[16520,16560],[16695,16775]]
    
    fig1, ax1 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        # Replace zeroes with nans
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax1.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax1.plot(models[i].wave, models[i].flux + .5*i, color='r', alpha=.8)
        ax1.text(plotbands[0][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
        ax1.text(plotbands[0][0]+10, .5*i + 1.1, str(specs[i].name), fontweight='bold', fontsize=12)
    for h in highlights1:
        ax1.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax1.set_xlim(plotbands[0])
    ax1.set_ylim([.5,6])
    ax1.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)
    ax1.set_title('APOGEE M dwarf Spectral Sequence', fontsize=20)

    fig1.savefig('1_' + str(out))
    plt.show()
    plt.close(fig1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax2.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax2.plot(models[i].wave, mdl_fl + .5*i, color='r', alpha=.8)
        ax2.text(plotbands[1][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
    for h in highlights2:
        ax2.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax2.set_xlim(plotbands[1])
    ax2.set_ylim([.5,6])
    ax2.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax2.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)

    fig2.savefig('2_' + str(out))
    plt.show()
    plt.close(fig2)

    fig3, ax3 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax3.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax3.plot(models[i].wave, mdl_fl + .5*i, color='r', alpha=.8)
        ax3.text(plotbands[2][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
    for h in highlights3:
        ax3.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax3.set_xlim(plotbands[2])
    ax3.set_ylim([.5,6])
    ax3.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax3.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)

    fig3.savefig('3_' + str(out))
    plt.show()
    plt.close(fig3)


def plotBands(**kwargs):

    """
    Plot OH-bands which are sensitive to Teff.
    Unless specified, 1) black = spectrum, 2) red = cannon model, 3) blue = btsettl, 4) green = other
    """

    # required
    specs = kwargs.get('specs')
    plot_band = kwargs.get('bands', 'OH')

    # optional 
    title  = kwargs.get('title', None)
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])
    save   = kwargs.get('save', True)
    output = kwargs.get('out', 'Plot_Bands.pdf')

    if plot_band == 'OH':
        bands = [[15400,15450], [16350,16360], [16860,16890]]
    elif plot_band == 'Ca':
        bands = [[16131,16141], [16145,16155], [16152,16162]]
    elif plot_band == 'K':
        bands = [[15158,15168], [15163,15173]]
    elif plot_band == 'Mg':
        bands = [[15735,15745], [15743,15753], [15760, 15770]]
    elif plot_band == 'Al':
        bands = [[16713, 16723], [16745,16755], [16758,16768]]
    elif plot_band == 'Cannon':
        bands = [[15650,15780], [16150,16280]]
    elif plot_band == 'Full':
        bands = [[15200,15800],[15870,16420],[16490,16940]]
    #Regions from Rajpurohit paper:
    elif plot_band == 'R1':
        bands = [[15150,15450]]
    elif plot_band == 'R2':
        bands = [[15450,15800]]
    elif plot_band == 'R3':
        bands = [[15850,16420]]
    elif plot_band == 'R4':
        bands = [[16500,16910]]
    nbands = len(bands)

    fig, axs = plt.subplots(1, nbands, figsize=(16,4))

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j], alpha=.8)
        if i==0:
            ax.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]', fontsize=20)
        ax.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=20)
        ax.set_xlim(bands[i])
        ax.set_ylim([0.7, 1.15])
 
    if title != None:
        plt.suptitle(title, fontsize=25)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()


def plotRajBands(**kwargs):

    #required
    specs = kwargs.get('specs')
    title = kwargs.get('title', 'Spectrum_Raj_Bands')
    
    #optional
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])
    save   = kwargs.get('save', True)
    # output = kwargs.get('out', 'Plot_Bands.pdf')

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']
    bands  = [[15200,15450],[15450,15800],[15850,16420],[16500,16910]]

    fig, axs = plt.subplots(4, 1, figsize=(20,16), sharey=True) 

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j], alpha=.8, linewidth=1)
        ax.set_xlim(bands[i])
        if i == 0:
            ax.set_title(title, fontsize=20)

    if save == True:
        plt.savefig('plots/'+title+'.pdf')

    plt.show()
    plt.close()


# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# majorLocator = MultipleLocator(major)
#     majorFormatter = FormatStrFormatter('%d')
#     minorLocator = MultipleLocator(minor)

# ax.xaxis.set_major_locator(majorLocator)
#     ax.xaxis.set_major_formatter(majorFormatter)


