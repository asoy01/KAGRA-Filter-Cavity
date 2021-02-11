
import copy
import scipy

import numpy as np
from numpy import pi, sqrt

from gwinc.struct  import Struct
from gwinc.precomp import dhdl

def precomp_kagra(f, ifoin, verbose = 'None'):
    """Add precomputed data to the IFO model.

    This function DOES NOT modify IFO. Its return value is a copy of
    IFO with precomputations filled out.
    """
    ifo = copy.deepcopy(ifoin)

    if 'gwinc' not in ifo:
        ifo.gwinc = Struct()

    ##############################
    # optics values

    # calculate optics' parameters
    ifo.Materials.MirrorVolume = pi*ifo.Materials.MassRadius**2 * \
                                 ifo.Materials.MassThickness
    ifo.Materials.MirrorMass   = ifo.Materials.MirrorVolume * \
                                 ifo.Materials.Substrate.MassDensity
    ifo.Optics.ITM.Thickness   = ifo.Materials.MassThickness

    mass = ifo.Materials.MirrorMass

    #beam radius at ITM (m); scale size with respect to mass
    ifo.Optics.ITM.BeamRadius = 3.5e-2*(mass/22.8)**(1./3)

    #beam radius at ETM (m); scale size with respect to mass
    ifo.Optics.ETM.BeamRadius = 3.5e-2*(mass/22.8)**(1./3)

    ##############################
    # calc power and IFO parameters

    ## calc TM fiber suspension radius
    precompFiber(ifo)

    ## DETERMINATION OF LASER POWER AT BS ##########
    temp_to_power(ifo)

    ## calc finesse and arm power
    precompPower(ifo)

    ## calc optimal filter cavity parameters for FD squeezing
    precompFDSqz(ifo)


    ##############################
    # strain to length conversion
    #
    # used for converting displacement to strain in all noise calculations

    L = ifo.Infrastructure.Length
    ifo.gwinc.dhdl_sqr = 1/L**2
    ifo.gwinc.sinc_sqr = 1
    #ifo.gwinc.dhdl_sqr, ifo.gwinc.sinc_sqr = dhdl(f, L)

    if not verbose == 'None':
        print
        print(verbose)
        print('Input power:       %.0f W' % ifo.Laser.Power)
        print('Mirror mass:       %.2f kg' % ifo.Materials.MirrorMass)
        print('Power at BS:       %.0f W'  % ifo.gwinc.pbs)
        print('Fiber diameter:    %.2f mm' % \
            (ifo.Suspension.Stage[0].WireRadius*2000))

    return ifo

def params_kagra(ifoin, conf, SQZlvl,RTLlvl):
    """Update IFO parameters with the specified configuration

    This function DOES NOT modify IFO. Its return value is a copy of
    IFO with parameters filled out.
    """

    ifo = copy.deepcopy(ifoin)

    ifo.Materials.MirrorVolume = pi*ifo.Materials.MassRadius**2 * \
                                 ifo.Materials.MassThickness
    ifo.Materials.MirrorMass   = ifo.Materials.MirrorVolume * \
                                 ifo.Materials.Substrate.MassDensity

    phidet = pi/2-ifo.Optics.SRM.Tunephase/2    # SRC Detunning
    xi     =      ifo.Optics.Quadrature.dc      # Homodyne Readout phase
    tau    = sqrt(ifo.Optics.SRM.Transmittance) # SRM Transmittance [amplitude]
    rSRM   = sqrt(1 - tau**2)                   # SRM Reflectivity  [amplitude]
    Tm     = ifo.Materials.Substrate.Temp       # temperature at TM [k]
    I0attn = ifo.Laser.I0attenuation
    len3   = ifo.Suspension.Stage[0].Length     # wire2 [m]
    wsafe  = ifo.Suspension.Sapphire.wsafe
    mass   = ifo.Materials.MirrorMass
    LF     = ifo.Suspension.LF
    SQZ    = 0
    FC     = 0
#    if   conf == 'DRSE':        #detuned RSE
#    	phidet = pi/180*(90-3.5)
#        xi     = pi/180*135.1
    
    if conf == 'DRSE':
        phidet = pi/180*(90-3.5)
        xi     = pi/180*135.1  
        
    elif conf == 'BRSEhomo':    #broadband RSE with homodyne
        phidet = pi/180*90
        xi     = pi/180*119.1
    elif conf == 'BRSE':        #broadband RSE without homodyne
      	phidet = pi/180*90
      	xi     = pi/180*90
    elif conf == 'BRSEFIS':        #broadband RSE with FIS
        phidet = pi/180*90
        xi     = pi/180*90
        SQZ    = SQZlvl #4
        RTL    = RTLlvl 
        # FC     = 60 #30
        # rSRM   = sqrt(1 - 0.1536) 
#        len3   = 50.0e-2
#        wsafe  = 5
    elif conf == 'BRSEAFC':     #broadband RSE with FD squeezing
        phidet = pi/180*90
        xi     = pi/180*90
        SQZ    = SQZlvl #4
        FC     = 60 #30
        RTL    = RTLlvl 
        # rSRM   = sqrt(1 - 0.1536) 
#        len3   = 50.0e-2
#        wsafe  = 5
    elif conf == 'BRSEFC':     #broadband RSE with FD squeezing
        phidet = pi/180*90
        xi     = pi/180*90
        SQZ    = SQZlvl #4
        FC     = 60 #30
        RTL    = RTLlvl 
        # rSRM   = sqrt(1 - 0.1536) 
#        len3   = 50.0e-2
#        wsafe  = 5
    elif conf == 'BRSEsqz':     #broadband RSE with FD squeezing
        phidet = pi/180*90
        xi     = pi/180*90
        SQZ    = 14 #4
        FC     = 60 #30
        rSRM   = sqrt(1 - 0.1536) 
#        len3   = 50.0e-2
#        wsafe  = 5
    elif conf == 'LF':          #KAGRA+ LF
        phidet = pi/180*(90-28.5)
        xi     = pi/180*133.6
        LF     = 1
        Tm     = 23.6
        I0attn = 0.928
        rSRM   = sqrt(0.955)
        len3   = 99.8e-2
        wsafe  = 0.995
    elif conf == 'HF':          #KAGRA+ HF
        phidet = pi/180*(90-0.1)
        xi     = pi/180*97.1
        Tm     = 20.8
        rSRM   = sqrt(0.907)
        len3   = 20.1e-2
        wsafe  = 30.25
        SQZ    = 10
    elif conf == '40kg':        #KAGRA+ 40kg
        phidet = pi/180*(90-3.5)
        xi     = pi/180*123.2
        Tm     = 21.0
        rSRM   = sqrt(0.922)
        len3   = 28.6e-2
        mass   = 40
        wsafe  = 13.45
    elif conf == 'FDSQZ':       #KAGRA+ FDSQZ
        phidet = pi/180*(90-0.2)
        xi     = pi/180*93.1
        Tm     = 21.3
        rSRM   = sqrt(0.832)
        len3   = 23.0e-2
        wsafe  = 17.8
        SQZ    = 10
        FC     = 30
    else:
        print('Unknown configuration: ',conf,' DRSE is assumed')

    scale = (mass/ifo.Materials.MirrorMass)**(1./3)

    # Update IFO parameters
    ifo.Optics.SRM.Tunephase       = pi-phidet*2
    ifo.Optics.Quadrature.dc       = xi

    ifo.Optics.SRM.Transmittance   = 1-rSRM**2
    ifo.Materials.Substrate.Temp   = Tm
    ifo.Laser.I0attenuation        = I0attn
    ifo.Suspension.Stage[0].Length = len3
    ifo.Suspension.Sapphire.wsafe  = wsafe
    ifo.Suspension.LF              = LF
    ifo.Suspension.Stage[0].Mass   = mass

    ifo.Materials.MirrorMass       = mass
    ifo.Materials.MirrorVolume     = mass /ifo.Materials.Substrate.MassDensity
    ifo.Materials.MassRadius       = scale*ifo.Materials.MassRadius
    ifo.Materials.MassThickness    = scale*ifo.Materials.MassThickness

    if SQZ:
        if 'Squeezer' not in ifo:
            ifo.Squeezer = Struct()

        ifo.Squeezer.Type        = 'Freq Independent'
        ifo.Squeezer.AmplitudedB = SQZ
        ifo.Squeezer.lambsq1     = 0.05   #between squeezer and filter cavity
        ifo.Squeezer.lambsq2     = RTL  #inside the filter cavity
        ifo.Squeezer.lambsq3     = 0.1   #between filter cavity and SRM

        if FC:
            if 'FilterCavity' not in ifo.Squeezer:
                ifo.Squeezer.FilterCavity = Struct()

            ifo.Squeezer.Type           = 'Freq Dependent'
            ifo.Squeezer.FilterCavity.L = FC
            ifo.Squeezer.lambsq1     = 0.05   #between squeezer and filter cavity
            ifo.Squeezer.lambsq2     = RTL  #inside the filter cavity
            ifo.Squeezer.lambsq3     = 0.1   #between filter cavity and SRM
            ifo.Squeezer.LOAngleRMS  = 30e-3  #quadrature noise [radians]
            tbp = ifo.Squeezer.lambsq2*1e6
            print('FC losses:       %.2f ppm' % tbp)
            tbp = ifo.Squeezer.LOAngleRMS*1e3
            print('RMS phase noise: %.2f mrad' % tbp)
            print('FC length:       %.2f m' % ifo.Squeezer.FilterCavity.L)

        # parameters for shotrad (original GWINC)
        ifo.Squeezer.SQZAngle = 0
        ifo.Squeezer.InjectionLoss = ifo.Squeezer.lambsq1
        if FC:
            ifo.Squeezer.InjectionLoss += ifo.Squeezer.lambsq3

    return ifo


def precompPower(ifo):
    #
    #Followings are imported from gwinc.precomp.precompPower
    #pbs, parm, finesse, prfactor, Tpr = precompPower(ifo, PRfixed)
    #
    c       = scipy.constants.c
    t1      = sqrt(ifo.Optics.ITM.Transmittance)
    r1      = sqrt(1 - ifo.Optics.ITM.Transmittance)
    r2      = sqrt(1 - ifo.Optics.ETM.Transmittance)
    loss    = ifo.Optics.Loss                          # single TM loss
    #
    # Finesse, effective number of bounces in cavity, power recycling factor
    finesse = 2*pi / (t1**2 + 2*loss)        # arm cavity finesse
    neff    = 2 * finesse / pi               # effective number of bounces
    #
    # Arm cavity reflectivity with finite loss
    garm = t1/(1-r1*r2*sqrt(1-2*loss))  # amplitude gain wrt input field
    rarm = r1   -t1*r2*sqrt(1-2*loss) *garm
    #
    pbs  = ifo.gwinc.pbs
    parm = pbs * garm**2 / 2       # arm power from BS power
    #
    ifo.gwinc.parm    = parm
    ifo.gwinc.finesse = finesse
    #ifo.gwinc.prfactor = prfactor
    #ifo.Optics.PRM.Transmittance = Tpr
    #


def precompFiber(ifo):
    ##############################
    # calc TM fiber suspension radius

    wsafe = ifo.Suspension.Sapphire.wsafe
    LF    = ifo.Suspension.LF
    mass  = ifo.Materials.MirrorMass

    #radius of wire (m)
    r1 = ifo.Suspension.Stage[2].WireRadius #wire1: between IM and MN
    l1 = ifo.Suspension.Stage[2].Length     #wire1
    m1 = ifo.Suspension.Stage[2].Mass       #IM

    #wire2 between blade spring and TM; scale with mirror mass
    # (wsafe for safety factor; default value calculated
    #  with modulus of rupture to be 350 MPa (360e6*pi*r2**2/g*4/22.8))
    r2 = 1.6e-3/2*(mass/22.8)**(1./2)*(wsafe/12.5765)**(1./2)

    if LF:
        m1 = 20.5*4     # for KN KAGRA-LF
        r1 = 0.2e-3/2   # for KN KAGRA-LF
#        r2 = 3.2e-4/2   # for KN KAGRA-LF
        l1 = 26.1e-2*3  #wire1 for KN KAGRA-LF

    # Update ifo parameters
    ifo.Suspension.Stage[2].Mass       = m1 #IM
    ifo.Suspension.Stage[2].Length     = l1 #wire1
    ifo.Suspension.Stage[2].WireRadius = r1 #wire1: between IM and MN
    ifo.Suspension.Stage[0].WireRadius = r2 #wire2: blade spring and TM


def precompFDSqz(ifo):
    ##############################
    # calc optimal filter cavity parameters for FD squeezing

    if 'Squeezer' not in ifo:
#        print('caccaculo')
        return
    if 'FilterCavity' not in ifo.Squeezer:
#        print('caccaculo2')
        return

    ifo.Squeezer.FilterCavity.fdetune = 0
    ifo.Squeezer.FilterCavity.Ti      = 0

    sqzType = ifo.Squeezer.get('Type', 'Freq Independent')

    if not sqzType == 'Freq Dependent':
        return

    c       = scipy.constants.c                  # SOL [m/s]
    L       = ifo.Infrastructure.Length          # Length of arm cavities [m]
    lamb    = ifo.Laser.Wavelength               # Laser Wavelength [m]
    I0      = ifo.gwinc.pbs                      # power at BS
    mass    = ifo.Materials.MirrorMass           # Mirror mass [kg]
    trans   = ifo.Optics.ITM.Transmittance       # ITM Transmittance [Power]
    tau     = sqrt(ifo.Optics.SRM.Transmittance) # SRM Transmittance [amp.]
    rSRM    = sqrt(1 - tau**2)                   # SRM Reflectivity  [amp.]
    gamma   = trans*c/(4*L)                      # cavity pole (Hz)
    lambsq2 = ifo.Squeezer.lambsq2
    Lfilter = ifo.Squeezer.FilterCavity.L

    # for filter cavity parameters, see JGW-T1809537
    gammaRSE = (1+rSRM)/(1-rSRM)*gamma	# RSE linewidth
    I0RSE    = (1+rSRM)/(1-rSRM)*I0  # effective input power at BS for RSE (W)

    # SQL reaching frequency for BRSE (rad/s)
    omegaSQL = sqrt((-gammaRSE**2+sqrt(gammaRSE**4 \
                                     +64*pi*c*I0RSE/(mass*lamb*L**2)))/2)
    # loss parameter
    epsilonfilter = 4/(2+sqrt(2+2*sqrt(1+(4*omegaSQL*Lfilter/(c*lambsq2))**4)))

    if epsilonfilter == 1.0:
	# to avoid zero division in kappafilter when omegaSQL is too small
    	epsilonfilter = 1.0-1e-6

    kappafilter = sqrt(2/(2-epsilonfilter)/sqrt(1-epsilonfilter)) \
                *omegaSQL/sqrt(2)  # line width of filter cavity (rad/s)

    #detuning of filter cavity (Hz)
    detune = sqrt(1-epsilonfilter)*kappafilter
    Ttot   = kappafilter*4*Lfilter/c #transmittance of filter cavity
    Tfc    = Ttot-lambsq2            # filter cavity input mirror transmittance

    # detune = 0
    # Tfc = lambsq2

    ifo.Squeezer.FilterCavity.fdetune = -detune/(2*pi)
    ifo.Squeezer.FilterCavity.Ti      = Tfc

    # Parameters for sqzFilterCavity (original GWINC)
    ifo.Squeezer.FilterCavity.Te      = 0
    ifo.Squeezer.FilterCavity.Rot     = 0
    ifo.Squeezer.FilterCavity.Lrt     = ifo.Squeezer.lambsq2
    
    print('FC half bandwidth: %.2f Hz' % (kappafilter/(2*pi)))
    print('FC detuning:       %.2f Hz' % (detune/(2*pi)))
    print('FC input mirror transmissivity: %.4f' % (Tfc*100))


def temp_to_power(ifo):
    ## DETERMINATION OF LASER POWER AT BS ##########
    #Power at BS is determined from sapphre mirror temperature 
    # and np.absorption.

    mass = ifo.Materials.MirrorMass

    #radius of mirror (m);    scale with respect to mass (fixed aspect ratio)
    #radius = 11e-2*(mass/22.8)**(1./3)

    #thickness of mirror (m); scale with respect to mass (fixed aspect ratio)
    height = 15e-2*(mass/22.8)**(1./3)

    LF     = ifo.Suspension.LF
    trans  = ifo.Optics.ITM.Transmittance    #ITM Transmittance [Power]
    Tm     = ifo.Materials.Substrate.Temp    #temperature at TM [k]
    T2     = ifo.Suspension.Temp             #temperature at blade spring [K]
    r2     = ifo.Suspension.Stage[0].WireRadius #wire2: blade spring and TM
    l3     = ifo.Suspension.Stage[0].Length     #wire2 [m]
    N3     = ifo.Suspension.Stage[0].NWires

    subloss = 50e-6*1e2       #energy loss rate of substrate (/m)
    #see JGW-G1706429-v1 (2017)
    coaloss = 0.5e-6          #energy loss rate of coating
    outside = 0.05            #radiation energy from outside (W)
    if LF:
        outside = 0.003	      # for KN KAGRA-LF
    #see CQG 31, 224003 (2014)
    
    alphak = 7.98*(r2/0.8e-3)	# included r2-dependance on March 29, 2018
    betak  = 2.2      #kappa=alphak*T**betak ,thermal conductivity (W/K/m)
    #see CQG 31, 105004 (2014)
    
    #total heatflow
    Kabs   = N3*alphak/(betak+1)*pi*r2**2/l3*(Tm**(betak+1)-T2**(betak+1))
    
    #power between BS and ITM
    Pmich  = (Kabs-outside)/(2*subloss*height+coaloss*4/trans)

    I0attenuation = ifo.Laser.I0attenuation

    I0     = 2*Pmich           #power at BS
    I0     = I0*I0attenuation  #attenuated from maximum power
    if  I0 < 0:
    	I0 = 1e-3              # make it very small if minus
    #print "Power at BS: %d W" % (I0)

    ifo.gwinc.pbs = I0
