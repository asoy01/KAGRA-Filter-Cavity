
from .precomp import precomp_kagra
from .precomp import params_kagra

from . import noise_kagra as noise

from collections import OrderedDict

def noise_calc(f, ifo):
    """Calculate all IFO noises and return as dict

    Assumes ifo has already been run through precompIFO().

    """
    noises = OrderedDict()

    noises['Quantum Vacuum']       = noise.quantum.shotrad             (f, ifo)
    noises['Suspension Thermal']   = noise.suspensionthermal.susptherm (f, ifo)
    noises['Coating Brownian']     = noise.coatingthermal.coatbrownian (f, ifo)
   #noises['Coating Thermo-Optic'] = noise.coatingthermal.thermooptic  (f, ifo)
    noises['Substrate Thermo-Elastic']= noise.substratethermal.subtherm(f, ifo)
    noises['Substrate Brownian']   = noise.substratethermal.subbrownian(f, ifo)
    noises['Seismic']              = noise.seismic.seismic             (f, ifo)
    noises['Newtonian Gravity']    = noise.newtonian.gravg             (f, ifo)
   #noises['Excess Gas']  = noise.residualgas.gas(f, ifo)

    noises['Total'] = sum(noises[curve] for curve in noises)
    noises['Freq']  = f

    return noises

# def gwinc_kagra(freq, ifoin, conf = 'DRSE', verbose = False):
def gwinc_kagra(freq, ifoin, conf = 'DRSE', verbose = False, SQZ=0,RTL=0):
    """Calculate strain noise budget for a specified interferometer model.

    Argument `freq` is the frequency array for which the noises will
    be calculated, and `ifoin` is the IFO model (see the `load_ifo()`
    function).

    Returns tuple of (noises, ifo)

    """
    ifo = params_kagra(ifoin, conf, SQZ, RTL)

    title = 'None'
    if verbose:
        title = 'bKAGRA ('+conf+')' if conf.find('RSE') >= 0 \
           else 'KAGRA+ ('+conf+')'

    # add some precomputed info to the ifo struct
    #this implicitly deepcopies and the return value is the copy
    ifo = precomp_kagra(freq, ifo, title)

    noises = noise_calc(freq, ifo)

    return noises #, ifo
