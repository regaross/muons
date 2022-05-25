#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'

'''
muons.py


Contact:
Regan Ross
rross@laurentian.ca

A module for simulating muons underneath SNOLAB's overburden for nEXO's Outer Detector.

'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################

import numpy as np
from time import time
#from stl import mesh

### Seed for random number generation from system time
SEED = np.random.seed(time)


#################################################
#                   CONSTANTS                   }
#             Referenced Throughout             }
#################################################

# COORDINATE SYSTEM

ORIGIN = (0,0,0) # m

# nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 13.333      # m
OD_CENTER = ORIGIN      # m         defines coordinate system with respect to literal centre of OD
OC_RADIUS = 2.770       # m
OC_POSITION = (0,0,0)   # m         positions OC with respect to OD centre

# MUON FLUX PARAMETERS AT SNOLAB
SNOLAB_MU_FLUX = 3.31e-10       # \pm (0.01 (stat) \pm 0.09 (sys))e-10 mu/cm^2/s # arXiv:0902.2776v1
SNOLAB_MU_E_AVG = 363.0         # \pm 1.2 GeV # arXiv:1909.11728v1
SNOLAB_DEPTH = 5.890    #km.w.e        # \pm 94 km.w.e.  # arXiv:1909.11728v1


#################################################
#                    CLASSES                     }
#                                                }
#################################################


class OuterCryostat:
    '''
    A class for the spherical outer cryostat within nEXO's outer detector
        - radius                Sphere radius [meters]
        - center      *         Sphere center [meters, meters, meters]

        * Center positions the outer cryostat w.r.t the outer detector center.

        '''

    def __init__(self, radius=OC_RADIUS, position=OC_POSITION) -> 'OuterCryostat':
        '''A basic constructor; see class docstring'''
        self.radius = radius
        self.position = position


class OuterDetector:
    '''
    A class for nEXO's *cylindrical* outer detector.
        - radius                Cylinder radius [meters]
        - height                Cylinder height [meters]
        - center      *         Cylinder center [meters, meters, meters]
        - vertical depth        Attenuation by overburden [km.w.e]
        - fill_height   **      Vertical level to which tank is filled with water [meters]

    * Subordinate components will be positioned w.r.t the center which is the literal center of the cylinder,
    NOT the center of the bottom of the cylinder.

    ** fill_height specifies to which vertical level the cylinder is filled with liquid. This number is used to
        define the COVER GAS region.

    '''

    def __init__(self, radius=OD_RADIUS, height=OD_HEIGHT, center=OD_CENTER, vertical_depth=SNOLAB_DEPTH, fill_height=OD_HEIGHT-0.2) -> 'OuterDetector':
        ''' A basic constructor; see class docstring'''
        self.radius = radius
        self.height = height
        self.center = center
        self.vertical_depth = vertical_depth
        self.fill_height = fill_height

class Muon:
    '''
    A class for muons parameterized by empirical functions.

        - zenith                [rad]
        - azimuth               [rad]
        - energy                [GeV]
        - initial_position      [meters, meters, meters]
        - 
    '''

    def __init__(self, zenith=0, azimuth=0, energy=SNOLAB_MU_E_AVG, initial=(0,0,OD_HEIGHT)) -> 'Muon':
        ''' A constructor for the muon. Defaults to vertical muon at average SNOLAB energy'''

        self.zenith = zenith
        self.azimuth = azimuth
        self.energy = energy
        self.initial = initial


    ### Instance Functions

    def get_unit_vec(self) -> tuple:
        ''' Returns a unit vector for the muon direction (x,y,z)'''

        x = np.sin(self.zenith)*np.cos(self.azimuth)
        y = np.sin(self.zenith)*np.sin(self.azimuth)
        z = np.cos(self.zenith)

        return (x,y,z)

    def get_cartesian_track(self) -> tuple:
        ''' Returns the muon track (x, y, z, x0, y0, z0)'''

        return self.get_unit_vec()+self.initial