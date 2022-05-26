#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'

'''
muon_functions.py


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
import time
#from stl import mesh

### Seed for random number generation from system time
SEED = int(time.time())
np.random.seed(SEED)

#################################################
#                   CONSTANTS                   }
#             Referenced Throughout             }
#################################################


# nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 12.800      # m
OD_CENTER = (0,0,0)     # m         defines coordinate system with respect to literal centre of OD
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
        - outer_cryo            The outer cryostat within the detector - a sphere for path lengths

    * Subordinate components will be positioned w.r.t the center which is the literal center of the cylinder,
    NOT the center of the bottom of the cylinder.

    ** fill_height specifies to which vertical level the cylinder is filled with liquid. This number is used to
        define the COVER GAS region.

    '''

    def __init__(self, radius=OD_RADIUS, height=OD_HEIGHT, vertical_depth=SNOLAB_DEPTH,\
         fill_height=OD_HEIGHT-0.2, outer_cryo = OuterCryostat()) -> 'OuterDetector':
        ''' A basic constructor; see class docstring'''
        self.radius = radius
        self.height = height
        self.vertical_depth = vertical_depth
        self.fill_height = fill_height
        self.outer_cryo = outer_cryo

class Muon:
    '''
    A class for muons parameterized by empirical functions.

        - zenith                [rad]
        - azimuth               [rad]
        - energy                [GeV]
        - initial_position      [meters, meters, meters]
        - 
    '''

    def __init__(self, zenith=0, azimuth=0, energy=SNOLAB_MU_E_AVG, initial=(0,0,0)) -> 'Muon':
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

    def closest_approach(self, point):
        ''' Finds the muon's closest distance to the provided point in 3-space (cartesian coords)'''

        unit_vector = self.get_unit_vec()

        # Parameterize for simplicity

        i, j, k = point[0], point[1], point[2]
        mx, my, mz = unit_vector[0], unit_vector[1], unit_vector[2]
        x0, y0, z0 = self.initial[0], self.initial[1], self.initial[2]

        t = ((mx*i + my*j + mz*k) - (mx*x0 + my*y0 + mz*z0))/(mx**2 + my**2 + mz**2)

        muon_point = (mx*t + x0, my*t + y0, mz*t + z0)

        x_dist = muon_point[0] - i
        y_dist = muon_point[1] - j
        z_dist = muon_point[2] - k

        distance = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)

        return distance



#################################################
#                   FUNCTIONS                    }
#                                                }
#################################################


def mei_hime_intensity(zenith_angles, vert_depth = SNOLAB_DEPTH):
    ''' Function from Mei & Hime's zenith angle intensity relation- exactly the same as Eqn. 3 from the paper.'''
    # Parameters same as in paper (Equation No. 3)
    I1 = 8.60e-6 #  /sec/cm^2/sr
    I2 = 0.44e-6 #  /sec/cm^2/sr
    lam1 = 0.45 #  km.w.e.
    lam2 = 0.87 #  km.w.e.

    intensity = 2*np.pi*(I1*np.exp(-vert_depth/(lam1*np.cos(zenith_angles)))+I2*np.exp(-vert_depth/\
        (lam2*np.cos(zenith_angles))))/np.cos(zenith_angles)

    return intensity

def mei_hime_normed_discrete(zenith_angles, vert_depth = SNOLAB_DEPTH) -> np.ndarray:
    ''' Returns the normalized through-going muon flux through horizontal surface for an array of angles 
        "angles" at a specified vertical depth vert_depth in km.w.e. NORMALIZED for the given array of angles
        by making array sum to 1. Differential solid angle factor included.'''

    normed_flux = mei_hime_intensity(zenith_angles)* np.cos(zenith_angles)*np.sin(zenith_angles)  # Horizonal projection * differential solid angle

    #Normalizes the distribution function such that array sums to 1 for any number of elements
    normed_flux = normed_flux / np.sqrt(normed_flux.sum()**2)

    return normed_flux


def mei_hime_normed_continuous(zenith_angles, vert_depth = SNOLAB_DEPTH):
    ''' Returns the normalized through-going muon flux through horizontal surface 
        for an array of angles "angles" at a specified vertical depth vert_depth in km.w.e
        NORMALIZED based on the integral, not the discrete array'''

    norm_const = 6.006989507403272e-11 #Integrated from 0 to pi/2

    #Normalizes the distribution function
    normed_flux = mei_hime_intensity(zenith_angles)*np.cos(zenith_angles)*np.sin(zenith_angles)

    return normed_flux(zenith_angles)/norm_const


def gaisser_normed_discrete(energies, zenith):
    ''' Surface muon energy distribution based on Gaisser's formalism '''

    # Making Gaisser's distribution as a lambda function [cm^{-2}s^{-1}sr^{-1}GeV^{-1}]
    dNdE = lambda E, theta : 0.14*(E**(-2.7))*((1/(1+(1.1*E*np.cos(theta)/115)))+(0.054/(1+(1.1*E*np.cos(theta)/850))))

    energy_array = dNdE(energies, zenith)

    # Normalize it for use as PDF
    return energy_array/np.sum(energy_array)


def generate_muons(how_many, outer_detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
    ''' Generates an array of Muons (instances of class) given an OuterDetector cylinder
        - how_many          The number of muons to generate
        - outer_detector    The OuterDetector volume to be used in simulations
        - gen_radius        The radius of the concentric disk on which muons will be instantiated
        - gen_offset        The distance from the top of the outer_detector to the concentric disk
        '''

    muons = np.empty(how_many, dtype=object)

    sampling_size = int(10*how_many)

    # If the user has not set the generator offset and or radius
    if gen_offset == 0:
        gen_offset = outer_detector.height

    if gen_radius == 0:
            gen_radius = np.tan(1)*(outer_detector.height + gen_offset) + outer_detector.radius

    # Define muon initial positions
    rhos = np.random.random(size = how_many)*(gen_radius**2)
    gen_angles = np.random.random(size = how_many)*np.pi*2

    # Keeping in mind the coordinate transformations 
    initial_x = np.sqrt(rhos)*np.cos(gen_angles)
    initial_y = np.sqrt(rhos)*np.sin(gen_angles)
    initial_z = np.ones(how_many)*(gen_offset + outer_detector.height/2)

    # Determining zenith and azimuthal angles
    theta_radians = np.linspace(0, np.pi/2, sampling_size)
    zenith_probabilities = mei_hime_normed_discrete(theta_radians, outer_detector.vertical_depth)
    zeniths = np.random.choice(theta_radians, p = zenith_probabilities, size = how_many)
    
    azimuths = np.random.random(size = how_many)*np.pi*2

    # Selecting Energies
    surface_energy_range = np.linspace(6500, 100000, sampling_size)
    surface_energies = np.random.choice(surface_energy_range, p = gaisser_normed_discrete(surface_energy_range, 0), size = how_many)

    # Attenuate based on zenith angles
    b = 0.4 #km.w.e.^-1
    energies_underground = (surface_energies - SNOLAB_MU_E_AVG*(np.exp(b*np.cos(zeniths)*SNOLAB_DEPTH)-1))/np.exp(b*np.cos(zeniths)*SNOLAB_DEPTH)
    for e in range(len(energies_underground)): 
        if energies_underground[e] < 0: energies_underground[e] = 0
    
    
    # Instantiate muons
    for i in range(how_many):
        muons[i] = Muon(zeniths[i], azimuths[i], energies_underground[i], (initial_x[i], initial_y[i], initial_z[i]))

    return muons

def intersection_points(muon, outer_detector = OuterDetector(), labels = True, tolerance = 0.001):
    ''' A function for analytically determining the intersection points of a muon with an outer detector cylinder.  '''
    
    entryPoint, exitPoint = False, False
    entryLabel, exitLabel = '',''

    detRadius = outer_detector.radius
    detHeight = outer_detector.height
    #det_z_translation = detector.position[2]

    #We can parametrize the muon for simplicity:
    zenith, azimuth = muon.zenith, muon.azimuth
    mx = np.sin(zenith)*np.cos(azimuth)
    my = np.sin(zenith)*np.sin(azimuth)
    mz = -np.cos(zenith)        # Does this need to be negative?

    x0, y0, z0 = muon.initial[0], muon.initial[1], muon.initial[2] 

    # Quadratic equation parameters
    a = (mx**2 + my**2)
    b = 2*(mx*x0 + my*y0)
    c = (x0**2 + y0**2 - detRadius**2)
    det_squared = (b**2 - 4*a*c)
    
    # Initial setting
    qhighCheck = detRadius + 1
    qlowCheck = qhighCheck
    
    if det_squared > 0:
        qlow = float(-b/(2*a) - np.sqrt(det_squared)/(2*a))
        zlow = mz*qlow + z0
        qlowCheck = a*qlow**2 + b*qlow + (x0**2 + y0**2)
        qhigh = float(-b/(2*a) + np.sqrt(det_squared)/(2*a))
        zhigh = mz*qhigh + z0
        qhighCheck = a*qhigh**2 + b*qhigh + (x0**2 + y0**2)


    # Find the first point
    ptop = (detHeight/2 - z0)/mz
    xtop = mx*ptop + x0
    ytop = my*ptop + y0

    # ENTRY POINT POSSIBILITIES
    entryLabel = ''
    if xtop**2 + ytop**2 <= detRadius**2:
        # Hits top
        entryPoint = (xtop, ytop, detHeight/2)
        entryLabel = 'TOP'
    elif qhighCheck < (detRadius + tolerance)**2 and qhighCheck > (detRadius - tolerance)**2 and zhigh**2 < (detHeight/2)**2:
        # Hits the side at the higher z point
        entryPoint = (mx*qhigh + x0, my*qhigh + y0, zhigh)
        entryLabel = 'SIDE'
        qhigh = False
    elif qlowCheck < (detRadius + tolerance)**2 and qlowCheck > (detRadius - tolerance)**2 and zlow**2 < (detHeight/2)**2:
        # Hits the side at the lower z point
        entryPoint = (mx*qlow + x0, my*qlow + y0, zlow)
        entryLabel = 'SIDE'
        qlow = False

    if type(entryPoint) is tuple: #If it does actually enter the cylinder
        # Bottom Point Parameters
        pbottom = (-detHeight/2 - z0)/mz
        xbottom = mx*pbottom + x0
        ybottom = my*pbottom + y0

        exitLabel = ''
        # EXIT POINT POSSIBILITIES
        if xbottom**2 + ybottom**2 <= detRadius**2:
            # Hits the bottom of the cylinder
            exitPoint = (xbottom, ybottom, -detHeight/2)
            exitLabel = 'BOT'

        elif qhighCheck < (detRadius + tolerance)**2 and qhighCheck > (detRadius - tolerance)**2 and zhigh**2 < (detHeight/2)**2 and type(qhigh) is float:
            exitPoint = (mx*qhigh + x0, my*qhigh + y0, zhigh)
            exitLabel = 'SIDE'

        elif qlowCheck < (detRadius + tolerance)**2 and qlowCheck > (detRadius - tolerance)**2 and zlow**2 < (detHeight/2)**2 and type(qlow) is float:
            exitPoint = (mx*qlow + x0, my*qlow + y0, zlow)
            exitLabel = 'SIDE'
    
    if type(entryPoint) is bool:
        return False

    elif labels:
        return (entryPoint, entryLabel, exitPoint, exitLabel)
    else:
        return (entryPoint, exitPoint)


def path_length(muon, outer_detector = OuterDetector(), cryostat = False):
    ''' Returns the path length of a muon through the Outer Detector. False if it doesn't hit.'''

    points = intersection_points(muon, outer_detector, labels = False)
    path_length = 0

    if type(points) is not bool:
        x = points[1][0] - points[0][0]
        y = points[1][1] - points[0][1]
        z = points[1][2] - points[0][2]

        path_length = np.sqrt(x**2 + y**2 + z**2)

        return path_length

    else:
        return False

#################################################
#                   PLOTTING                     }
#                  FUNCTIONS                     }
#################################################

def plot_path_lengths(points, savefile = False):
    ''' A function designed to plot a histogram of the points from the previous intersection_points function. Pass this function an array
    argument equivalent to the output from the previous.
    '''
    import matplotlib.pyplot as plt
    plt.rcParams.update({
    "text.usetex": True,})

    pathLengths = []
    top_bottom = [] #Muons that pass through both top and bottom
    top_side = []
    side_side = []
    ssCount = 0
    side_bottom = []

    ##(entryPoint, entryLabel, exitPoint, exitLabel)
    for point in points:
        if not type(point) is bool: #If there is an entry point
            if type(point[1]) is str:
                pathLength = np.sqrt((point[0][0]-point[2][0])**2 + (point[0][1] - point[2][1])**2 + (point[0][2] - point[2][2])**2)
                #If points are both top and bottom, append there
                if point[1] == 'TOP' and point[3] == 'BOT':
                    top_bottom.append(pathLength)
                #If points are top side, append there
                elif point[1] == 'TOP' and point[3] == 'SIDE':
                    top_side.append(pathLength)
                #If points are side side, append there
                elif point[1] == 'SIDE' and point[3] == 'SIDE':
                    side_side.append(pathLength)
                    ssCount += 1
                #If points are side bottom, append there
                elif point[1] == 'SIDE' and point[3] == 'BOT':
                    side_bottom.append(pathLength)
                pathLengths.append(pathLength)
            else:
                print(' Points input for path_lengths function requires labels being switched on ')

    bins = 100
    counts_tb, bins_tb = np.histogram(top_bottom, bins = bins, density=False)
    counts_ts, bins_ts = np.histogram(top_side, bins = bins, density=False)
    counts_ss, bins_ss = np.histogram(side_side, bins = bins, density=False)
    counts_sb, bins_sb = np.histogram(side_bottom, bins = bins, density=False)

    counts, bins = np.histogram(pathLengths, bins = 50, density=False)
    plt.figure()
    plt.hist(bins[:-1], bins, weights=counts, histtype='stepfilled', alpha=0.4, color='orange', label = 'Total')
    plt.xlabel('Path Length [m]', size = 'large'); plt.ylabel('Count', size = 'large')
    #plt.text(1, counts[0]*1.5, 'Mean = ' + str(np.average(pathLengths)), size = 12)
    plt.title('Path Length Distribution: '+ str(len(pathLengths)) + ' paths', size = 'x-large')

    #Plotting other subordinate hists
    top_bot = str(len(top_bottom)*100/len(pathLengths))[0:4]
    plt.hist(bins_tb[:-1], bins, weights=counts_tb, histtype='step', alpha=1.0, color='blue', label = r'Top $\rightarrow$ Bottom: ' + top_bot + '\%')
    top_sides = str(len(top_side)*100/len(pathLengths))[0:4]
    plt.hist(bins_ts[:-1], bins, weights=counts_ts, histtype='step', alpha=1.0, color='green', label = r'Top $\rightarrow$ Side: ' + top_sides + '\%')
    sides = str(len(side_side)*100/len(pathLengths))[0:4]
    plt.hist(bins_ss[:-1], bins, weights=counts_ss, histtype='step', alpha=1.0, color='red', label = r'Side $\rightarrow$ Side: ' + sides + '\%')
    side_bot = str(len(side_bottom)*100/len(pathLengths))[0:4]
    plt.hist(bins_sb[:-1], bins, weights=counts_sb, histtype='step', alpha=1.0, color='purple', label = r'Side $\rightarrow$ Bottom: ' + side_bot + '\%')

    plt.yscale('log')
    plt.legend(loc = 8, fontsize = 'large')
    #plt.grid()
    #print('Mean:', np.mean(pathLengths))
    if savefile:
        plt.savefig('pathLengths.png', facecolor = 'white')

    plt.show()