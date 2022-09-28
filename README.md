# muons

##  muon_functions.py
 This is a module for simulating muons passing through nEXO's outer detector (OD hereafter). It comprises a few classes, each with particular attributes and instance methods. While the tools within are specifically for determining intersection points, path lengths and fluxes for nEXO, they are written to be generalizeable for any cylindrical underground apparatus.
 
### Classes
#### Muon
The muon class is the basis for this module. Each generated muon is an instance of this class. The price paid in memory and processing speed is made up for in clarity. One can call upon the attributes of a muon's track without needing to memorize array indices or coordinate transformations. It also makes the code easy to change and simplifies the addition of further functions.

##### Attributes
| Attribute | Name | Description | Type |
|-----------|------|-------------|------|
| Zenith angle |  ``` zenith ```| Muon track zenith angle w.r.t vertical (radians) | float |
| Azimuthal angle| ``` azimuth ```| Muon track azimuthal angle w.r.t x-axis (radians) | float |
| Energy | ```energy``` | Muon energy (GeV) | float |
| Initial Position | ```initial``` | Muon instantiation position in Cartesian coord's (meters) | tuple |

##### Instance Functions
The following are the headers and descriptions of the muon class instance functions. 

```python
def get_unit_vec(self) -> tuple:
```

This function returns a tuple unit vector for the muon in the format ($x_\mu$  $y_\mu$ $z_\mu$).  The elements sum to 1 in quadrature.

```python
def get_cartesian_track(self) -> tuple:
```

This function returns a tuple of the muon's <em>track</em> given by its unit vector elements followed by the coordinates of its initial position. Namely, (x, y, z, x0, y0, z0). 

```python
def closest_approach(self, point) -> float:
```

This function returns the distance of closest approach of the muon to the argument ```point```. Where ```point``` is an iterable object containing 3 cartesian coordinates. 


#### Outer Detector
The Outer Detector class is a class meant to encapsulate the relevant details of nEXO's Outer Detector. Those relevant details at present are few, but must be mutable. While the OD is in planning stages, dimensions may change. It's easier to create a new instance of a class than to re-write a large number of hard-coded values. The OD here is modeled as a basic cylinder filled to a certain level with water. Its axis is oriented with the vertical.

##### Attributes:
| Attribute | Name | Description | Type |
|-----------|------|-------------|------|
| Radius |  ``` radius ```| The radius of the OD cylinder (meters) | float |
| Height | ``` height ```| The height of the OD cylinder (meters) | float |
| Depth Underground | ```vertical_depth``` | The depth of the detector underground in attenuation units (km.w.e) | float |
| Position | ```center``` | The position of OD centre with respect to the coordinate system | tuple |
| Water Fill Level |```fill_height``` | The distance from the bottom of the OD to the fill level (meters) | float |
| Outer Cryostat | ```outer_cryo``` | An instance of the Outer Cryostat Class | OuterCryostat|

The only instance function for this class at present is the constructor.

#### OuterCryostat
This class is basically just a sphere with a position in space. It is used to subtract from muon path lengths- no Cherenkov light would be produced in this volume within the Outer Detector.

##### Attributes 
| Attribute | Name | Description  | Type |
|----------|-------|--------------|------|
| Radius | ```radius``` | The radius of the spherical OC (meters) | float |
| Center | ``center`` | The center of the OC with respect to the OD centre (meters) | tuple |

Once more, the only function for this class is the constructor.

### Functions
This module makes use of a number of equations, some empirical, some phenomenological, some parametric. There are a few functions whose initial purpose was the same, but only one is used by default. The remaining unused functions are left because they may still be useful to a future user of this code. Where possible, the functions were written to allow for array use. That is, the functions work just as well for a `numpy` array full of values as they do for one single value. This makes the code faster for larger jobs.

The following is a listing of the function headers and a brief description of their arguments, outputs and use cases.

```python 
def mei_hime_intensity(zenith_angles, vert_depth = SNOLAB_DEPTH) -> np.ndarray:
``` 

This function is an equation from Mei & Hime's paper: [Muon Induced Background Study for Underground Laboratories ](https://arxiv.org/abs/astro-ph/0512125). More specifically, it is the zenith angle intensity relation- exactly the same as Eqn. 3 from the paper. The arguments are an array of zenith angles and the vertical depth at which to evaluate the muon intensity from the angles. 

**Returns** a numpy array of intensities [cm$^{-2}$ sr$^{-1}$ s$^{-1}$  ]

```python
def mei_hime_normed_discrete(zenith_angles, vert_depth = SNOLAB_DEPTH) -> np.ndarray:
```

This function is an equation from Mei & Hime's paper: [Muon Induced Background Study for Underground Laboratories ](https://arxiv.org/abs/astro-ph/0512125). More specifically, it is the zenith angle intensity relation **projected** onto the horizontal and **normalized** discretely (as opposed to continuously by integration). The arguments are an array of zenith angles and the vertical depth at which to evaluate the muon intensity from the angles. 

**Returns** a numpy array of zenith angle probabilities 


```python
def mei_hime_normed_continuous(zenith_angles, vert_depth = SNOLAB_DEPTH) -> np.ndarray:
```

This function is an equation from Mei & Hime's paper: [Muon Induced Background Study for Underground Laboratories ](https://arxiv.org/abs/astro-ph/0512125). More specifically, it is the zenith angle intensity relation **projected** onto the horizontal and **normalized continously** by integrating the function numerically over the zenith angle range $\theta \in [0, \pi/2)$. The arguments are an array of zenith angles and the vertical depth at which to evaluate the muon intensity from the angles. 

**Returns** a numpy array of zenith angle probabilities 

```python
def mh_energy_probs(energies, zenith = 0, sample = True) -> np.ndarray:
```

This function uses an equation from Mei & Hime's paper: [Muon Induced Background Study for Underground Laboratories ](https://arxiv.org/abs/astro-ph/0512125). More specifically, it uses equation 8. It takes an array of energies and returns a normalized array of relative probabilities; one for each energy. If `sample = True` the array is discretely normalized so that the returned array sums to 1. Otherwise, it is normalized by integrating the function (which is more useful for plotting). 

**Returns** a numpy array of energy probabilities

```python
def gaisser_normed_discrete(energies, zenith) -> np.ndarray:
```

This function uses an equation from Thomas Gaisser's book [Cosmic Rays and Particle Physics](https://doi.org/10.1017/CBO9781139192194). It parameterizes the intensity of high-energy muons in the atmosphere. This function determines the relative probabilities of energies using Gaisser's equation and the argument zenith angle. The equation is normalized discretely such that the sum of all elements becomes 1. It can be deployed to sample high energy muons for propogation underground. It is currently used in the `generate_muons()` function to sample energies.

**Returns** a numpy array of energy probabilities

```python
def generate_muons(how_many, outer_detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
```

This function does most of the work in this entire module. It uses the `gaisser_normed_discrete` and `mei_hime_normed_discrete` functions to sample muon energy and zenith angle respectively. It uses the provided OuterDetector argument to determine the maximum generation radius based on a criterion outlined in [this document](#). If `gen_radius` and `gen_offset` are provided as arguments, they override this default behaviour based on the OuterDetector object. If they're left to each be 0 by default, the generation parameters will scale with the OuterDetector.

**Returns** a numpy array of Muon objects

```python
def muons_per_square_meter(rate, outer_detector, z_offset = 0, gen_radius = 0):
```

This function calls `generate_muons()` to generate an array of muons having an initial constant areal flux on the generating disk. The argument `rate` is the number of muons to be generated per meter squared of the calculated (or given) disk size. The array of muons returned will be of size `rate`$\times A$ where $A$ is the area of the generating disk.

**Returns** a numpy array of Muons

```python
def intersection_points(muon, outer_detector = OuterDetector(), labels = True, tolerance = 0.001) -> tuple:
```

This function takes as arguments *a single Muon* object, an OuterDetector, a boolean, `labels`, and a tolerance level. The muon is converted into a parametric line in 3-space where its intersection points with the argument OuterDetector are determined by analytically solving parametric equations. `labels` , if true, will return a label with each entry or exit point. Each label is either 'TOP', 'BOT', 'SIDE' corresponding to where the muon enters and exits the OD. 

Also of note, if no entry point is determined, the function doesn't bother looking for an exit point (even though there may be one due to tiny floating point error) and just returns `(False, False)` as the entry and exit points respectively. 

If the muon does intersect the OD, the function returns the entry and exit points, and  `if labels`, the tuple takes the form ((x1,y1,z1), entry_label, (x2,y2,z2), exit_label)

**Returns** a tuple of tuple points, OR (False, False) OR tuple of tuple points with labels

```python
def path_length(muon, outer_detector = OuterDetector(), labels=True, ignore_cover_gas=True, ignore_cryostat=True) -> float:
```

This function uses `intersection_points()` to determine, by the Pythagorean theorem, and some parametric equations (if cyrostat is True) to determine the path length of a muon through the OD. `if ignore_cryostat == True`, the path length is the Pythagorean distance minus the distance traversed through the cryostat. The same goes for the cover gas. A muon would not produce light visible to PMTs within the cryostat or the cover gas.

**Returns** `False` if the muon doesn't intersect, a float distance if it does.

```python
def hits_detector(muon, outer_detector) -> bool:
```

This function is pretty self-explanatory. It checks whether or not the muon hits the detector using the `intersection_points()` function.

**Returns** a boolean; `True` if it hits, `False` if it doesn't

```python
def intersecting_muons(how_many, outer_detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
```

This function uses `generate_muons()` to generate an array of muons but discards those that do not intersect the OuterDetector. It uses `hits_detector()` to determine whether or not to discard the muon in question.

**Returns** a numpy array of Muons

```python
def path_through_covergas(muon, outer_detector) -> float:
```

This function determines the path of an intersecting muon through the cover gas portion of the OuterDetector. No Cherenkov light is produced in this region, albeit small, but this function determines that path length.

**Returns** a float of this path length 

```python
def path_through_cryostat(muon, outer_detector) -> float:
```

This function determines the path length through the cryostat of `outer_detector` of a muon intersecting `outer_detector` . This is done by solving a series of parametric equations.

**Returns** a float of this path length

```python
def get_cherenkov(muon, outer_detector = OuterDetector(), photons_per_meter = False) -> np.ndarray:
```

This function returns the Cherenkov yield of photons (per meter if `photons_per_meter = True`) and the angle at which they're emitted for a muon based on its path length through the OD. It automatically removes the path length through the cryostat and the cover gas. The function deployes the Frank Tamm Formula and the relativistic kinetic energy to calculate these values.

**Returns** a numpy array: [\<number of photons>, \<angle of emission>]

## generate_muons.py
This is a simple python script that uses `muon_functions.py` to write an array of muons to a csv file. It pulls parameters from `muon_params.yaml`
