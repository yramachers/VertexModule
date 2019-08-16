# Vertex Extrapolator module readme

Yorck Ramachers (Warwick)
Last updated August 13, 2019

The Vertex Extrapolator module is a SuperNEMO reconstruction module. It ought to use the TTD data bank in Falaise in the current
event data model to extract trajectory data and extrapolate the solutions to create vertices. Without fit errors in the TTD bank,
this inevitably does not work. This is hence not(!) a working Falaise module.

## Files:

- vertex_module.cpp
- vertex_module.h
- vertex_library.cpp
- vertex_library.h
- CMakeLists.txt
- ve.conf.in


## Description

Add to an flreconstruct pipeline to extrapolate trajectories to the foil and calorimeter walls in order to obtain vertex areas. To build it, do

``` console
$ ls
CMakeLists.txt  ve.conf.in vertex_library.cpp vertex_library.h vertex_module.cpp  vertex_module.h  README.md  testing

$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=$(brew --prefix) ..
...
$ make
...
... If you are developing the module, you can test it by calling
$ make test
...
... or obtain more detail on the tests and launch in the build directory
$ ctest -V
```

Note: if you get a QT5 error, you may need to specify the QT5 path when you run the cmake line, as given by `brew --prefix qt5-base`. For example, you can run:
``` console
$ cmake -DCMAKE_PREFIX_PATH="$(brew --prefix qt5-base);$(brew --prefix)" ..
``` 

The build will create the `libVertex.so` shared library. Assuming that you have an `input.brio` file that contains a `TTD` bank (i.e. tracker trajectories) from the reconstruction, this could be run after editing the configuration file to point at the library location, see above for why this is not possible at the moment. No module configuration parameter are required.:

``` console
...
[name="flreconstruct.plugins" type="flreconstruct::section"]
plugins : string[1] = "Vertex"
Vertex.directory : string = "/my/path/to/build/directory"

# Define the modules in the pipeline:
[name="pipeline" type="vertex_module"]
...
$ flreconstruct -i /path/to/input.brio -p ve.conf 
```

## Implementation
Given a set of fitted models, lines, helices and broken line models, this module obtains the border planes of the inside of the calorimeter walls and the foil and extrapolates the models to those planes. That results in 
intersection points which form a vertex region. Here the region is assumed conservatively to be a rectangle as opposed to an ellipse. This simplifies the vertex area geometry and results from the assumption of error independence of 
perpendicular fit model parameters. The latter is considered conservative, increasing the vertex area, i.e. the places where the true vertex may be, to a rectangle that includes the expected error ellipse.  

The difficulty with vertex extrapolation is less the intersection calculation of a model with a plane but the possibility that the vertex area may touch up to three calorimeter planes, if for instance the best fit points towards a 
box corner. In that case, the extrapolated error variations of the model may well touch different planes to the left and right and again other planes higher and lower than the centre intersection point. All these areas are valid 
vertex areas even if they are on different calorimeter planes. That complicates matters a little.  

Another item of consideration is the fundamental difference between line models and helix models. The former can by definition only show a unique intersection point with a valid calorimeter wall. All other walls will be pierced 
outside the defined tracker limits. Helices, however, can have up to four permissible intersections. Finding the correct one is hence an additional process. 

On wire vertices, with no possibility to single out preferred regions in a cell as the origin of a vertex on a wire, the entire cell is declared as a vertex region in form of a box with sides the predefined +-2.2 cm and height the 
z-resolution. It is not clear how useful such a definition is but this should cover all cases and refinements can be implemented as required.

For subsequent modules, these vertex areas should be considered as characertising a 2D PDF of potential intersection points from fitting a model to a cluster. The area is spanned by varying model parameters by a single sigma error 
around the best fit parameter, independent from other fit parameters. A calorimeter association module should take that into account when taking the decision to associate a calo hit with an overlapping vertex area in order to 
declare a particle solution.  


## Data model
Rectangle structure:
- the two axes through the centre point (as Interval objects)
- a third axis, normally empty, to allow definition of a box region (Interval object)
- the fraction of the total area (double) for this rectangle vertex area
- the id of the plane this area sits on (int)
- which side of the tracker (int)

Interval structure (helper object):
- holds two doubles with useful methods.
