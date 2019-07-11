# Vertex Extrapolator module readme

Yorck Ramachers (Warwick)
Last updated July 11, 2019

The Vertex Extrapolator module is a SuperNEMO reconstruction module. It attempts to use the TTD data bank in Falaise in the current
event data model to extract trajectory data and extrapolate the solutions to create vertices. 

## Files:

- vertex_module.cpp
- vertex_module.h
- vertex_library.cpp
- vertex_library.h
- CMakeLists.txt
- fit.conf.in


## Description

Add to an flreconstruct pipeline to fit clustered tracker hits for reconstruction. To build it, do

``` console
$ ls
CMakeLists.txt  fit.conf vertex_module.cpp  vertex_module.h  README.md  testing

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

The build will create the `libVertex.so` shared library. Assuming that you have an `input.brio` file that contains a `TCD` bank (i.e. clustered data) from the reconstruction, this can be run after editing the configuration file to 
point at the library location. No module configuration parameter are required.:

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


## Data model
