# OpenFOAM repository

OpenFOAM is the open source software I used during my PhD to perform Computational Fluid Dynamic (CFD) simulations.
It is written in C++ and a clone of the repo can be found in the [OpenFOAM-dev](OpenFOAM-dev) directory.

In our research group at Politecnico di Milano, we worked on a fork of that library that is not public as it is currently being used by several companies. As we worked on the development version of OpenFOAM, where changes are pushed on a daily basis, part of my work was to ensure that our fork was up to date with the official public version.

PoliMi's library mainly focuses on extensions to deal with mesh motion based on topology changes, adding and removing cells where it is necessary. This involves dealing with the core of the library, as the number of cells in the domain changes and, therefore, the algebraic system of equations that needs to be solved. 

In terms of code, we mainly modified the following parts of the code:
* [finiteVolume](OpenFOAM-dev/src/finiteVolume/): core of the solver
* [fvMeshDistributors](OpenFOAM-dev/src/fvMeshDistributors/): for load balancing on parallel computations
* [fvMeshMovers](OpenFOAM-dev/src/fvMeshMovers/): for mesh motion
* [fvMeshStitchers](OpenFOAM-dev/src/fvMeshStitchers/): to correct mesh fluxes
* [fvMeshTopoChangers](OpenFOAM-dev/src/fvMeshTopoChangers/): to dynamically refine or un-refine the mesh
* [meshTools](OpenFOAM-dev/src/meshTools/): added some tools used to change the topology of the mesh
* [polyTopoChange](OpenFOAM-dev/src/polyTopoChange/): to modify the topology of the mesh

In addition, for the second part of my PhD I also implemented a novel hybrid spray model, modelling the fuel injection inside a combustion chamber. This part was an enhancement of the [lagrangian](OpenFOAM-dev/src/lagrangian/) standard code.

In the [dani-dev/src](dani-dev/src/) directory you can find an example of some code I did during my PhD. This code extends [Adaptive Mesh Refinement (AMR)](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) to work with 1, 2 or 3 dimensional geometries.