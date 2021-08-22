# OpenCSLSolver

<img src="http://csl.etu.edu.tr/wp-content/uploads/2021/01/CSL_Logo.png" width="200">

OpenCSLSolver is an open-source, cross-platform, parallel, pressure based CFD solver. OpenCSLSolver is a side project of me. OpenCSLSolver is refactored version of some parts from the solver I developed for my [M.Sc. thesis](https://www.researchgate.net/publication/332719238_Development_of_a_Pressure_Based_Unstructured_GPU_Accelerated_CFD_Solver_for_Compressible_Reacting_Flows_at_all_MACH_Numbers) for [Combustion Systems Laboratory](http://csl.etu.edu.tr/) at [TOBB University of Economics and Technology](https://www.etu.edu.tr/).

Features
* Segregated Pressure Based
* 3D unstructured mesh using [.su2](https://su2code.github.io/docs/Mesh-File/) format (1D and 2D can be simulated using EMPTY BC)
* [Paraview](https://www.paraview.org/) VTK file generation for post processing
* [OpenMP](https://www.openmp.org/) parallel
* Fine-grained parallel CG and ILUBICGSTAB matrix solvers
* Developed using modern C++ techniques
* No dependencies, everyting is pure C++17, Standart Library and OpenMP

## Quickstart
Windows binary executable and validation cases can be found in last stable release on this project's GitHub page. Just download the executable. To run Shocktube test case, open Command Window at this folder and run this command:
```bash
OpenCSLSolver.exe config_shocktube.txt
```
This binaries compiled without OpenMP so parallel running is not possible. If you want to run parallel simulations from optimized executable, build project from source.

## Building From Source
This project is using [CMAKE](https://cmake.org/) for building and installation. Just build and install this project using CMAKE in Windows, Linux or macOS.

## Usage
Using config_template.txt create a configuration file for your simulation. A 1D mesh can be generated inside OpenCSLSolver or a .su2 mesh file can be read. Simulation can be run from terminal or command window:
```bash
<path to solver executable> <path to configuration file>
```
Or you can use meshes and configuration file from validation cases from the last stable release on this project's GitHub page.

## Version Control
There will be no version numbering. Latest stable release with compilation date can be found on release page.

## Support
Please open an issue to get help and support about this project.

## Roadmap
K-epsilon turbulence model and simple chemistry solver will be integrated in near future. Code documentation, user guide and theoretical guide is planned.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Acknowledgements
This solver is possible thanks to Dr. Sıtkı USLU. And I want to thank all current and past researchers of Combustion Systems Laboratory for their help and contribution. Also I want to thank CMPS development team for sharing their knowledge. And I want to thank Peren Turan for her patience and support.

## License

Distributed under the [GPLv3](http://www.gnu.org/licenses/) License. See `LICENSE.txt` for more information.
