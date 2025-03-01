# GED_2024

Research on ILPs for Graph Edit Distance

## Installing

### Following software is needed:

- GUROBI
- GEDLIB https://github.com/dbblumenthal/gedlib
- LIBLSAP: https://forge.greyc.fr/projects/liblsap/ (header only, add LIBLSAP_ROOT=~/liblsap to cmake options)
need to change ```static const DataType _zero = 0;``` in liblsap/cpp/include/lsap.hh, line 50 to ```static constexpr DataType _zero = 0;```

- (maybe install Gurobi as well)

### Variables to add to path/cmake:

#### add the following to your cmake options:

- -DGUROBI_DIR=/path/to/gurobi1200/linux64/

#### add the following to your environment:

- LIBLSAP_ROOT=/path/to/liblsap;GEDLIB_ROOT=/path/to/gedlib;GUROBI_HOME="/path/to/gurobi1200/linux64/"
- for a screenshot of how my CLion configuration is set up look at "cmake_settings.png"
