This folder contains some patches that are required for compiling some of the required libraries in Ubuntu 18.04 and cmake.

Issues:

* To be able to compile old releases of CGAL (<4.7, I think...) a small change is needed in Installation/src/CMakeLists.txt. You need to find and replace "//CMakeLists.txt" by "/CMakeLists.txt" (one slash less!). Solution found here: https://github.com/vigente/gerardus/issues/107
* We also need to set the gcc/g++ 4.8 version in boost configuration files.
