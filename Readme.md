# Splats Surface Reconstruction

Legacy project including all the surface reconstruction from unorganized points methods presented in my [thesis](https://www.tesisenred.net/handle/10803/145380) (and a some more).

This repository is kept as reference for what was implemented at the time.

## Requirements

This project depends on the following libraries (tested with these SPECIFIC versions):

* [CGAL 4.0.2](https://doc.cgal.org/Manual/4.0.2/doc_html/cgal_manual/packages.html)
* [Boost 1.52](https://www.boost.org/users/history/version_1_52_0.html)
* [levmar 2.6](http://users.ics.forth.gr/~lourakis/levmar/)
* [ARPACK/ARPACK++](https://www.caam.rice.edu/software/ARPACK/)
  
Also, all sources must be compiled with gcc/g++ version **4.8** (may work with other versions, but has not been tested, and definetely not working with the latest ones...).

## Methods

The methods/executables provided in this project are:

* `mesh_point_set`: direct meshing of point sets using a Delaunay refinement mesher. It implements the method presented in:
  > R. Campos, R. Garcia, P. Alliez, and M. Yvinec. A Surface Reconstruction Method for In-Detail Underwater 3D Optical Mapping. International Journal of Robotics Research, 34(1):64-89, 2015. eISSN: 1741-3176, ISSN: 0278-3649.

* `compute_splats`: computes the splats structure. The rest of methods work with this kind of structure as input, so executing this method is the first step of the rest of methods in the project (see below).

* `mesh_splats`: meshes a set of splats using a Delaunay refinement mesher. Directly uses the splats to compute the intersection with the query segments required by the mesher. Due to the approximations performed, depending on the level of noise in the splats the resulting surface may be non-manifold (cleaning steps are also provided to minimize this issue), and the global orientation of the mesh may be flipped (we do not care what is the "inside" or the "outside" of the object, so we cannot orient it consequently). In conjuntion with `compute_splats` they implement the method presented in:
  > Ricard Campos, Rafael Garcia, Pierre Alliez, and Mariette Yvinec. Splat-based surface reconstruction from defect-laden point sets. Graphical Models, Volume 75, Issue 6, November 2013, Pages 346-361, ISSN 1524-0703

* `orient_splats`: the splats computed with `compute_splats` do not take into account the orientation of the resulting local surfaces. This method tries to orient them coherently by propagating an arbitrary seed orientation following a Minimum Spanning Tree.

* `mesh_oriented_splats`: given a set of oriented splats (output of `orient_splats`), it creates a Signed Distance Function (SDF). This SDF is then meshed by extracting the zero iso-surface via a Delaunay refinement mesher. 
  Note that, while not being published anywhere, using `orient_splats` + `mesh_oriented_splats` closely resembles the approach followed in the [classic article by H. Hoppe](https://hhoppe.com/proj/recon/), but allowing for higher-order local surfaces other than *tangent planes*:
  > Hugues Hoppe, Tony DeRose, Tom Duchamp, John McDonald, and Werner Stuetzle. 1992. Surface reconstruction from unorganized points. SIGGRAPH Comput. Graph. 26, 2 (July 1992), 71–78. https://doi.org/10.1145/142920.134011

* `mesh_splats_gcut`/`mesh_splats_ncut`: both methods create an unsigned distance function (UDF) from the splats discretely evaluated on an adaptive tetrahedral grid. In order to transform the UDF to a SDF, so as to be able to extract the zero isosurface, both methods try to partition the tetrahedra in the grid containing the UDF into being *inside*/*outside* of the object. Both methods mainly differ in the graph partitioning method used, being Graph Cuts for `mesh_splats_gcut` and Normalize Cuts for `mesh_splats_ncut`. 
  As a main difference between them, `mesh_splats_gcut` requires the object to be "watertight", while `mesh_splats_ncut` does not. Both methods are detailed in the following paper:
  > R. Campos and R. Garcia. Surface meshing of underwater maps from highly defective point sets. Journal of Field Robotics 35 (4), 491-515.

The parameters of each method can be obtained by calling the executables with the `--help` option. For instance:

```
compute_splats --help
```

Aside from the brief description provided in the help messages, the parameters of each method are not documented. However, you can find a detailed description and usage for most of them in all the articles listed above.

## Docker

Assuming you are in the folder containing this file, to build the image from the Dockerfile provided in this project you should run:

```
docker build -t <the_image_tag> .
```

Note that the Dockerfile contains a detailed description of the configuration and installation steps required for compiling this project in Ubuntu 18.04.

We provide an already-compiled docker image at [DockerHub](https://hub.docker.com/r/ricardcd/splats_surface_reconstruction), so that you can use the tools in this project without having to compile it yourself.

In order to run the tools with your point sets, just mount the directory containing the data when running the image:

```
docker run -it -v <your_local_dir>:<mount_dir_inside_container> docker pull ricardcd/splats_surface_reconstruction:latest
```

From within the image, you can run directly the tools listed in the previous section of this document, as they are already in the path. 
For instance, you can run the demos provided by:

* Downloading the docker image:

```
docker pull ricardcd/splats_surface_reconstruction:1.0
```

* Running a bash inside a container (note the `-v` option to mount a local drive to `./demo/results`, so you can access the results easily):

```
docker run -it -v <your_local_dir:/opt/splats_surface_reconstruction/demo/results ricardcd/splats_surface_reconstruction:1.0 /bin/bash
```

* And executing the demos script from within the container:

```
cd demo
./run_demos_sphere.sh
```

