# pyfm

Python script that creates a mesh around an airfoil to be used in CFD analysis. The mesh is generated using the Gmsh api.

## How to install:

To install the package from pip, use the code below:

### Linux
```python
pip3 install git+https://https://github.com/lucasralves/pyfm
```

### Windows (not tested)
```python
pip3 install -e git+https://https://github.com/lucasralves/pyfm#egg=pyfm
```

## How to use it:

```python
import pyfm

file = 'file.txt'         # file containing the airfoil points
ns = 100                  # number of points on the profile surface
nt = 5                    # number of points on the trailing edge
nf = 30                   # number of layers in the boundary layer
delta = 0.1               # boundary layer height
exp1 = 1.1                # expansion ratio of points on the surface towards the leading edge
exp2 = 1.1                # expansion ratio of points on the surface towards the trailing edge
exp3 = 1.1                # boundary layer expansion ratio
ext_radius = 10.0         # outer countor radius
ext_cell_size = 1.0       # size of the element on the outer counter
out_file = 'out_file.msh' # output file
gmsh_view = True

pyfm.init(file, ns, nt, nf, delta, exp1, exp2, exp3, ext_radius, ext_cell_size)
pyfm.build_surface()
pyfm.show_surface()
pyfm.build_mesh(out_file, gmsh_view=gmsh_view)
```

## Mesh

The surface mesh, created by the function 'pyfm.buid_surface(...)', can be seen by calling 'pyfm.show_surface()'.

<p float="left">
  <img src="./doc/surface_mesh.png" width="45%" />
  <img src="./doc/surface_mesh_zoom.png" width="45%" />
</p>

The volume mesh contain two regions, one close to the wall made of quadrangular elements that comprises the boundary layer, and a second region made of triangular elements.

<img src="./doc/mesh_zoom_2.png" width="45%" />

<img src="./doc/mesh_full.png" width="45%" />

The triangular region is divided into four parts, in order to maintain a refined region closer to the airfoil.

<img src="./doc/mesh_zoom_1.png" width="45%" />

There are 4 physical groups, foil, external, laterals and volume, that can be used to set the boundary conditions.