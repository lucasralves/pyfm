import sys
sys.path.append('./src/')

import pyfm

if __name__ == '__main__':

    file = './examples/data/foils/foil-3.txt'
    ns = 100
    nt = 5
    nf = 30
    delta = 0.1
    exp1 = 1.1
    exp2 = 1.1
    exp3 = 1.1
    ext_radius = 10.0
    ext_cell_size = 1.0
    out_file = './examples/data/meshes/mesh-3.msh'
    gmsh_view = True

    pyfm.init(file, ns, nt, nf, delta, exp1, exp2, exp3, ext_radius, ext_cell_size)
    pyfm.build_surface()
    pyfm.show_surface()
    pyfm.build_mesh(out_file, gmsh_view=gmsh_view)