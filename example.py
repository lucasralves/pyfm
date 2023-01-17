from src.airfoil_mesh import airfoil_mesh

if __name__ == '__main__':

    mesh = airfoil_mesh(ns=100, nt=5, nf=30, delta=0.1, exp1=1.1, exp2=1.1, exp3=1.1, ext_radius=10.0, ext_cell_size=1.0)
    mesh.build_surface(file='./data/foils/foil-3.txt')
    mesh.show_surface()
    mesh.build_mesh('./data/meshes/mesh-1.msh', gmsh_view=True)