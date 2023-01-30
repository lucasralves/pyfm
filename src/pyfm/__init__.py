from pyfm.utils.airfoil_mesh_impl import _AirfoilMeshImpl


def init(file: str = None,
         ns: int = None, nt: int = None, nf: int = None,
         delta: float = None,
         exp1: float = None, exp2: float = None, exp3: float = None,
         ext_radius: float = None, ext_cell_size: float = None) -> None:
    """
    pyfm.init(...)

    Initialize global parameters. This must be called before any call to the other
    functions. If the parameters are None, they must be given when pyfm.build_surface(...)
    and pyfm.build_mesh(...) are called.

    Parameters
    ----------
    - file: airfoil points.
    - ns: number of points on the profile surface.
    - nt: number of points on the trailing edge.
    - nf: number of layers in the boundary layer.
    - delta: boundary layer height.
    - exp1: expansion ratio of points on the surface towards the leading edge.
    - exp2: expansion ratio of points on the surface towards the trailing edge.
    - exp3: boundary layer expansion ratio.
    - ext_radius: outer countor radius.
    - ext_cell_size: size of the element on the outer counter.
    """
    global mesh
    mesh = _AirfoilMeshImpl(file=file, ns=ns, nt=nt, nf=nf, delta=delta, exp1=exp1, exp2=exp2, exp3=exp3, ext_radius=ext_radius, ext_cell_size=ext_cell_size)
    return

def build_surface(file: str = None,
                  ns: int = None, nt: int = None,
                  exp1: float = None, exp2: float = None) -> None:
    """
    pyfm.build_surface(...)
    
    Divide the surface into ns elements.

    Parameters
    ----------
    - ns: number of points on the profile surface.
    - nt: number of points on the trailing edge.
    - exp1: expansion ratio of points on the surface towards the leading edge.
    - exp2: expansion ratio of points on the surface towards the trailing edge.
    """
    mesh.build_surface(file=file, ns=ns, nt=nt, exp1=exp1, exp2=exp2)
    return

def show_surface() -> None:
    """
    pyfm.show_surface()
    
    Show a plot of the surface with the grid points.
    """
    mesh.show_surface()
    return

def build_mesh(out_file: str,
               gmsh_view: bool = False,
               nf: int = None,
               delta: float = None,
               exp3: float = None,
               ext_radius: float = None, ext_cell_size: float = None) -> None:
    """
    pyfm.build_mesh(...)
    
    Create the mesh.

    Parameters
    ----------
    -out_file:
    - nf: number of layers in the boundary layer.
    - delta: boundary layer height.
    - exp3: boundary layer expansion ratio.
    - ext_radius: outer countor radius.
    - ext_cell_size: size of the element on the outer counter.
    """

    mesh.build_mesh(out_file, gmsh_view=gmsh_view, nf=nf, delta=delta, exp3=exp3, ext_radius=ext_radius, ext_cell_size=ext_cell_size)
    return
