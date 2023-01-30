from abc import ABC, abstractclassmethod

##################################################################
##################################################################

class AirfoilMeshAbs(ABC):
    """
    Creates a mesh around and airfoil.

    Parameters
    ----------
    - ns: number of points on the profile surface.
    - nt: number of points on the trailing edge.
    - nf: number of layers in the boundary layer.
    - delta: boundary layer height.
    - exp1: expansion ratio of points on the surface towards the leading edge.
    - exp2: expansion ratio of points on the surface towards the trailing edge.
    - exp3: boundary layer expansion ratio.
    - ext_radius: outer countor radius.
    - cell_ratio: ratio between the size of the element on the outer counter and on the
      surface of the boundary layer.
    """
    
    def __init__(self, file: str = None, ns: int = None, nt: int = None, nf: int = None, delta: float = None, exp1: float = None, exp2: float = None, exp3: float = None, ext_radius: float = None, cell_ratio: float = None) -> None:
        ...

    @abstractclassmethod
    def build_surface(self, file: str = None, ns: int = None, nt: int = None, exp1: float = None, exp2: float = None) -> None:
        """Divide the surface into ns elements"""
        ...
    
    @abstractclassmethod
    def show_surface(self) -> None:
        """Show a plot of the surface with the grid points"""
        ...
    
    @abstractclassmethod
    def build_mesh(self) -> None:
        """Create the mesh"""
        ...