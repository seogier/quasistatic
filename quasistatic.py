import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0, pi
import pyqtgraph as pg

class Grid:
    """Rectilinear grid class
    Grid axes are along X, Y, and Z, and units are in meters
    Grid doesn't have to have more than one point in each axis
    Grid points don't have to be evenly spaced.
    """
    def __init__(self, Xs, Ys, Zs) -> None:
        """Initialize Grid
        X, Y, and Z are lists of points along those axes"""
        self.Xs = Xs
        self.Ys = Ys
        self.Zs = Zs
    
    @property
    def points(self):
        """Get list of (x, y, z) points
        Result is actually 2D, but the second dimension is the XYZ components
        [nX*nY*nZ, 3]"""
        X, Y, Z = np.meshgrid(self.Xs,self.Ys,self.Zs, indexing='ij')
        points = np.stack([X.flatten().T, Y.flatten().T, Z.flatten().T],1)
        return points

    @property
    def shape(self):
        """Get shape of grid"""
        return [len(self.Xs), len(self.Ys), len(self.Zs)]
    
    @classmethod
    def evenly_spaced_3D(cls, X_range=0, dX=0, Y_range=0, dY=0, Z_range=0, dZ=0):
        """Generate grid of points evenly spaced in X, Y, Z
        X/Y/Z_range are [start, stop] for each axis
        dx/dY/dZ are step size for each axis
        
        If X, Y, or Z are single values or only 1 long
        OR
        if dX, dY, or dZ == 0
        The the first (or only) value given for X, Y, or Z will be the only
        value in that dimension for that grid.
        """
        # Make inputs into lists if they're not indexible
        if type(X_range) not in (list, set, tuple):
            X_range = [X_range]

        if type(Y_range) not in (list, set, tuple):
            Y_range = [Y_range]
        
        if type(Z_range) not in (list, set, tuple):
            Z_range = [Z_range]

        # If input is of length 1 or dX == 0, X is X_range[0]
        # Similarly for Y and Z
        if (len(X_range) == 1) or (dX == 0):
            Xs = [X_range[0]]
        else:
            Xs = np.arange(X_range[0], X_range[1]+dX, dX)

        if (len(Y_range) == 1) or (dY == 0):
            Ys = [Y_range[0]]
        else:
            Ys = np.arange(Y_range[0], Y_range[1]+dY, dY)

        if (len(Z_range) == 1) or (dZ == 0):
            Zs = [Z_range[0]]
        else:
            Zs = np.arange(Z_range[0], Z_range[1]+dZ, dZ)

        return cls(Xs, Ys, Zs)
    
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(*self.points.T)
 
class Result(Grid):
    def __init__(self, grid, values) -> None:
        Grid.__init__(self, grid.Xs, grid.Ys, grid.Zs)
        self.values = values
    
    @classmethod
    def from_1D(cls, grid: Grid, values):
        """Reshape [nX*nY*nZ, 3] list of points to [nX, nY, nZ, 3]
        and return a results objkect."""
        
        nX = len(grid.Xs)
        nY = len(grid.Ys)
        nZ = len(grid.Zs)
        shape = (nX, nY, nZ, 3)
        # print('New Shape:', shape)
        values = values.reshape(shape)
        return cls(grid, values)
        
    @property
    def B1_mag(self):
        """Calculate magnitude of XY component of values"""
        return np.linalg.norm(self.values[..., :-1], axis=-1)
    
    @property
    def B1_phase(self):
        """Calculate phase of B1"""
        return np.arctan2(self.values[...,1], self.values[...,0])
    
    def view_B1_mag(self):
        p = pg.plot(self.B1_mag.T,axes={'t':1,'x':0,'y':2})
        p.view.invertY(False)

class Coil:
    """MRI coil for quasistatic analysis.

    Coil is loaded from specified CSV file
    """

    def __init__(self, points) -> None:
        self.points = points

        # Compute individual line segment vectors
        self.l_vecs = self.points[1:,:]-self.points[0:-1,:] # vectors between l_n and l_n+1
        self.l_mags = np.linalg.norm(self.l_vecs, axis=1) # magnitude of l vectors
        self.l_hats = self.l_vecs/self.l_mags[:,None] # unit l vectors
    
    @classmethod
    def from_csv(cls, points_csv):
        # Load coil points from csv file
        points = np.genfromtxt(points_csv, delimiter=',', skip_header=1)
        return cls(points)

    def plot_coil(self, ax = None):
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1, projection='3d')
        ax.plot(*self.points.T)
        ax.scatter(*self.points.T)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_aspect('equal')
        ax.view_init(elev=172, azim=-28, roll=90)
        if ax == None:
            fig.tight_layout()
        plt.show()
        return ax
    
    def compute(self, grid: Grid) -> Result:
        obs_points = grid.points
        result_array = np.empty_like(obs_points)
        for i in range(obs_points.shape[0]):
            P = obs_points[i,:] # Observation Point
            l0s = self.points[0:-1,:]
            l1s = self.points[1:,:]

            # Note: einsum('ij,ij->i', A, B) is the row-wise dot product between corresponding rows in A and B
            A = P-l0s
            Proj = np.einsum('ij,ij->i', A, self.l_hats)
            Rs = (A) - Proj[:,None]*self.l_hats # Shortest vector from line to observation point
            R_mags = np.linalg.norm(Rs, axis=1)
            R_hats = Rs/R_mags[:,None]

            cos_theta_0 = np.einsum('ij,ij->i', P-l0s, l1s-l0s)/(np.linalg.norm(P-l0s, axis=1)*np.linalg.norm(l1s-l0s, axis=1))
            cos_theta_1 = np.einsum('ij,ij->i', P-l1s, l0s-l1s)/(np.linalg.norm(P-l1s, axis=1)*np.linalg.norm(l0s-l1s, axis=1))

            # print(f'R:     {Rs[0,:]}')
            # print(f'R mag: {R_mags[0]}')
            # print(f'R hat: {R_hats[0,:]}')
            # print(f'cos(theta_zero): {cos_theta_0[0]}')
            # print(f'cos(theta_one):  {cos_theta_1[0]}')

            phis = np.cross(self.l_hats, R_hats)
            summand = phis * (mu_0/(4*pi*R_mags)*(cos_theta_0+cos_theta_1))[:,None]
            result_array[i,:] = np.sum(summand,0)
        # print(result_array.shape)
        result = Result.from_1D(grid, result_array)
        # print(result.B1_mag.shape)
        return result
