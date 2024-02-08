import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0, pi

debug = False

class ObsPoints:
    def __init__(self, points) -> None:
        self.points = points
    
    @classmethod
    def from_csv(cls, points_csv):
        points = np.genfromtxt(points_csv, delimiter=',', skip_header=1)
        return cls(points)

    @classmethod
    def plane(cls, width: tuple[float, float], step: float, orientation: str='XY', offset: float =0):
        """Generate a rectangular grid in the XY plane at offset Z

        width is (X, Y) dimensions in mm
        step is sampling interval in mm
        orientation is orientation of plane XY, XZ, YZ
        offset is offset in mm
        """
        width_m = [w/1000 for w in width]
        offset_m = offset/1000
        step_m = step/1000
        
        if orientation == 'XY':
            x = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
            y = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)
            z = np.array(offset_m)
        elif orientation == 'XZ':
            x = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
            y = np.array(offset_m)
            z = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)
        elif orientation == 'YZ':
            x = np.array(offset_m)
            y = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
            z = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)


        X, Y, Z = np.meshgrid(x,y,z)
        points = np.stack([X.ravel().T, Y.ravel().T, Z.ravel().T],1)

        return cls(points)

    @classmethod
    def cube(cls, width: tuple[float, float, float], step: float):
        """Generate a cubic rectilinear grid

        width is (X, Y, Z) dimensions in mm
        step is sampling interval in mm
        """
        width_m = [w/1000 for w in width]
        step_m = step/1000
        
        x = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
        y = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)
        z = np.arange(-width_m[2]/2, width_m[2]/2+step_m, step_m)

        X, Y, Z = np.meshgrid(x,y,z)
        points = np.stack([X.ravel().T, Y.ravel().T, Z.ravel().T],1)

        return cls(points)
    
    @classmethod
    def line(start, stop, points):
        """Generate a line of evenly spaced points from start to stop

        width is (X, Y, Z) dimensions in mm
        step is sampling interval in mm
        """
        start_m = [s/1000 for s in start]
        stop_m = [s/1000 for s in stop]
        X = np.linspace(start_m[0], stop_m[0], points)
        Y = np.linspace(start_m[1], stop_m[1], points)
        Z = np.linspace(start_m[2], stop_m[2], points)
        points = np.stack([X,Y,Z], axis=1)
        return points

    def plot_points(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(*self.points.T)

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

    def plot_coil(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        ax.plot(*self.points.T)
        ax.scatter(*self.points.T)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.show()
        return ax
    
    def compute(self, obs: ObsPoints):
        result = np.empty_like(obs.points)
        for i in range(obs.points.shape[0]):
            P = obs.points[i,:] # Observation Point
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

            if debug:
                print(f'R:     {Rs[0,:]}')
                print(f'R mag: {R_mags[0]}')
                print(f'R hat: {R_hats[0,:]}')
                print(f'cos(theta_zero): {cos_theta_0[0]}')
                print(f'cos(theta_one):  {cos_theta_1[0]}')

            phis = np.cross(self.l_hats, R_hats)
            summand = phis * (mu_0/(4*pi*R_mags)*(cos_theta_0+cos_theta_1))[:,None]
            result[i,:] = np.sum(summand,0)
        return result
