import numpy as np
from quasistatic import Coil, ObsPoints

debug = False

if __name__ == '__main__':
    line = Coil.from_csv('line_2pts_1m.csv')

    point = ObsPoints(np.array([[.1,0,.6]]))
    plane = ObsPoints.from_csv('plane_X_500_Y_500_10mm.csv')

    results = line.compute(point)
    print(results)

    # fig, ax = plt.subplots(1,1)
    # c1 = ax.scatter(plane.points[:,0], plane.points[:,1], c=np.linalg.norm(results,axis=1), vmin=0, vmax=1e-6)
    # fig.colorbar(c1)
    # plt.show()