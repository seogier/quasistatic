import numpy as np
import argparse

def save_csv(points, fname):
    np.savetxt(fname, points, header='X,Y,Z', delimiter=',')

def plane_XY(width: tuple[float, float], step: float, fname=None, Z_offset: float =0):
    """Generate and save a rectangular grid in the XY plane at offset Z

    width is (X, Y) dimensions in mm
    step is sampling interval in mm
    Z is Z offset in mm
    """
    width_m = [w/1000 for w in width]
    Z_offset_m = Z_offset/1000
    step_m = step/1000
    
    if fname == None:
        fname = f'plane_X_{width[0]:.0f}_Y_{width[1]:.0f}_{step:.0f}mm.csv'
    
    x = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
    y = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)
    z = np.array(Z_offset_m)

    X, Y, Z = np.meshgrid(x,y,z)
    points = np.stack([X.ravel().T, Y.ravel().T, Z.ravel().T],1)
        
    save_csv(points, fname)

def cube(width: tuple[float, float, float], step: float, fname=None):
    """Generate and save a rectilinear grid in cube

    width is (X, Y, Z) dimensions in mm
    step is sampling interval in mm
    """
    width_m = [w/1000 for w in width]
    step_m = step/1000

    if fname == None:
        fname = f'plane_X_{width[0]:.0f}_Y_{width[1]:.0f}_Z_{width[2]:.0f}_{step:.0f}mm.csv'
    
    x = np.arange(-width_m[0]/2, width_m[0]/2+step_m, step_m)
    y = np.arange(-width_m[1]/2, width_m[1]/2+step_m, step_m)
    z = np.arange(-width_m[2]/2, width_m[2]/2+step_m, step_m)

    X, Y, Z = np.meshgrid(x,y,z)
    points = np.stack([X.ravel().T, Y.ravel().T, Z.ravel().T],1)
        
    save_csv(points, fname)

def line(start, stop, points, fname=None):
    # start_m = [s/1000 for s in start]
    # stop_m = [s/1000 for s in stop]
    if fname == None:
        fname = f'line.csv'
    X = np.linspace(start[0], stop[0], points)
    Y = np.linspace(start[1], stop[1], points)
    Z = np.linspace(start[2], stop[2], points)
    points = np.stack([X,Y,Z], axis=1)
    return points

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate basic observation point arrays')
    parser.add_argument('pattern', type=str, help='Layout of observation points to generate.', choices=['plane_XY', 'cube', 'line'])
    parser.add_argument('save_name', type=str, nargs='?', help='Filename to save as', default=None)
    parser.add_argument('-x', '--x-width', type=float, help='Width in X in mm')
    parser.add_argument('-y', '--y-width', type=float, help='Width in Y in mm')
    parser.add_argument('-z', '--z-width', type=float, help='Width in Z in mm')
    parser.add_argument('-zo', '--z-offset', type=float, help = 'Z offset of sampling points in XY plane', default=0.0)
    parser.add_argument('-s', '--step', type=float, help='Sampling step size in mm')

    args = parser.parse_args()
    shape = args.pattern
    fname = args.save_name
    x_width = args.x_width
    y_width = args.y_width
    z_width = args.z_width
    z_offset=  args.z_offset
    step = args.step

    if shape == 'plane_XY':
        plane_XY((x_width, y_width), step, fname, z_offset)
    elif shape == 'cube':
        cube((x_width, y_width, z_width), step, fname)
    