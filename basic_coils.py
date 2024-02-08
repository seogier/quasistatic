import numpy as np
import argparse

def save_csv(points, fname):
    np.savetxt(fname, points, header='X,Y,Z', delimiter=',')

def circ_loop(radius: float, segments: int, fname=None, Z: float =0):
    """Generate and save a circular loop in the XY plane at offset Z
    """
    if fname == None:
        fname = f'circ_loop_{radius*2:.1f}mm.csv'
    
    angles = np.linspace(0, 2*np.pi, segments+1)
    
    points = np.zeros((segments+1,3))
    for i, angle in enumerate(angles):
        points[i,:] = [radius*np.cos(angle), radius*np.sin(angle), Z]
    
    return points
    # save_csv(points, fname)

def line(length: float, segments: int, fname=None):
    """Generate and save a line of current along the Z axis
    """
    if fname == None:
        fname = f'line_{length:.1f}mm.csv'

    length_m = length/1000

    points = np.zeros((segments+1,3))
    points[:,2] = np.linspace(-length_m/2, length_m/2, segments+1)

    save_csv(points, fname)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate basic coils')
    parser.add_argument('shape', type=str, help='Shape of coil to generate', choices=['circ_loop','line'])
    parser.add_argument('save_name', type=str, nargs='?', help='Filename to save as', default=None)
    parser.add_argument('-r', '--radius', type=float, help='Radius of loop in mm')
    parser.add_argument('-l', '--length', type=float, help='Length of line in mm')
    parser.add_argument('-s', '--segments', type=int, help='Number of segments to discritize coil')
    parser.add_argument('-z', '--z-offset', type=float, help = 'Z offset of coil in XY plane', default=0.0)

    args = parser.parse_args()
    shape = args.shape
    fname = args.save_name
    radius = args.radius
    segments = args.segments
    z_offset=  args.z_offset
    length = args.length

    if shape == 'circ_loop':
        circ_loop(radius, segments, fname, z_offset)
    elif shape == 'line':
        line(length, segments, fname)
    