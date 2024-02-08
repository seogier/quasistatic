# Quasi-magnetostatic EM Solver

This solver uses the Biot-Savart Law to implement a quasi-magnetostatic solver for the magnetic fields produced by a wire structure.

This solver makes the following assumptions:
 - The coil is a single infintessimally thin wire path
 - The current is a uniform 1 A
 - All distance units are in meters

# File Formats
## Coil
The coil is a CSV file of (x,y,z) coordinates with the following format:

| X   | Y   | Z   |
|-----|-----|-----|
| x0  | y0  | z0  |
| x1  | y1  | z1  |
| ... | ... | ... |
| xN  | yN  | zN  |

Each wire segment is a straight line between (x_n, y_n, z_n) and (x_n+1, y_n+1, z_n+1).

The coil definition makes the present assumptions:
 - The coil has a constant current of 1.0 A
 - The coil is a single continuous wire from the (x0, y0, z0) to (xN, yN, zN)

If a closed loop is desired, the last point must be the same as the first.

## Observation Points
The obeservation points are defined in a file with the following format:

| X  | Y   | Z  |
|----|-----|----|
| x0 | y0  | z0 |
| x1 | y1  | z1 |
| ...| ... | ...|
| xN | yN  | zN |

## Results
The results for a given coil and set of observation points are retured in a file with the following format:
| X   | Y   | Z   | Bx  | By  | Bz  |
|-----|-----|-----|-----|-----|-----|
| x0  | y0  | z0  | Bx0 | By0 | Bz0 |
| x1  | y1  | z1  | Bx1 | By1 | Bz1 |
| ... | ... | ... | ... | ... | ... |
| xN  | yN  | zN  | BxN | ByN | BzN |