# @Author: Tanushree Gupta<tshreegupta@gmail.com>
# @Date:   2021-11-29T19:11:22-05:00
# @Last modified by:   tshree
# @Last modified time: 2021-11-29T19:33:47-05:00

import numpy as np
from voxel import write

## Creates Phantom ball with N*N*N size
# Params:
#   1)fname: output file to save Phantom
#   2) N: length of 1 dimension of 3-D phantom.
def create_ball(fname,N):
    radius = 0.3
    zshift = -0.5

    x0 = 0.0
    y0 = 0.0
    z0 = zshift + 0.0
    img = 0.0002*np.ones((N,N,N))
    _x = np.arange(N)
    _y = np.arange(N)
    _z = np.arange(N)
    [Z,Y,X] = np.meshgrid(_z,_y,_x) #N*N*N grid

    dist_from_center = np.sqrt((X-x0)**2 + (Y-y0)**2 +(Z-z0)**2)

    mask = dist_from_center <=radius
    img[mask] = 0.06

    xc = 0.0
    yc = 0.0
    zc = N/2
    Del_xy = 1.0
    Del_z = 1.0

    res=[Del_xy, Del_xy, Del_z]
    loc=[ xc-(N-1)*Del_xy/2,  yc-(N-1)*Del_xy/2, zc-(N-1)*Del_z/2 ]
    img=1000*img/(.02-.0000226)
    write(img,fname,res,loc)
    return img

if __name__=='__main__':
    fname = './data/phantom_ball.npy'  #numpy binary file
    N = 10  #NxNxN phantom
    create_ball(fname,N)
