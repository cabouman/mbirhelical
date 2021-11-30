# @Author: Tanushree Gupta<tshreegupta@gmail.com>
# @Date:   2021-11-29T19:29:50-05:00
# @Last modified by:   Tanushree Gupta<tshreegupta@gmail.com>
# @Last modified time: 2021-11-29T19:34:28-05:00

import numpy as np

## Write 3D voxel to the file in below format
# ASCII header will contain 3 lines
#   Nx Ny Nz
#   Dx Dy Dz  : voxel size (mm) (coming from res=[Dx Dy Dz])
#   x0 y0 z0  : location of 1st index (coming from loc=[x0 y0 z0])
#
# The rest is little endian 16-bit unsigned short representing Hounsfield
def write(img,fname,res,loc):
    with open(fname,'wb') as fp:
        np.save(fp,img.shape)
        np.save(fp,res)
        np.save(fp,loc)
        np.save(fp,img)

# ASCII header should contain 3 lines
#   Nx Ny Nz
#   Dx Dy Dz  : voxel size (mm) (coming from res=[Dx Dy Dz])
#   x0 y0 z0  : location of 1st index (coming from loc=[x0 y0 z0])
#
# The rest is little endian 16-bit unsigned short representing Hounsfield
#
# Input
#   fname : filename
# Outputs
#   img : 3d image array in order (z,y,x)
#   res : voxel size [Dx Dy Dz]
#   loc : coords of 1st voxel [x0 y0 z0]
def read(fname):
    with open(fname, 'rb') as fp:
        [Nz,Ny,Nx] = np.load(fp)
        res = np.load(fp)
        loc = np.load(fp)
        img = np.load(fp)

    return img, res, loc
