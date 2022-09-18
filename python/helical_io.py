import numpy as np
"""
File to mimic the io operations in mbirhelical io.c 
"""


def create_image_info(image, res, loc):

    image_info = dict()
    image_info['Nx'] = image.shape[0]
    image_info['Ny'] = image.shape[1]
    image_info['Nz'] = image.shape[2]
    image_info['Nz_mid'] = image.shape[2]
    image_info['xc'] = loc[0]
    image_info['yc'] = loc[1]
    image_info['zc'] = loc[2]
    image_info['Del_xy'] = res[0]
    image_info['Del_z'] = res[2]
    image_info['rI'] = np.amin(image.shape) // 2
    image_info['imgFile'] = 'NA'
    image_info['maskFile'] = 'NA'

    return image_info


def write_image_info(filename, img_info, image_file_name):
    """
    Args:
        filename: file name to save info
        img_info: dict containing relevant image info

    Returns:
        None
    """
    img_info['imgFile'] = image_file_name
    with open(filename, 'w') as fp:
        print("number of voxels in x", file=fp)
        print("%d\n" % img_info['Nx'], file=fp)
        print("number of voxels in y", file=fp)
        print("%d\n" % img_info['Ny'], file=fp)
        print("number of voxels in z (good slices)", file=fp)
        print("%d\n" % img_info['Nz_mid'], file=fp)
        print("number of voxels in z (total)", file=fp)
        print("%d\n" % img_info['Nz'], file=fp)
        print("x coordinate of the center voxel (mm)", file=fp)
        print("%f\n" % img_info['xc'], file=fp)
        print("y coordinate of the center voxel (mm)", file=fp)
        print("%f\n" % img_info['yc'], file=fp)
        print("z coordinate of the center voxel (mm)", file=fp)
        print("%f\n" % img_info['zc'], file=fp)
        print("voxel spacing in xy (mm)", file=fp)
        print("%f\n" % img_info['Del_xy'], file=fp)
        print("voxel spacing in z (mm)", file=fp)
        print("%f\n" % img_info['Del_z'], file=fp)
        print("radius of the reconstruction mask (mm)", file=fp)
        print("%f\n" % img_info['rI'], file=fp)
        print("initial reconstruction image location", file=fp)
        print("%s\n" % img_info['imgFile'], file=fp)
        print("mask file location", file=fp)
        print("%s\n" % img_info['maskFile'], file=fp)


def write_image(filename, img_info, image):
    with open(filename, 'w+b') as fp:
        image_dims = np.array([img_info[s] for s in ['Nx', 'Ny', 'Nz']], dtype=int)

        # TODO: check for exact data type size and endianness to match c code.
        fp.write(image_dims.astype(np.int32).tobytes())
        fp.write(image.astype(np.float64).tobytes())
