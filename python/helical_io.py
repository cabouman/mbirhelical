import numpy as np
# import vtk
# from vtk.util import numpy_support

"""
File to mimic the io operations in mbirhelical io.c 
"""


def create_image_info(image, res, loc):
    """
    Args:
        image:
        res:
        loc:

    Returns:

    """
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
        image_file_name: filename of image to be included in image_info (should have no extension)

    Returns:
        None
    """
    img_info['imgFile'] = image_file_name + '.mhimage'
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
    """
    Save an image in two formats:
        1. mhimage, which has 3 binary ints given the x, y, z dims, followed by data in 32-bit float form
        2. raw, which has just the data in 32 bit form.  (This is mostly for debugging).
        TODO: remove the raw format
    Args:
        filename: base name of file at which to store image - should have no extension
        img_info: an image_info dict as returned by create_image_info
        image: a 3D numpy array containing the image data

    Returns:
        None
    """
    with open(filename+'.mhimage', 'w+b') as fp:
        image_dims = np.array([img_info[s] for s in ['Nx', 'Ny', 'Nz']], dtype=int)

        # TODO: check for exact data type size and endianness to match c code.
        fp.write(image_dims.astype(np.int32).tobytes())
        fp.write(image.astype(np.float32).tobytes())

    with open(filename+'.raw', 'w+b') as fp:
        fp.write(image.transpose().astype(np.float32).tobytes())


def write_angle_list(filename, angle_list):
    """
    Save an angle list to a text file, one per line:
    Args:
        filename: file name for angle list
        angle_list: a 1D numpy array containing the angles

    Returns:
        None
    """
    with open(filename, 'w') as fp:
        for angle in angle_list:
            print("%f" % angle, file=fp)

# def numpy_array_as_vtk_image_data(source_numpy_array):
#     """
#     :param source_numpy_array: source array with 2-3 dimensions. If used, the third dimension represents the channel count.
#     Note: Channels are flipped, i.e. source is assumed to be BGR instead of RGB (which works if you're using cv2.imread function to read three-channel images)
#     Note: Assumes array value at [0,0] represents the upper-left pixel.
#     :type source_numpy_array: np.ndarray
#     :return: vtk-compatible image, if conversion is successful. Raises exception otherwise
#     :rtype vtk.vtkImageData
#     """
#
#     if len(source_numpy_array.shape) > 2:
#         channel_count = source_numpy_array.shape[2]
#     else:
#         channel_count = 1
#
#     output_vtk_image = vtk.vtkImageData()
#     output_vtk_image.SetDimensions(source_numpy_array.shape[1], source_numpy_array.shape[0], channel_count)
#
#     vtk_type_by_numpy_type = {
#         np.uint8: vtk.VTK_UNSIGNED_CHAR,
#         np.uint16: vtk.VTK_UNSIGNED_SHORT,
#         np.uint32: vtk.VTK_UNSIGNED_INT,
#         np.uint64: vtk.VTK_UNSIGNED_LONG if vtk.VTK_SIZEOF_LONG == 64 else vtk.VTK_UNSIGNED_LONG_LONG,
#         np.int8: vtk.VTK_CHAR,
#         np.int16: vtk.VTK_SHORT,
#         np.int32: vtk.VTK_INT,
#         np.int64: vtk.VTK_LONG if vtk.VTK_SIZEOF_LONG == 64 else vtk.VTK_LONG_LONG,
#         np.float32: vtk.VTK_FLOAT,
#         np.float64: vtk.VTK_DOUBLE
#     }
#     vtk_datatype = vtk_type_by_numpy_type[source_numpy_array.dtype.type]
#
#     source_numpy_array = np.flipud(source_numpy_array)
#
#     # Note: don't flip (take out next two lines) if input is RGB.
#     # Likewise, BGRA->RGBA would require a different reordering here.
#     if channel_count > 1:
#         source_numpy_array = np.flip(source_numpy_array, 2)
#
#     depth_array = numpy_support.numpy_to_vtk(source_numpy_array.ravel(), deep=True, array_type = vtk_datatype)
#     depth_array.SetNumberOfComponents(channel_count)
#     output_vtk_image.SetSpacing([1, 1, 1])
#     output_vtk_image.SetOrigin([-1, -1, -1])
#     output_vtk_image.GetPointData().SetScalars(depth_array)
#
#     output_vtk_image.Modified()
#    return output_vtk_image
