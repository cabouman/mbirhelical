import numpy as np
import python.helical_io as hio


# File to create a phantom with a single ball or cylinder, then write out relevant
# information for projection.

def create_cone(r, N):
    """
    Create an NxNxN array of floats, one constant value in a cylinder of given radius
    Args:
        r: radius of ball in pixels
        N: Side length of background cube, in pixels

    Returns:
        numpy array of floats of size NxNxN
        dictionary of image info
    """

    x0 = int(N / 2)
    y0 = int(N / 2)
    # z0 = int(N/2) #zshift + 0.0
    Nx = N
    Ny = N + 1
    Nz = N + 2
    img = 0.0 * np.ones((Nx, Ny, Nz), dtype=np.float32)
    _x = np.arange(img.shape[0])
    _y = np.arange(img.shape[1])
    [X, Y] = np.meshgrid(_x, _y, indexing='ij')  # N*N*N grid

    dist_from_center = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)

    for j in np.arange(Nz):
        radius = r * (j + 1) / Nz
        mask = dist_from_center <= radius
        img[mask, j] = 1.0

    xc = 0.0
    yc = 0.0
    zc = N / 2
    Del_xy = 1.0
    Del_z = 1.0

    res = [Del_xy, Del_xy, Del_z]
    loc = [xc - (N - 1) * Del_xy / 2, yc - (N - 1) * Del_xy / 2, zc - (N - 1) * Del_z / 2]
    # img=1000*img/(.02-.0000226)
    image_info = hio.create_image_info(img, res, loc)
    return img, image_info


def create_ball(r, N):
    """
    Create an NxNxN array of floats, one constant value in a ball of given radius
    Args:
        r:
        N:

    Returns:

    """
    radius = r
    # zshift = -0.5

    x0 = int(N / 2)
    y0 = int(N / 2)
    z0 = int(N / 2)  # zshift + 0.0
    img = 0.0002 * np.ones((N, N, N))
    _x = np.arange(N)
    _y = np.arange(N)
    _z = np.arange(N)
    [X, Y, Z] = np.meshgrid(_x, _y, _z)  # N*N*N grid

    dist_from_center = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2 + (Z - z0) ** 2)

    mask = dist_from_center <= radius
    img[mask] = 0.06

    xc = 0.0
    yc = 0.0
    zc = N / 2
    Del_xy = 1.0
    Del_z = 1.0

    res = [Del_xy, Del_xy, Del_z]
    loc = [xc - (N - 1) * Del_xy / 2, yc - (N - 1) * Del_xy / 2, zc - (N - 1) * Del_z / 2]
    img = 1000 * img / (.02 - .0000226)

    image_info = hio.create_image_info(img, res, loc)
    return img, image_info


if __name__ == '__main__':
    base_name = 'cylinder_phantom'
    fname = base_name + '_info.txt'
    N = 100
    r = 10
    # create_ball(fname,N)
    image, img_info = create_cone(r, N)

    hio.write_image_info(fname, img_info, base_name)
    hio.write_image(base_name, img_info, image)
    a = 0
