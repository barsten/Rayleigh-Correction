import numpy as np
import struct

class ADF:

    def __init__(self, file_path):
        self.file_path = file_path
        self.ray_coeff_matrix = _read_ray_coeff_matrix(file_path)


def _read_ray_coeff_matrix(file_path):
    """
    Atmosphere - MPH 1247
    Atmosphere - SPH 2898
    Atmosphere - GADS General 2772
    Atmosphere - GADS Optical Thicknesses 300
    Atmosphere - GADS H2O Transmission 912
    Atmosphere - GADS Rayleigh Scattering Function 3744
    Atmosphere - GADS Rayleigh Spherical Albedo 68
    Atmosphere - ADS O2 transmission around 779 nm 5985000
    Atmosphere - ADS Apparent Pressure Parameters 1161216
    Atmosphere - GADS spare 3120
    Atmosphere - ADS Rayleigh Reflectance over Ocean 1345500
    Atmosphere - GADS Photosynthetically Available Radiation 640000

    :param file_path: the file path
    :return: rayleigh transmission cooefficients
    """
    _OFFSET = 1247 + 2898 + 2772 + 300 + 912
    _NUM_THETA = 12
    _ORDER = 3
    _NUM_COEFFS = 4
    with open(file_path, 'rb') as fp:
        fp.seek(_OFFSET)
        matrix = np.zeros(shape=(_ORDER, _NUM_THETA, _NUM_THETA, _NUM_COEFFS), dtype=np.float32)
        for i_order in range(_ORDER):
            for i_theta_s in range(_NUM_THETA):
                for i_theta_v in range(i_theta_s, _NUM_THETA):
                    for i_coeff in range(_NUM_COEFFS):
                        data = fp.read(4)
                        elem = struct.unpack('>1f', data)
                        matrix[i_order, i_theta_s, i_theta_v, i_coeff] = elem[0]
                        if i_theta_v != i_theta_s:
                            matrix[i_order, i_theta_v, i_theta_s, i_coeff] = matrix[
                                i_order, i_theta_s, i_theta_v, i_coeff]
    return matrix


if __name__ == '__main__':
    AUX_FILE = 'C:\\Users\\carsten\\Dropbox\\Carsten\\Tagesordner\\20160104\\Rayleigh-Correction-Processor\\' \
               'ADF\\MER_ATP_AXVACR20091026_144725_20021224_121445_20200101_000000'

    adf = ADF(AUX_FILE)
    m = adf.ray_coeff_matrix
    print(m)
