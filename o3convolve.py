import numpy as np
import struct

class O3:

    def __init__(self, file_path):
        self.file_path = file_path
        self.coeffhighres = _read_O3_coeff(file_path)

    def convolve(self, lower, upper):
        return 123.4

def _read_O3_coeff(file_path):
    """
    :param file_path: the file path
    :return: Ozone absorption cooefficients at 1nm resolution, in Dobson Units

    :param file_path:
    :return:
    """
    # with open(file_path, 'r') as fp:
    #     fp.readline()
    #     fp.readline()
    coeff = np.loadtxt(file_path,skiprows=2)
    return coeff

if __name__ == '__main__':
    AUX_FILE = 'C:\\Users\\carsten\\Dropbox\\Carsten\\Tagesordner\\20160104\\Rayleigh-Correction-Processor\\ozone-highres.txt'

    ozone = O3(AUX_FILE)
    m = ozone.coeffhighres
    print(m)

    lower = 412.5 - 5.0
    upper = 412.5 + 5.0
    O3b1 = ozone.convolve(lower, upper)
    print(O3b1)
