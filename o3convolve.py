import numpy as np
import struct

class O3:

    def __init__(self, file_path):
        self.file_path = file_path
        self.coeffhighres = _read_O3_coeff(file_path)

    def convolve(self, lower, upper):
        O3absorption = self.coeffhighres[:,1]
        O3wavelength = self.coeffhighres[:,0]
        numval = O3wavelength.__len__()
        weight = np.zeros(numval)
        for i in range(numval):
            if (O3wavelength[i]>= lower and O3wavelength[i]<= upper):
                weight[i]=1
        O3value = np.average(O3absorption, weights=weight)
        return O3value

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
    # print(m)

    lower = 550 - 5.0
    upper = 550 + 5.0
    O3b1 = ozone.convolve(lower, upper)
    print(O3b1)
