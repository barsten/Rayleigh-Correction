import numpy as np


class O3:
    def __init__(self, file_path):
        self.file_path = file_path
        self.coeffhighres = _read_O3_coeff(file_path)

    def convolve(self, lower, upper):
        O3absorption = self.coeffhighres[:, 1]
        O3wavelength = self.coeffhighres[:, 0]
        numval = O3wavelength.__len__()
        weight = np.zeros(numval)
        for i in range(numval):
            if (O3wavelength[i] >= lower and O3wavelength[i] <= upper):
                weight[i] = 1
        O3value = np.average(O3absorption, weights=weight)
        return O3value

    def convolveInstrument(self, instrument):
        """
        :param instrument: name of sensor, so that the right bands can be chosen
        :return: numpy array with ozone absorption coefficients for each wavelength of instrumenet, in Dobson units
        """
        o3absorpInstrument = 0.0
        if (instrument == 'MERIS'):
            absorb_ozon = np.array([0.0002174, 0.0034448, 0.0205669, 0.0400134, 0.105446, 0.1081787, 0.0501634,
                                    0.0349671, 0.0187495, 0.0086322, 0.0000001, 0.0084989, 0.0018944, 0.0012369,
                                    0.000001])  # MERIS
            lamC = np.array(
                    [412.5, 442.0, 490.0, 510.0, 560.0, 620.0, 665.0, 681.25, 708.75, 753.0, 761.25, 779.0, 865.0,
                     885.0, 900])
            lamW = np.array([10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 3.75, 10., 20., 20., 40.])
            o3absorpInstrument = np.zeros(15, dtype=np.float32)
            for i in range(15):
                lower = lamC[i] - lamW[i] / 2
                upper = lamC[i] + lamW[i] / 2
                o3absorpInstrument[i] = O3.convolve(self, lower, upper)
                # print(i, absorb_ozon[i], o3absorpInstrument[i], 100*(absorb_ozon[i]-o3absorpInstrument[i])/absorb_ozon[i])

        if (instrument == 'OLCI'):
            lamC = np.array(
                    [400.0, 412.5, 442.0, 490.0, 510.0, 560.0, 620.0, 665.0, 673.75, 681.25, 708.75, 753.75, 761.25, 764.375, 767.5, 778.75, 865.0,
                     885.0, 900.0, 940.0, 1020.0])
            lamW = np.array([15., 10., 10., 10., 10., 10., 10., 10., 7.5, 7.5, 10., 7.5, 2.5, 3.75, 2.5, 15., 20., 10., 10., 20., 40.])
            o3absorpInstrument = np.zeros(21, dtype=np.float32)
            for i in range(21):
                lower = lamC[i] - lamW[i] / 2
                upper = lamC[i] + lamW[i] / 2
                o3absorpInstrument[i] = O3.convolve(self, lower, upper)

        return o3absorpInstrument


def _read_O3_coeff(file_path):
    """
    :param file_path: the file path
    :return: Ozone absorption cooefficients at 1nm resolution, in Dobson Units

    """

    coeff = np.loadtxt(file_path, skiprows=2)
    return coeff


if __name__ == '__main__':
    import os.path
    AUX_FILE = os.path.join(os.path.dirname(__file__), 'ozone-highres.txt')

    ozone = O3(AUX_FILE)
    m = ozone.coeffhighres
    # print(m)

    lower = 550 - 5.0
    upper = 550 + 5.0
    O3b1 = ozone.convolve(lower, upper)
    print(O3b1)

    print(ozone.convolveInstrument('MERIS'))
