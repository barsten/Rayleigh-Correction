import sys

import math
from collections import OrderedDict

import numpy as np
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
from snappy import FlagCoding
from snappy import jpy
from raycorr.readfile import readRayADF
from raycorr.read_RayCoeff import ADF
from raycorr.o3convolve import O3

import json

class JSONNumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def json_as_numpy(dct):
    return {key: np.array(value) for key, value in dct.items()}

def main():
    AUX_FILE = 'C:\\Users\\carsten\\Dropbox\\Carsten\\Tagesordner\\20160104\\Rayleigh-Correction-Processor\\' \
           'ADF\\MER_ATP_AXVACR20091026_144725_20021224_121445_20200101_000000'

    adf = ADF(AUX_FILE)
    ray_coeff_matrix = adf.ray_coeff_matrix
    rayADF = readRayADF(AUX_FILE)

    new_aux = OrderedDict()
    new_aux['tau_ray'] = rayADF['tR']
    new_aux['theta'] = rayADF['theta']
    new_aux['ray_albedo_lut'] = rayADF['rayAlbLUT']
    new_aux['ray_coeff_matrix'] = ray_coeff_matrix

    json_str = json.dumps(new_aux, cls=JSONNumpyEncoder, indent=2)
    # print(json_str)
    with open('raycorr_auxdata.json', 'w') as fp:
        fp.write(json_str)
        fp.close()
    # verify that file can be read correctly
    obj1 = json.loads(json_str, object_hook=json_as_numpy)
    with open('raycorr_auxdata.json', 'r') as fp:
        obj2 = json.load(fp, object_hook=json_as_numpy)
        fp.close()
    # print(obj1)
    # print(obj2)
    print(cmp(obj1, obj2))
    print('hello world')

if __name__ == '__main__':
    main()
