from unittest import TestCase

from raycorr.processor import main
import numpy as np
from snappy import Product
# from snappy import ProductData
from snappy import ProductIO
# from snappy import ProductUtils
# from snappy import FlagCoding
# from snappy import jpy


class RayCorrTest(TestCase):
    def test_end_to_end(self):
        PRODPATH = "C:\\Users\\carsten\\Dropbox\\Carsten\\SWProjects\\Rayleigh-Correction\\testdata\\"
        # PRODPATH = "D:\\Dropbox\\Carsten\\SWProjects\\Rayleigh-Correction\\testdata\\"

        # validate here
        numerr=0
        # read output
        # SENSOR = 'MERIS'
        SENSOR = 'OLCI'
        main([SENSOR])
        if (SENSOR=='MERIS'):
            REF_FILE =OUT_FILE = PRODPATH+'Reftestprodukt1_MER_RR_20050713.dim'
            TEST_FILE=OUT_FILE = PRODPATH+'Testprodukt1_MER_RR_20050713.dim'
            NBANDS=15
        if (SENSOR=='OLCI'):
            REF_FILE = PRODPATH+'Reftestproduct3_S3A_OL_1_EFR____20160509T103945.dim'
            TEST_FILE = PRODPATH+'Testproduct3_OL_1_EFR____20160509T103945.dim'
            NBANDS=21

        print("Opening reference product ...")
        refproduct = ProductIO.readProduct(REF_FILE)
        width = refproduct.getSceneRasterWidth()
        height = refproduct.getSceneRasterHeight()

        print("Opening test product ...")
        testproduct = ProductIO.readProduct(TEST_FILE)
        widthtest = testproduct.getSceneRasterWidth()
        heighttest = testproduct.getSceneRasterHeight()

        # compare
        print("Start comparing ...")
        try:
            self.assertEqual(width, widthtest)
            print("  widths agree: ",width)
        except:
            print("  widths error: ref=",width," test=",widthtest)
            numerr+=1

        try:
            self.assertEqual(height, heighttest)
            print("  height agree: ",height)
        except:
            print("  height error: ref=",height," test=",heighttest)
            numerr+=1

        refvalues=np.zeros((width,height),dtype=np.float32)
        testvalues=np.zeros((width,height),dtype=np.float32)
        bandname="rBRR_"
        for i in range(1,NBANDS+1):
            refsource = refproduct.getBand(bandname+str(i))
            refvalues = refsource.readPixels(0, 0, width, height, refvalues)
            testsource = testproduct.getBand(bandname+str(i))
            testvalues = testsource.readPixels(0, 0, width, height, testvalues)
            try:
                result=np.allclose(refvalues,testvalues, rtol=0.0, atol=1e-6, equal_nan=True)
                self.assertEqual(result,True)
                print("   ",bandname+str(i)," agrees")
            except:
                print("   ",bandname+str(i)," disagrees")
                numerr+=1
        print("toal number of tests failed =", numerr)