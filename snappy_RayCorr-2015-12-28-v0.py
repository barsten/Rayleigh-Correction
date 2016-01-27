import sys

import math
import numpy
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
from snappy import FlagCoding

if len(sys.argv) != 2:
    print("usage: %s <file>" % sys.argv[0])
    sys.exit(1)

file = sys.argv[1]

print("Reading...")
product = ProductIO.readProduct(file)
width = product.getSceneRasterWidth()
height = product.getSceneRasterHeight()
name = product.getName()
description = product.getDescription()
band_names = product.getBandNames()

print("Product:     %s, %s" % (name, description))
print("Raster size: %d x %d pixels" % (width, height))
print("Start time:  " + str(product.getStartTime()))
print("End time:    " + str(product.getEndTime()))
print("Bands:       %s" % (list(band_names)))

raycorProduct = Product('RayCorr', 'RayCorr', width, height)
writer = ProductIO.getProductWriter('BEAM-DIMAP')
raycorProduct.setProductWriter(writer)

# Create TOA reflectance bands
nbands = product.getNumBands()-2 # the last 2 bands are l1flags and detector index; we don't need them
for i in range(nbands):
    bsource = product.getBandAt(i)
    b_name = "rtoa_"+str(i+1)
    toareflBand = raycorProduct.addBand(b_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,toareflBand)

# create band for Rayleigh optical thickness
taurBand = raycorProduct.addBand('taur', ProductData.TYPE_FLOAT32)
raycorProduct.setAutoGrouping('rtoa')

# Create flag coding
raycorFlagsBand = raycorProduct.addBand('raycor_flags', ProductData.TYPE_UINT8)
raycorFlagCoding = FlagCoding('raycor_flags')
raycorFlagCoding.addFlag("testflag_1", 1, "Flag 1 for Rayleigh Correction, reflectance 1 > 0.3")
raycorFlagCoding.addFlag("testflag_2", 2, "Flag 2 for Rayleigh Correction, reflectance 13 < 0.05")
group = raycorProduct.getFlagCodingGroup()
group.add(raycorFlagCoding)
raycorFlagsBand.setSampleCoding(raycorFlagCoding)

# add geocoding and create the product on disk (meta data, empty bands)
# question to developers: why do I have to do this all at the beginning. I tried putting these two lines to the end of the code and the product was not created correctly
# ProductUtils.copyGeoCoding(product, raycorProduct) #geocoding is copied when tie point grids are copied,
# i.e. the next copy makes this one redundant (actually it leads to an error beciase lat,lon would be copied twice
ProductUtils.copyTiePointGrids(product, raycorProduct)
raycorProduct.writeHeader('raycor_output.dim')

# calculate and write toa reflectances
# ====================================
# first, we need some arrays and other variables to store stuff that we need

# Dummy code to create a band
b1 = product.getBand('radiance_1')
b2 = product.getBand('radiance_2')

radiance = numpy.zeros(width, dtype=numpy.float32)
r1 = numpy.zeros(width, dtype=numpy.float32)
r2 = numpy.zeros(width, dtype=numpy.float32)
taur = numpy.zeros(width, dtype=numpy.float32)

print("Processing ...")

for y in range(height):
    # start radiance to reflectance conversion
    for i in range(nbands):
        b_source = product.getBandAt(i)
        radiance = b_source.readPixels(0, y, width, 1, radiance)
        E0 = b_source.getSolarFlux()
        reflectance = radiance * math.pi / E0
        b_out = raycorProduct.getBandAt(i)
        b_out.writePixels(0, y, width, 1, reflectance)
        # this is dummy code to create a flag
        if i == 0:
            flag1 = reflectance > 0.3
        if i == 12:
            flag2 = reflectance < 0.05
    raycorFlags = flag1 + 2*flag2
    raycorFlagsBand.writePixels(0, y, width, 1, raycorFlags)
    # end flags dummy code
    # radiance to reflectance conversion completed

    # start Rayleigh optical thickness calculation
    print("processing line ", y, " of ", height)
    r1 = b1.readPixels(0, y, width, 1, r1)
    r2 = b2.readPixels(0, y, width, 1, r2)

    for x in range(width):
        z1 = y+x
        taur[x]=z1

    taurBand.writePixels(0, y, width, 1, taur)
    # Rayleigh optical thickness calculation completed

raycorProduct.closeIO()

print("Done.")