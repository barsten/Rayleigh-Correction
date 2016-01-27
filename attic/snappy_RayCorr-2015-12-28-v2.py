# coding=utf-8

# Rayleigh Correction Processor
# MERIS only
# best version on the evening of 28.12.2015
# creates TOA reflectances
# calculates Rayleigh Optical Thickness in all bands
# code is clean

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

nbands = product.getNumBands()-2 # the last 2 bands are l1flags and detector index; we don't need them
# Create TOA reflectance and Rayleig optical thickness bands
for i in range(nbands):
    bsource = product.getBandAt(i)
    btoa_name = "rtoa_"+str(i+1)
    toareflBand = raycorProduct.addBand(btoa_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,toareflBand)
    btaur_name = "taur_"+str(i+1)
    taurBand = raycorProduct.addBand(btaur_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,taurBand)

raycorProduct.setAutoGrouping('rtoa:taur')

# Create flag coding
raycorFlagsBand = raycorProduct.addBand('raycor_flags', ProductData.TYPE_UINT8)
raycorFlagCoding = FlagCoding('raycor_flags')
raycorFlagCoding.addFlag("testflag_1", 1, "Flag 1 for Rayleigh Correction")
raycorFlagCoding.addFlag("testflag_2", 2, "Flag 2 for Rayleigh Correction")
group = raycorProduct.getFlagCodingGroup()
group.add(raycorFlagCoding)
raycorFlagsBand.setSampleCoding(raycorFlagCoding)

# add geocoding and create the product on disk (meta data, empty bands)
# question to developers: why do I have to do this all at the beginning. I tried putting these two lines to the end of the code and the product was not created correctly
# ProductUtils.copyGeoCoding(product, raycorProduct) #geocoding is copied when tie point grids are copied,
# i.e. the next copy makes this one redundant (actually it leads to an error beciase lat,lon would be copied twice
ProductUtils.copyTiePointGrids(product, raycorProduct)
raycorProduct.writeHeader('raycor_output.dim')

# Calculate and write toa reflectances and Rayleigh optical thickness
# ===================================================================
# constants
AVO = 6.0221367E+23 # Avogadro's number
m_a_zero = 28.9595 # Mean molecular weight of dry ait (zero CO2)
g0_45 = 980.616 # Acceleration of gravity (sea level and 458 latitude)
Ns = 2.5469E19 # Molecular density of gas in molecules / cm3

# constants describing the state of the atmosphere and which we don't know; better values may be used if known
CO2 = 3E-4 # CO2 concentration at pixel; typical values are 300 to 360 ppm
C_CO2 = CO2*100 # CO2 concentration in ppm
m_a = 15.0556 * CO2 + m_a_zero # mean molecular weight of dry air as function of actual CO2

# arrays which are needed to store some stuff
radiance = numpy.zeros(width, dtype=numpy.float32)
taur = numpy.zeros((nbands, width), dtype=numpy.float32)
sigma = numpy.zeros(nbands, dtype=numpy.float32)

tp_alt = product.getTiePointGrid('dem_alt')
alt = numpy.zeros(width, dtype=numpy.float32)

tp_press = product.getTiePointGrid('atm_press')
press0 = numpy.zeros(width, dtype=numpy.float32)

tp_latitude = product.getTiePointGrid('latitude')
phi = numpy.zeros(width, dtype=numpy.float32)

print("Processing ...")
# Calculate the Rayleigh cross section, which depends only on wavelength but not on air pressure
for i in range(nbands):
    print("processing Rayleigh cross section of band", i)
    b_source = product.getBandAt(i)
    lam = b_source.getSpectralWavelength() # wavelength of MERIS band i in nm
    lam = lam / 1000.0 # wavelength in micrometer
    lam2 = lam / 10000 # wavelength in cm
    F_N2 = 1.034+0.000317/(lam**2) # King factor of N2
    F_O2 = 1.096+0.001385/(lam**2)+0.0001448/(lam**4) # King factor of O2
    F_air = (78.084*F_N2+20.946*F_O2+0.934*1+C_CO2*1.15)/(78.084+20.946+0.934+C_CO2) #depolarization ratio or King Factor, (6+3rho)/(6-7rho)
    n_ratio = 1+0.54*(CO2-0.0003)
    n_1_300 = (8060.51+(2480990/(132.274-lam**(-2)))+(17455.7/(39.32957-lam**(-2))))/100000000
    nCO2 = n_ratio*(1+n_1_300) # reflective index at CO2
    sigma[i] = (24*math.pi**3*(nCO2**2-1)**2)/(lam2**4*Ns**2*(nCO2**2+2)**2)*F_air

for y in range(height):
    print("processing line ", y, " of ", height)
    # start radiance to reflectance conversion
    for i in range(nbands):
        b_source = product.getBand("radiance_"+str(i+1))
        radiance = b_source.readPixels(0, y, width, 1, radiance)
        E0 = b_source.getSolarFlux()
        reflectance = radiance * math.pi / E0
        # btoa_name = "rtoa_"+str(i+1)
        b_out = raycorProduct.getBand("rtoa_"+str(i+1))
        b_out.writePixels(0, y, width, 1, reflectance)
    # radiance to reflectance conversion completed

    # this is dummy code to create a flag
    flag1 = numpy.zeros(width, dtype=numpy.bool_)
    flag2 = numpy.zeros(width, dtype=numpy.bool_)
    raycorFlags = flag1 + 2*flag2
    raycorFlagsBand.writePixels(0, y, width, 1, raycorFlags)
    # end flags dummy code

    # start Rayleigh optical thickness calculation
    alt = tp_alt.readPixels(0, y, width, 1, alt)
    press0 = tp_press.readPixels(0, y, width, 1, press0)
    phi = tp_latitude.readPixels(0, y, width, 1, phi)

    # Now calculate the pressure depending terms and finally the Rayleigh optical thickness
    for x in range(width):
        # Calculation to get the pressure
        z = alt[x] # altitude at pixel in meters, taken from MERIS tie-point grid
        z = max(z, 0) # clip to sea level
        Psurf0 = press0[x] # pressure at sea level in hPa, taken from MERIS tie-point grid
        Psurf = Psurf0 * (1 - 0.0065*z/288.15)**5.255 # air pressure at the pixel (i.e. at altitude) in hPa, using the international pressure equation
        P = Psurf * 1000 # air pressure at pixel location in dyn / cm2, which is hPa * 1000
        # calculation to get the constant of gravity at the pixel altitude, taking the air mass above into account
        dphi = math.radians(phi[x]) # latitude in radians
        cos2phi = math.cos(2*dphi)
        g0 = g0_45*(1-0.0026373*cos2phi+0.0000059*cos2phi**2)
        zs = 0.73737*z + 5517.56 # effective mass-weighted altitude
        g = g0-(0.0003085462 + 0.000000227*cos2phi)*zs + (0.00000000007254 + 0.0000000000001 * cos2phi)*zs**2 - (1.517E-17 + 6E-20*cos2phi)*zs**3
        # calculations to get the Rayeigh optical thickness
        factor = (P * AVO) / (m_a * g)
        for i in range(nbands):
            taur[(i,x)]=sigma[i] * factor

    for i in range(nbands):
        btaur_name = "taur_"+str(i+1)
        taurBand = raycorProduct.getBand(btaur_name)
        taurBand.writePixels(0, y, width, 1, taur[i])
    # Rayleigh optical thickness calculation completed

raycorProduct.closeIO()

print("Done.")