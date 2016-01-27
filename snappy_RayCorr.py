# coding=utf-8

# Rayleigh Correction Processor
# MERIS only
# experimental version on the evening of 28.12.2015
# creates TOA reflectances
# calculates Rayleigh Optical Thickness in all bands
# calculates the Rayleigh reflectances

import sys

import math
import numpy as np
from scipy.interpolate import interpn
from scipy.interpolate import interp1d
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
from snappy import FlagCoding
from snappy import jpy
from readfile import readRayADF
from read_RayCoeff import ADF

DEMFactory = jpy.get_type('org.esa.snap.dem.dataio.DEMFactory')
Resampling = jpy.get_type('org.esa.snap.core.dataop.resamp.Resampling')
GeoPos = jpy.get_type('org.esa.snap.core.datamodel.GeoPos')

AUX_FILE = 'C:\\Users\\carsten\\Dropbox\\Carsten\\Tagesordner\\20150814\\Rayleigh-Correction-Processor\\' \
               'ADF\\MER_ATP_AXVACR20091026_144725_20021224_121445_20200101_000000'

adf = ADF(AUX_FILE)
ray_coeff_matrix = adf.ray_coeff_matrix
rayADF = readRayADF()

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

    brhor_name = "rRay_"+str(i+1)
    rhorBand = raycorProduct.addBand(brhor_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rhorBand)
    # Fourier Terms, during debugging only
    brhorF1_name = "rRayF1_"+str(i+1)
    rhorF1Band = raycorProduct.addBand(brhorF1_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rhorF1Band)
    brhorF2_name = "rRayF2_"+str(i+1)
    rhorF2Band = raycorProduct.addBand(brhorF2_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rhorF2Band)
    brhorF3_name = "rRayF3_"+str(i+1)
    rhorF3Band = raycorProduct.addBand(brhorF3_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rhorF3Band)
    rayTransS_name = "transSRay_"+str(i+1)
    rayTransSBand = raycorProduct.addBand(rayTransS_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rayTransSBand)
    rayTransV_name = "transVRay_"+str(i+1)
    rayTransVBand = raycorProduct.addBand(rayTransV_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rayTransVBand)
    sARay_name = "sARay_"+str(i+1)
    sARayBand = raycorProduct.addBand(sARay_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,sARayBand)
    rtoaR_name = "rtoaRay_"+str(i+1)
    rtoaRBand = raycorProduct.addBand(rtoaR_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rtoaRBand)
    rBRR_name = "rBRR_"+str(i+1)
    rBRRBand = raycorProduct.addBand(rBRR_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,rBRRBand)
    spf_name = "sphericalAlbedoFactor_"+str(i+1)
    spfBand = raycorProduct.addBand(spf_name, ProductData.TYPE_FLOAT32)
    ProductUtils.copySpectralBandProperties(bsource,spfBand)

raycorProduct.setAutoGrouping('rtoa:taur:rRay:rRayF1:rRayF2:rRayF3:transSRay:transVRay:sARay:rtoaRay:rBRR:sphericalAlbedoFactor')

airmassBand = raycorProduct.addBand('airmass', ProductData.TYPE_FLOAT32)
azidiffBand = raycorProduct.addBand('azidiff', ProductData.TYPE_FLOAT32)
altBand = raycorProduct.addBand('altitude', ProductData.TYPE_FLOAT32)

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
raycorProduct.writeHeader('raycor_output_5.2.dim')

# Calculate and write toa reflectances and Rayleigh optical thickness
# ===================================================================
# some stuff needed to get the altitude from an external DEM; can be omitted if altitude is used from the product
# resamplingMethod = 'NEAREST_NEIGHBOUR'  # Resampling.NEAREST_NEIGHBOUR.getName()
resamplingMethod = Resampling.NEAREST_NEIGHBOUR.getName()
demName = 'GETASSE30'  # alternative 'SRTM 3Sec'
dem = DEMFactory.createElevationModel(demName, resamplingMethod)

# constants
AVO = 6.0221367E+23 # Avogadro's number
m_a_zero = 28.9595  # Mean molecular weight of dry ait (zero CO2)
g0_45 = 980.616     # Acceleration of gravity (sea level and 458 latitude)
Ns = 2.5469E19      # Molecular density of gas in molecules / cm3

# constants describing the state of the atmosphere and which we don't know; better values may be used if known
CO2 = 3.E-4 # CO2 concentration at pixel; typical values are 300 to 360 ppm
C_CO2 = CO2*100 # CO2 concentration in ppm
m_a = 15.0556 * CO2 + m_a_zero # mean molecular weight of dry air as function of actual CO2

# other constants
PA = 0.9587256      # Rayleigh Phase function, molecular asymetry factor 1
PB = 1.-PA           # Rayleigh Phase function, molecular asymetry factor 2
tpoly = rayADF['tR']  # Polynomial coefficients for Rayleigh transmittance

# arrays which are needed to store some stuff
radiance = np.zeros(width, dtype=np.float32)
reflectance = np.zeros((nbands, width), dtype=np.float32)
taur = np.zeros((nbands, width), dtype=np.float32)
sigma = np.zeros(nbands, dtype=np.float32)
airmass = np.zeros(width, dtype=np.float32)
azidiff = np.zeros(width, dtype=np.float32)
PR = np.zeros(3, dtype=np.float32)  # Fourier coefficients of the Rayleigh Phase function
rho_Rf = np.zeros(3, dtype=np.float32)  # Fourier terms of the Rayleigh primary scattering reflectance
rho_Rm = np.zeros((3, nbands, width), dtype=np.float32)  # Fourier terms of the Rayleigh scattering reflectance, corrected for multiple scattering
rho_R = np.zeros((nbands,width), dtype=np.float32)  # first approximation of Rayleigh reflectance
rho_toaR = np.zeros((nbands,width), dtype=np.float32)  # toa reflectance corrected for Rayleigh scattering
rho_BRR = np.zeros((nbands,width), dtype=np.float32)  # top of aerosol reflectance, which is equal to bottom of Rayleigh reflectance
sphericalFactor = np.zeros((nbands,width), dtype=np.float32)  # spherical Albedo Correction Factor (for testing only, can be integrated into the equation later)

tp_alt = product.getTiePointGrid('dem_alt')
alt = np.zeros(width, dtype=np.float32)

tp_press = product.getTiePointGrid('atm_press')
press0 = np.zeros(width, dtype=np.float32)

tp_latitude = product.getTiePointGrid('latitude')
lat = np.zeros(width, dtype=np.float32)
tp_longitude = product.getTiePointGrid('longitude')
lon = np.zeros(width, dtype=np.float32)

tp_theta_s = product.getTiePointGrid('sun_zenith')
theta_s = np.zeros(width, dtype=np.float32)

tp_theta_v = product.getTiePointGrid('view_zenith')
theta_v = np.zeros(width, dtype=np.float32)

tp_azi_s = product.getTiePointGrid('sun_azimuth')
azi_s = np.zeros(width, dtype=np.float32)

tp_azi_v = product.getTiePointGrid('view_azimuth')
azi_v = np.zeros(width, dtype=np.float32)

# Rayleigh multiple scattering
# - Coefficients LUT
dimTheta = 12
dimThetaS = dimThetaV = dimTheta
gridThetaS = rayADF['theta']
gridThetaV = rayADF['theta']
gridGeometry = [gridThetaS, gridThetaV]
# RayScattCoeffA = rayADF['RayScattCoeffA']
# RayScattCoeffB = rayADF['RayScattCoeffB']
# RayScattCoeffC = rayADF['RayScattCoeffC']
# RayScattCoeffD = rayADF['RayScattCoeffD']
RayScattCoeffA = ray_coeff_matrix[:, :, :, 0]
RayScattCoeffB = ray_coeff_matrix[:, :, :, 1]
RayScattCoeffC = ray_coeff_matrix[:, :, :, 2]
RayScattCoeffD = ray_coeff_matrix[:, :, :, 3]
# - Fourier terms
a = np.zeros(3,dtype=np.float32)
b = np.zeros(3,dtype=np.float32)
c = np.zeros(3,dtype=np.float32)
d = np.zeros(3,dtype=np.float32)
rayMultiCorr = np.zeros(3,dtype=np.float32)

# Rayleigh transmittances and spherical albedo
tR_thetaS = np.zeros((nbands,width), dtype=np.float32)  # Rayleigh Transmittance sun - surface
tR_thetaV = np.zeros((nbands,width), dtype=np.float32)  # Rayleigh Transmittance surface - sun
dimTaur = 17
taurTab = np.linspace(0.0, 1.0, num=dimTaur)
rayAlb_f = interp1d(taurTab, rayADF['rayAlbLUT'])
sARay = np.zeros((nbands,width), dtype=np.float32)  # Rayleigh spherical albedo

print("Processing ...")
# Calculate the Rayleigh cross section, which depends only on wavelength but not on air pressure
for i in range(nbands):
    print("processing Rayleigh cross section of band", i)
    b_source = product.getBandAt(i)
    lam = b_source.getSpectralWavelength() # wavelength of MERIS band i in nm
    lam = lam / 1000.0 # wavelength in micrometer
    lam2 = lam / 10000.0 # wavelength in cm
    F_N2 = 1.034+0.000317/(lam**2) # King factor of N2
    F_O2 = 1.096+0.001385/(lam**2)+0.0001448/(lam**4) # King factor of O2
    F_air = (78.084*F_N2+20.946*F_O2+0.934*1+C_CO2*1.15)/(78.084+20.946+0.934+C_CO2) #depolarization ratio or King Factor, (6+3rho)/(6-7rho)
    n_ratio = 1+0.54*(CO2-0.0003)
    n_1_300 = (8060.51+(2480990./(132.274-lam**(-2)))+(17455.7/(39.32957-lam**(-2))))/100000000.0
    nCO2 = n_ratio*(1+n_1_300) # reflective index at CO2
    sigma[i] = (24*math.pi**3*(nCO2**2-1)**2)/(lam2**4*Ns**2*(nCO2**2+2)**2)*F_air

for y in range(height):
# for y in range(120,129):
    print("processing line ", y, " of ", height)
    # start radiance to reflectance conversion
    theta_s = tp_theta_s.readPixels(0, y, width, 1, theta_s)  # sun zenith angle in degree
    for i in range(nbands):
        b_source = product.getBand("radiance_"+str(i+1))
        radiance = b_source.readPixels(0, y, width, 1, radiance)
        E0 = b_source.getSolarFlux()
        reflectance[i] = radiance * math.pi / (E0 * np.cos(np.radians(theta_s)))
        b_out = raycorProduct.getBand("rtoa_"+str(i+1))
        b_out.writePixels(0, y, width, 1, reflectance[i])
    # radiance to reflectance conversion completed

    # this is dummy code to create a flag
    flag1 = np.zeros(width, dtype=np.bool_)
    flag2 = np.zeros(width, dtype=np.bool_)
    raycorFlags = flag1 + 2*flag2
    raycorFlagsBand.writePixels(0, y, width, 1, raycorFlags)
    # end flags dummy code

    lat = tp_latitude.readPixels(0, y, width, 1, lat)
    lon = tp_longitude.readPixels(0, y, width, 1, lon)

    # start Rayleigh optical thickness calculation
    # alt = tp_alt.readPixels(0, y, width, 1, alt)  # using the tie-point DEM in a MERIS product
    # get the altitude from an external DEM
    for x in range(width): alt[x] = dem.getElevation(GeoPos(lat[x], lon[x]))

    press0 = tp_press.readPixels(0, y, width, 1, press0)

    theta_s = tp_theta_s.readPixels(0, y, width, 1, theta_s)  # sun zenith angle in degree
    theta_v = tp_theta_v.readPixels(0, y, width, 1, theta_v)   # view zenith angle in degree
    azi_s = tp_azi_s.readPixels(0, y, width, 1, azi_s)  # sun azimuth angle in degree
    azi_v = tp_azi_v.readPixels(0, y, width, 1, azi_v)  # view azimuth angle in degree

    # Now calculate the pixel dependent terms (like pressure) and finally the Rayleigh optical thickness
    for x in range(width):
        # Calculation to get the pressure
        z = alt[x] # altitude at pixel in meters, taken from MERIS tie-point grid
        z = max(z, 0) # clip to sea level
        Psurf0 = press0[x] # pressure at sea level in hPa, taken from MERIS tie-point grid
        Psurf = Psurf0 * (1. - 0.0065*z/288.15)**5.255 # air pressure at the pixel (i.e. at altitude) in hPa, using the international pressure equation
        P = Psurf * 1000. # air pressure at pixel location in dyn / cm2, which is hPa * 1000
        # calculation to get the constant of gravity at the pixel altitude, taking the air mass above into account
        dphi = math.radians(lat[x]) # latitude in radians
        cos2phi = math.cos(2*dphi)
        g0 = g0_45*(1-0.0026373*cos2phi+0.0000059*cos2phi**2)
        zs = 0.73737*z + 5517.56 # effective mass-weighted altitude
        g = g0-(0.0003085462 + 0.000000227*cos2phi)*zs + (0.00000000007254 + 0.0000000000001 * cos2phi)*zs**2 - (1.517E-17 + 6E-20*cos2phi)*zs**3
        # calculations to get the Rayeigh optical thickness
        factor = (P * AVO) / (m_a * g)
        for i in range(nbands):
            taur[(i,x)]=sigma[i] * factor

        # Calculate Rayleigh Phase function
        ts = math.radians(theta_s[x])  # sun zenith angle in radian
        cts = math.cos(ts)  # cosine of sun zenith angle
        sts = math.sin(ts)  # sinus of sun zenith angle
        tv = math.radians(theta_v[x])  # view zenith angle in radian
        ctv = math.cos(tv)  # cosine of view zenith angle
        stv = math.sin(tv)  # sinus of view zenith angle
        airmass[x] = 1/ctv + 1/ctv  # air mass
        # Rayleigh Phase function, 3 Fourier terms
        PR[0] =  3. * PA / 4.  * ( 1. + cts**2 * ctv**2 + (sts**2 * stv**2)/2.)+PB
        PR[1] = -3. * PA / 8.  * cts * ctv * sts * stv
        PR[2] =  3. * PA / 16. * sts**2 * stv**2
        # Calculate azimuth difference
        azs = math.radians(azi_s[x])
        azv = math.radians(azi_v[x])
        cosdeltaphi = math.cos(azv - azs)
        azidiff[x] = math.acos(cosdeltaphi) # azimuth difference in radian
        # Fourier components of multiple scattering
        for j in [0, 1, 2]:
            a[j] = interpn(gridGeometry, RayScattCoeffA[j, :,:], [theta_s[x], theta_v[x]], method='linear', bounds_error=False, fill_value=None)
            b[j] = interpn(gridGeometry, RayScattCoeffB[j, :,:], [theta_s[x], theta_v[x]], method='linear', bounds_error=False, fill_value=None)
            c[j] = interpn(gridGeometry, RayScattCoeffC[j, :,:], [theta_s[x], theta_v[x]], method='linear', bounds_error=False, fill_value=None)
            d[j] = interpn(gridGeometry, RayScattCoeffD[j, :,:], [theta_s[x], theta_v[x]], method='linear', bounds_error=False, fill_value=None)

        for i in range(nbands):
            # Fourier series, loop
            for j in [0,1,2]:
                # Rayleigh primary scattering
                rho_Rf[j] = (PR[j]/(4.0 * (cts + ctv))) * (1. - math.exp(-airmass[x]*taur[(i,x)]))
                # correction for multiple scattering
                rayMultiCorr[j] = a[j] + b[j] * taur[(i,x)] + c[j] * taur[(i,x)]**2 + d[j] * taur[(i,x)]**3
                rho_Rm[(j, i, x)]  = rho_Rf[j] * rayMultiCorr[j]
            # rho_Rm[(0, i, x)]  = rho_Rf[0]
            # rho_Rm[(1, i, x)]  = 0.
            # rho_Rm[(2, i, x)]  = 0.
            # Fourier sum to get the Rayleigh Reflectance
            rho_R[(i,x)] = rho_Rm[(0, i, x)] + 2.0 * rho_Rm[(1, i, x)] * math.cos(azidiff[x]) + 2. * rho_Rm[(2, i, x)] * math.cos(2.*azidiff[x])
            # complete the Rayleigh correction: see MERIS DPM PDF-p251 or DPM 9-16
            # polynomial coefficients tpoly0, tpoly1 and tpoly2 from MERIS LUT
            tRs = ((2./3. + cts) + (2./3. - cts) * math.exp(-taur[(i,x)]/cts))/(4./3. + taur[(i,x)])
            tR_thetaS[(i,x)] = tpoly[0] + tpoly[1] * tRs + tpoly[2] * tRs**2  # Rayleigh Transmittance sun - surface
            tRv = ((2./3. + ctv) + (2./3. - ctv) * math.exp(-taur[(i,x)]/ctv))/(4./3. + taur[(i,x)])
            tR_thetaV[(i,x)] = tpoly[0] + tpoly[1] * tRv + tpoly[2] * tRv**2  # Rayleigh Transmittance surface - sensor

            sARay[(i,x)] = rayAlb_f(taur[(i,x)])  # Rayleigh spherical albedo

            rho_toaR[(i,x)] = (reflectance[(i,x)] - rho_R[(i,x)])/(tR_thetaS[(i,x)]*tR_thetaV[(i,x)])  # toa reflectance corrected for Rayleigh scattering
            sphericalFactor[(i,x)] = 1.0/(1.0+sARay[(i,x)]*rho_toaR[(i,x)]) # factor used in the next equation to account for the spherical albedo
            rho_BRR[(i,x)] = rho_toaR[(i,x)]*sphericalFactor[(i,x)]  # top of aerosol reflectance, which is equal to bottom of Rayleigh reflectance

    # Write bands to product
    airmassBand.writePixels(0, y, width, 1, airmass)
    azidiffBand.writePixels(0, y, width, 1, azidiff)
    altBand.writePixels(0, y, width, 1, alt)

    for i in range(nbands):
        taurBand = raycorProduct.getBand("taur_"+str(i+1))
        taurBand.writePixels(0, y, width, 1, taur[i])
        rhorBand = raycorProduct.getBand("rRay_"+str(i+1))
        rhorBand.writePixels(0, y, width, 1, rho_R[i])
        rhorF1Band = raycorProduct.getBand("rRayF1_"+str(i+1))
        rhorF1Band.writePixels(0, y, width, 1, rho_Rm[0, i])
        rhorF2Band = raycorProduct.getBand("rRayF2_"+str(i+1))
        rhorF2Band.writePixels(0, y, width, 1, rho_Rm[1, i])
        rhorF3Band = raycorProduct.getBand("rRayF3_"+str(i+1))
        rhorF3Band.writePixels(0, y, width, 1, rho_Rm[2, i])
        rayTransSBand = raycorProduct.getBand("transSRay_"+str(i+1))
        rayTransSBand.writePixels(0, y, width, 1, tR_thetaS[i])
        rayTransVBand = raycorProduct.getBand("transVRay_"+str(i+1))
        rayTransVBand.writePixels(0, y, width, 1, tR_thetaV[i])
        sARayBand = raycorProduct.getBand("sARay_"+str(i+1))
        sARayBand.writePixels(0, y, width, 1, sARay[i])
        rtoaRBand = raycorProduct.getBand("rtoaRay_"+str(i+1))
        rtoaRBand.writePixels(0, y, width, 1, rho_toaR[i])
        rBRRBand = raycorProduct.getBand("rBRR_"+str(i+1))
        rBRRBand.writePixels(0, y, width, 1, rho_BRR[i])
        spfBand = raycorProduct.getBand("sphericalAlbedoFactor_"+str(i+1))
        spfBand.writePixels(0, y, width, 1, sphericalFactor[i])

    # Rayleigh calculation completed

raycorProduct.closeIO()

print("Done.")