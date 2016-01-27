import sys
import numpy as np
from struct import *
import math

def readRayADF():

    rayADF = {}
    with open('C:\\Users\\carsten\\Dropbox\\Carsten\\Tagesordner\\20150814\\Rayleigh-Correction-Processor\\ADF\\MER_ATP_AXVACR20091026_144725_20021224_121445_20200101_000000','rb') as f:
        # skip MPH  = 1247 bytes
        f.read(1247)
        # skip SPH = 2898 bytes
        f.read(2898)

        # GADS general = 2772 bytes
        # tR Rayleigh Transmittance coefficients
        tR = np.zeros(3, dtype=np.float32)
        str_tR = f.read(3*4)
        tR = unpack('>3f',str_tR)
        print('Rayleigh transmission ', tR)
        rayADF['tR']=tR
        # tau_R tabulated values = 68bytes
        dimtauR = 17
        taurTab = np.zeros(dimtauR, dtype=np.float32)
        str_taurTab = f.read(dimtauR*4)
        taurTab = unpack('>17f', str_taurTab)
        rayADF['taurTab'] = taurTab
        print('tau R tabulated values ', taurTab)
        # lam tabulated values = 15 floats = 60 bytes
        lam = np.zeros(15, dtype=np.float32)
        str_lam = f.read(60)
        lam = unpack('>15f',str_lam)
        rayADF['lam'] = lam
        print('lam ', lam)
        # skip some variables that we don't need (92+52+100)
        f.read(92+52+100)
        # theta zenith angles used for Rayleigh scattering
        thetaL = np.zeros(12, dtype=np.uint32)
        str_theta = f.read(48)
        thetaL = unpack('>12L',str_theta)
        theta = np.zeros(12, dtype=np.float32)
        for i in range(12):
            theta[i]=thetaL[i]/1e6
        rayADF['theta'] = theta
        print('theta: ', theta)
        # thetaS and thetaV index combinations of Rayleigh scattering terms
        z = np.zeros(156, dtype=np.byte)
        theta_indices = np.zeros((78,2), dtype=np.byte)
        str_z = f.read(156)
        z = unpack('>156B',str_z)
        n = 0
        for i in range(78):
            theta_indices[(i,0)] = z[n]
            theta_indices[(i,1)] = z[n+1]
            n += 2
        rayADF['theta_indices'] = theta_indices
        # print('theta indices ', theta_indices)
        # Constants for Rayleigh Phase function
        rPhaseCoeff = np.zeros(2, dtype=np.float32)
        str_z = f.read(8)
        rPhaseCoeff = unpack('>2f', str_z)
        rayADF['rPhaseCoeff'] = rPhaseCoeff
        print('Constants for Rayleigh Phase function ', rPhaseCoeff)
        # skip airmass tabulated values
        f.read(24)
        # reference wavelengths for O2 transmission
        lamO2 = np.zeros(21, dtype = np.float32)
        str_lamO2 = f.read(21*4)
        lamO2 = unpack('>21f', str_lamO2)
        rayADF['lamO2'] = lamO2
        print( 'reference wavelengths for O2 transmission ', lamO2)
        # normalised radiance values at 779nm
        Lnorm779 = np.zeros(25, dtype=np.float32)
        str_Lnorm779 = f.read(25*4)
        Lnorm779 = unpack('>25f', str_Lnorm779)
        rayADF['lnorm779'] = Lnorm779
        print('normalised radiance values at 779nm ', Lnorm779)
        # thetaS values for O2 transmissions
        thetaS_O2 = np.zeros(15, dtype=np.float32)
        str_z = f.read(15*4)
        z = unpack('>15L', str_z)
        for i in range(15):
            thetaS_O2[i] = z[i]/1E6
        rayADF['thetaS_O2'] = thetaS_O2
        print('thetaS values for O2 transmissions ', thetaS_O2)
        # thetaV values for O2 transmissions
        thetaV_O2 = np.zeros(10, dtype=np.float32)
        str_z = f.read(10*4)
        z = unpack('>10L', str_z)
        for i in range(10):
            thetaV_O2[i] = z[i]/1E6
        rayADF['thetaV_O2'] = thetaV_O2
        print('thetaV values for O2 transmissions ', thetaV_O2)
        # delta phi values for O2 transmissions
        dphi_O2 = np.zeros(19, dtype=np.float32)
        str_z = f.read(19*4)
        z = unpack('>19L', str_z)
        for i in range(19):
            dphi_O2[i] = z[i]/1E6
        rayADF['dphi_O2'] = dphi_O2
        print('delta phi values for O2 transmissions ', dphi_O2)
        # skip rest of GADS general = 2772 - 12 - 68 - 60 - 92 - 52 - 100 - 48 - 156 - 8 - 24 - 84 - 100 - 60 - 40 - 76
        f.read(2772 - 12 - 68 - 60 - 92 - 52 - 100 - 48 - 156 - 8 - 24 - 84 - 100 - 60 - 40 - 76)

        # skip GADS optical thickness = 300 bytes
        f.read(300)

        # skip GADS H2O transmission = 912 byte
        f.read(912)

        # GADS Rayleigh Scattering function = 3744 bytes
        # Rayleigh scattering Fourier term 1st polynomial coeff. (thetaS x thetaV, order)
        dimTheta = 12
        dimThetaS = dimThetaV = dimTheta
        RayScattCoeffA = np.zeros( (dimThetaS, dimThetaV, 3), dtype=np.float32)
        RayScattCoeffB = np.zeros( (dimThetaS, dimThetaV, 3), dtype=np.float32)
        RayScattCoeffC = np.zeros( (dimThetaS, dimThetaV, 3), dtype=np.float32)
        RayScattCoeffD = np.zeros( (dimThetaS, dimThetaV, 3), dtype=np.float32)
        rSA = np.zeros( (78, 3), dtype=np.float32)
        rSB = np.zeros( (78, 3), dtype=np.float32)
        rSC = np.zeros( (78, 3), dtype=np.float32)
        rSD = np.zeros( (78, 3), dtype=np.float32)
        rS = np.zeros( (4,78, 3), dtype=np.float32)
        z = np.zeros(78*4*3, dtype=np.float32)
        str_z = f.read(4*78*4*3)
        z = unpack('>936f', str_z)

        # data are stored in 3 records, for each of the 3 Fourier terms.
        # for each Fourier term, there are 4 coefficients. In the equation later on, these are called RayScattCoeffA,
        # RayScattCoeffB, ..C and ..D.
        # Reading is implemented here in this way
        # - 1st step: read the whole dataset into one linear array. This has already been done with the code above. This
        #   big linear array is called z
        # - 2nd step: split the big linear array into 3 arrays, each for one of the 4 coefficients A, B, C and D
        #   note: each of these arrays has 2 dimensions: 78 combinations thetaS x thetaV, and 3 Fourier terms
        n = 0
        for i in range(3):
            for j in range(4):
                rS[j, :, i] = z[78*(n+0):78*(n+1)]
                n += 1
        z0 = np.zeros(78*4, dtype=np.float32)
        z0 = z[0:78*4]
        z1 = np.zeros(78*4, dtype=np.float32)
        z1 = z[78*4: 78*4*2]
        z2 = np.zeros(78*4, dtype=np.float32)
        z2 = z[78*4*2:78*4*3]

        fourierTerm0A = np.zeros(78, dtype=np.float32)
        fourierTerm0B = np.zeros(78, dtype=np.float32)
        fourierTerm0C = np.zeros(78, dtype=np.float32)
        fourierTerm0D = np.zeros(78, dtype=np.float32)

        fourierTerm1A = np.zeros(78, dtype=np.float32)
        fourierTerm1B = np.zeros(78, dtype=np.float32)
        fourierTerm1C = np.zeros(78, dtype=np.float32)
        fourierTerm1D = np.zeros(78, dtype=np.float32)

        fourierTerm2A = np.zeros(78, dtype=np.float32)
        fourierTerm2B = np.zeros(78, dtype=np.float32)
        fourierTerm2C = np.zeros(78, dtype=np.float32)
        fourierTerm2D = np.zeros(78, dtype=np.float32)

        n=0
        for i in range(78):
            fourierTerm0A[i] = z0[n]
            fourierTerm0B[i] = z0[n+1]
            fourierTerm0C[i] = z0[n+2]
            fourierTerm0D[i] = z0[n+3]
            n=n+4

        n=0
        for i in range(78):
            fourierTerm1A[i] = z1[n]
            fourierTerm1B[i] = z1[n+1]
            fourierTerm1C[i] = z1[n+2]
            fourierTerm1D[i] = z1[n+3]
            n=n+4

        n=0
        for i in range(78):
            fourierTerm2A[i] = z2[n]
            fourierTerm2B[i] = z2[n+1]
            fourierTerm2C[i] = z2[n+2]
            fourierTerm2D[i] = z2[n+3]
            n=n+4

        for i in range(78):
            i_theta_s = theta_indices[i,0]
            i_theta_v = theta_indices[i,1]

            RayScattCoeffA[i_theta_s, i_theta_v,0]= fourierTerm0A[i]
            RayScattCoeffA[i_theta_v, i_theta_s,0]= fourierTerm0A[i]

            RayScattCoeffA[i_theta_s, i_theta_v,1]= fourierTerm1A[i]
            RayScattCoeffA[i_theta_v, i_theta_s,1]= fourierTerm1A[i]

            RayScattCoeffA[i_theta_s, i_theta_v,2]= fourierTerm2A[i]
            RayScattCoeffA[i_theta_v, i_theta_s,2]= fourierTerm2A[i]

            RayScattCoeffB[i_theta_s, i_theta_v,0]= fourierTerm0B[i]
            RayScattCoeffB[i_theta_v, i_theta_s,0]= fourierTerm0B[i]

            RayScattCoeffB[i_theta_s, i_theta_v,1]= fourierTerm1B[i]
            RayScattCoeffB[i_theta_v, i_theta_s,1]= fourierTerm1B[i]

            RayScattCoeffB[i_theta_s, i_theta_v,2]= fourierTerm2B[i]
            RayScattCoeffB[i_theta_v, i_theta_s,2]= fourierTerm2B[i]

            RayScattCoeffC[i_theta_s, i_theta_v,0]= fourierTerm0C[i]
            RayScattCoeffC[i_theta_v, i_theta_s,0]= fourierTerm0C[i]

            RayScattCoeffC[i_theta_s, i_theta_v,1]= fourierTerm1C[i]
            RayScattCoeffC[i_theta_v, i_theta_s,1]= fourierTerm1C[i]

            RayScattCoeffC[i_theta_s, i_theta_v,2]= fourierTerm2C[i]
            RayScattCoeffC[i_theta_v, i_theta_s,2]= fourierTerm2C[i]

            RayScattCoeffD[i_theta_s, i_theta_v,0]= fourierTerm0D[i]
            RayScattCoeffD[i_theta_v, i_theta_s,0]= fourierTerm0D[i]

            RayScattCoeffD[i_theta_s, i_theta_v,1]= fourierTerm1D[i]
            RayScattCoeffD[i_theta_v, i_theta_s,1]= fourierTerm1D[i]

            RayScattCoeffD[i_theta_s, i_theta_v,2]= fourierTerm2D[i]
            RayScattCoeffD[i_theta_v, i_theta_s,2]= fourierTerm2D[i]


        # - 3rd step: construct a 2D matrix from the 78 thetaS x thetaV combinations
        #   this is described in IODD section 6.11.4, page 6.11.-3 (PDF page 167)
        # k = 0
        # for i_theta_s in range(dimTheta):
        #     for i_theta_v in range(i_theta_s,dimTheta):
        #         for fourier in range(3):
        #             RayScattCoeffA[i_theta_s, i_theta_v, fourier] = rS[0, k, fourier]
        #             RayScattCoeffB[i_theta_s, i_theta_v, fourier] = rS[1, k, fourier]
        #             RayScattCoeffC[i_theta_s, i_theta_v, fourier] = rS[2, k, fourier]
        #             RayScattCoeffD[i_theta_s, i_theta_v, fourier] = rS[3, k, fourier]
        #             if (i_theta_v != i_theta_s):
        #                 RayScattCoeffA[i_theta_v, i_theta_s, fourier] = RayScattCoeffA[i_theta_s, i_theta_v, fourier]
        #                 RayScattCoeffB[i_theta_v, i_theta_s, fourier] = RayScattCoeffB[i_theta_s, i_theta_v, fourier]
        #                 RayScattCoeffC[i_theta_v, i_theta_s, fourier] = RayScattCoeffC[i_theta_s, i_theta_v, fourier]
        #                 RayScattCoeffD[i_theta_v, i_theta_s, fourier] = RayScattCoeffD[i_theta_s, i_theta_v, fourier]
        #         k += 1

        # print ' Rayleigh Scattering Coefficient A', RayScattCoeffA
        # print ' Rayleigh Scattering Coefficient B', RayScattCoeffB
        # print ' Rayleigh Scattering Coefficient C', RayScattCoeffC
        # print ' Rayleigh Scattering Coefficient D', RayScattCoeffD
        rayADF['RayScattCoeffA'] = RayScattCoeffA
        rayADF['RayScattCoeffB'] = RayScattCoeffB
        rayADF['RayScattCoeffC'] = RayScattCoeffC
        rayADF['RayScattCoeffD'] = RayScattCoeffD
        # skip GADS Rayleigh Spherical Albedo = 68 bytes
        rayAlbLUT = np.zeros(dimtauR, dtype=np.float32)  # Rayleigh spherical albedo LUT
        str_rayAlbLUT = f.read(dimtauR*4)
        rayAlbLUT = unpack('>17f', str_rayAlbLUT)
        rayADF['rayAlbLUT'] = rayAlbLUT
        # print 'Rayleigh spherical albedo LUT ', rayAlbLUT

        # skip GADS O2 transmission around 779nm = 5985000 bytes
        f.read(5985000)

    return rayADF

    # # now write all the variable that we need into a json file
    # import json
    # with open('auxdata_Rayleigh_v1.txt','w') as RayleighADF:
    #     json.dump(lam, RayleighADF, ensure_ascii=False)
    #     json.dump(taurTab, RayleighADF, ensure_ascii=False)
