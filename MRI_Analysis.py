from dicom_parser import Image
from dicom_parser.utils.siemens.csa.header import CsaHeader
import pydicom
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
import os
import pathlib
import matplotlib.pyplot as plt
import argparse
import MRI_Data
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
import numpy as np 
import math

import Visuals
import warnings

warnings.simplefilter('error')

# supply a folder containing the sequence images:
def GetADC_ls(mriListOrig):
    mriList = mriListOrig.copy()
    
    xLen, yLen = mriList[0][0].pixelData.shape #get the size of the pixel array

    bVals = []
    b_avg = 0
    
    for image in mriList[0]: #get the bvals list
        bVals.append(image.bVal)
        b_avg = b_avg + bVals[-1]
    num_bVals = len(bVals)   
    b_avg = b_avg / num_bVals

    slice_progress = 0

    for slice in mriList: #for each slice 
        adcImage = np.zeros((128,128))
        for x_idx in range(xLen):
            for y_idx in range(yLen): #iterate through all pixels
                #Now can perform the ADC calculation using least squares
                pixels = []
                log_S_avg = 0
                for image in slice:
                    pixels.append(image.pixelData[x_idx, y_idx])
                    log_S_avg += np.log(pixels[-1])
                log_S_avg = log_S_avg / num_bVals

                sum_numerator = 0
                sum_denominator = 0
                for idx in range(len(pixels)):
                    b = bVals[idx]
                    log_S = np.log(pixels[idx])
                    sum_numerator = sum_numerator + (b-b_avg)*(log_S - log_S_avg)
                    sum_denominator = sum_denominator + (b-b_avg)**2

                adc = -sum_numerator / sum_denominator
                adcImage[x_idx, y_idx] = adc
        adcMRI = MRI_Data.MRI_Data(adcImage)

        slice.append(adcMRI)     
        slice_progress = slice_progress + 1
        # plt.imshow(adcImage, cmap=plt.cm.bone)  
        # plt.show()
        print("Calculated ADC for " + str(slice_progress) + " slices")   
        

    #also organize the mris into a 4d array with third index being slice and 4th being the b val and also adc at end
    len_bVals_and_adc = len(mriList[0])
    numSlices = len(mriList)
    array = np.zeros((xLen, yLen, numSlices, len_bVals_and_adc))
    for slice_idx in range(len(mriList)):
        for bval_idx in range(len(mriList[slice_idx])):
            array[:,:, slice_idx, bval_idx] = mriList[slice_idx][bval_idx].pixelData



    return mriList, array
    print("Finished Calculating ADCs")

def GetADC_ls_single(bVals, signals):
    b_avg = 0
    for b in bVals:
        b_avg = b_avg + b
    b_avg = b_avg / len(bVals)    

    log_S_avg = 0
    for val in signals:
        if val == 0: 
            return math.nan
        log_S_avg += np.log(val)
    log_S_avg = log_S_avg / len(bVals)
    sum_numerator = 0
    sum_denominator = 0
    for idx in range(len(signals)):
        b = bVals[idx]
        log_S = np.log(signals[idx])
        sum_numerator = sum_numerator + (b-b_avg)*(log_S - log_S_avg)
        sum_denominator = sum_denominator + (b-b_avg)**2

    adc = -sum_numerator / sum_denominator
    return adc                


def Sort(path, addBValues = True, rel=True):
    # if path given is relative to working directory, make it absolute
    if rel == True:
        path = os.path.join(pathlib.Path(__file__).parent.absolute(), path)

    # Need to sort images according to slice location and b value.
    sliceLocations = []
    bValues = []
    totalImages = 0
    # first get list of slice locations:
    files = os.listdir(path)
    for file in files:
        if "5" not in file:
            continue
        absoluteFile = os.path.join(path, file)
        metaData = pydicom.dcmread(absoluteFile)
        try:  # Some images in the July 21 set don't have a slice location encoded
            sliceLocation = metaData[0x0020, 0x1041].value
        except:
            continue
        if sliceLocation not in sliceLocations:
            sliceLocations.append(sliceLocation)
        totalImages = totalImages + 1

    sliceLocations.sort()

    # Get the list of b values used in the series
    image = Image(absoluteFile).header.get("CSAImageHeaderInfo", parsed=False)
    csa_header_image = CsaHeader(image)

    raw_csa = Image(absoluteFile).header.get("CSASeriesHeaderInfo", parsed=False)

    csa_header = CsaHeader(raw_csa).parse()
    bValues = csa_header['Diffusion']['BValue'].copy()
    # Need to add b value 0 to the list
    bValues.insert(0, 0)

    # Now go through the files again, creating a list with a list for each slice location containing all the CTs at that location.
    mriList = []
    for i in range(len(sliceLocations)):  # Make an empty list for each slice location
        mriList.append([])

    for file in files:
        if "5" not in file:
            continue
        absoluteFile = os.path.join(path, file)
        metaData = pydicom.dcmread(absoluteFile)

        try:
            sliceLocation = metaData[0x0020, 0x1041].value
            timeStamp = float(metaData['ContentTime'].value)
        except:
            continue
        # get idx of slicelocation in slicelocations to match that index in mriList
        sliceIdx = [x for x in range(
            len(sliceLocations)) if sliceLocations[x] == sliceLocation][0]

        # now make a MRI object for the image
        mri = MRI_Data.MRI_Data(metaData.pixel_array,
                                file, timeStamp, sliceLocation)
        mriList[sliceIdx].append(mri)

        #now sort mris according to image time:
    for sliceList in mriList:
        sliceList.sort(key=lambda x: x.timeStamp)


    #The July 21 data now needs all mriList elements that dont have len(bvalues)+2 images in them removed (the first and last image are adc and something else)
    num_bVals = len(bValues)
    tempList = mriList.copy()
    mriList = []
    for sliceList in tempList:
        if len(sliceList) == num_bVals + 2:
            sliceList.pop(0) 
            sliceList.pop(len(bValues))
            mriList.append(sliceList.copy())

    #Now can resave the files with the b value added to the metadata

    processedPath = os.path.join(path, "Modified")
    if not os.path.isdir(processedPath):
        os.mkdir(processedPath)
    for sliceList in mriList:
        b_idx = 0
        for mri in sliceList:
            mri.bVal = bValues[b_idx]
            if addBValues == True:
                fileName = os.path.join(path,mri.fileName)
                data = pydicom.dcmread(fileName)
                saveName = processedPath + "/" + mri.fileName

                diffusionSequence = Sequence()

                diffusionDirectionality = Dataset()
                diffusionDirectionality.add_new([0x0018, 0x9075], "CS", "ISOTROPIC")
                diffusionSequence.append(diffusionDirectionality)

                diffusionGradientOrientation = Dataset()
                diffusionGradientOrientation.add_new([0x0018, 0x9089], "FD", [0.0, 0.0, 0.0])
                diffusionGradientDirection = Dataset()
                diffusionGradientDirection.add_new([0x0018, 0x9076], "SQ", Sequence([diffusionGradientOrientation]))
                diffusionSequence.append(diffusionGradientDirection)

                diffusionBValue = Dataset()
                diffusionBValue.add_new([0x0018, 0x9087], "FD", bValues[b_idx])
                diffusionSequence.append(diffusionBValue)

                diffusionAnisotropyType = Dataset()
                diffusionAnisotropyType.add_new([0x0018, 0x9147], "CS", "RELATIVE")
                diffusionSequence.append(diffusionAnisotropyType)


                data.add_new([0x0018, 0x9117], "SQ", diffusionSequence)

                # data[0x0018, 0x9087].value = bVals[b_idx]
                data.save_as(saveName)

            b_idx = b_idx + 1

    # for file in glob.glob(os.path.join(processedPath, "*")):
    #     dataa = pydicom.dcmread(file)

    # for image in mriList[7]:
    #     plt.imshow(image.pixelData, cmap=plt.cm.bone)  
    #     plt.show()
    print("Finished sorting")
    return mriList
def GetHessianAndGradient(bVals, pixels, f, pseudoD, D):
    derivative_f = 0
    derivative_ff = 0
    derivative_pseudoD = 0
    derivative_pseudoD_pseudoD = 0
    derivative_fpseudoD = 0
    derivative_pseudoDf = 0
    num_bVals = len(bVals)
    for i in range(num_bVals):
        c = np.exp(-bVals[i] * D)
        b = bVals[i]
        expon = np.exp(-b * pseudoD)
        signal = pixels[i]

        derivative_f = derivative_f +  2.0 * (signal - f*expon - (1.0-f)*c) * (-expon + c)
        derivative_pseudoD = derivative_pseudoD + 2.0 * ( signal - f*expon - (1.0-f)*c ) * (b*f*expon)

        derivative_ff = derivative_ff + 2.0 * (c - expon)**2
        derivative_pseudoD_pseudoD = derivative_pseudoD_pseudoD +  2.0 * (b*f*expon) - 2.0 * (signal - f*expon-(1.0-f)*c)*(b**2*f*expon)

        derivative_fpseudoD = derivative_fpseudoD +  (2.0 * (c - expon)*b*f*expon) + (2.0 * (signal - f*expon - (1.0-f)*c) * b*expon )
        derivative_pseudoDf = derivative_pseudoDf + (2.0 * (b*f*expon)*(-expon + c)) + (2.0*(signal - f*expon - (1.0-f)*c)*(b*expon))

    H = np.zeros((2,2))
    H[0,0] = derivative_ff
    H[0,1] = derivative_fpseudoD
    H[1,0] = derivative_pseudoDf
    H[1,1] = derivative_pseudoD_pseudoD    
    
    grad = np.zeros((2,1))
    grad[0] = derivative_f
    grad[1] = derivative_pseudoD

    return H, grad

def GetBiExp_Params(mriListOrig, numIterations=100):
    '''This function will use the a segmented approach with Marquardts method of squared residuals minimization to fit the signal to a biexponential 
    '''
    mriList = mriListOrig.copy()
    xLen, yLen = mriList[0][0].pixelData.shape #get the size of the pixel array

    bVals = []
    
    for image in mriList[0]: #get the bvals list
        bVals.append(image.bVal)
    num_bVals = len(bVals)   

    slice_progress = 0

    for slice in mriList: #for each slice 
        paramImage = np.zeros((3,128,128))
        for x_idx in range(xLen):
            for y_idx in range(yLen): #iterate through all pixels
                #Now can perform the ADC calculation using least squares
                pixels = []
                log_S_avg = 0
                S0 = slice[0].pixelData[x_idx, y_idx]
                for image in slice:
                    pixels.append(image.pixelData[x_idx, y_idx] / S0)
                #First we use a cutoff of b values greater than 200 to have signals of the form S(b) = S(0) * exp(-b D) to obtain the diffusion coefficient
                #make a vector for the high b values and their signals
                bValsH = []
                signalsH = []
                for idx in range(len(pixels)):
                    b = bVals[idx]
                    if b > 200: 
                        bValsH.append(b)
                        signalsH.append(pixels[idx])

                #Now use least squares regression to obtain D from the high b value signals
                D = GetADC_ls_single(bValsH, signalsH)    
                if math.isnan(D):
                    paramImage[0, x_idx, y_idx] = 0
                    paramImage[1, x_idx, y_idx] = 0
                    paramImage[2, x_idx, y_idx] = 0
                    continue

                if D < 0:
                    D = 0
                #Now we can use this D value and fit f and D* with Marquardts method
                #Cost function is sum (Signal_i - (fexp(-bD*) + (1-f)exp(-bD)) )^2
                lam = 50
                pseudoD = 10.0 * D
                f = 0.5

                #Get the current cost
                cost = 0
                for i in range(len(bVals)):
                    bTemp = bVals[i]        
                    signalTemp = pixels[i]
                    cost = cost +  ((signalTemp) - f*math.exp(-bTemp * pseudoD) - (1-f)*math.exp(-bTemp * D))**2
                     
                for i in range(numIterations):
                    H, grad = GetHessianAndGradient(bVals, pixels, f, pseudoD, D)
                    H = H + np.identity(2) * lam 
                    inverse = np.linalg.inv(H)

                    #now update params
                    newParams = np.matmul(inverse, grad)
                    new_f = f - newParams[0][0]
                    new_pseudoD = pseudoD - newParams[1][0]

                    if new_f < 0: 
                        new_f = 0
                    elif new_f > 1:
                        new_f = 1

                    if new_pseudoD < 0:
                        new_pseudoD = 0 

                    newCost = 0     
                    for i in range(len(bVals)):
                        bTemp = bVals[i]        
                        signalTemp = pixels[i]
                        newCost = newCost +  ((signalTemp) - new_f*math.exp(-bTemp * new_pseudoD) - (1-new_f)*math.exp(-bTemp * D))**2       
                    if newCost < cost:
                        cost = newCost
                        lam = lam * 0.8
                        f = new_f
                        pseudoD = new_pseudoD
                    else:
                        lam = lam * 2    
                        

                paramImage[0, x_idx, y_idx] = D
                paramImage[1, x_idx, y_idx] = pseudoD
                paramImage[2, x_idx, y_idx] = f
        D_MRI = MRI_Data.MRI_Data(paramImage[0,:,:])
        pseudoD_MRI = MRI_Data.MRI_Data(paramImage[1,:,:])
        f_MRI = MRI_Data.MRI_Data(paramImage[2,:,:])

        slice.append(D_MRI)
        slice.append(pseudoD_MRI)
        slice.append(f_MRI)

        slice_progress = slice_progress + 1
        # plt.imshow(adcImage, cmap=plt.cm.bone)  
        # plt.show()
        print("Calculated ADC for " + str(slice_progress) + " slices")   
        
    #also organize the mris into a 4d array with third index being slice and 4th being the b val and also param images at end
    len_bVals_and_params = len(mriList[0])
    numSlices = len(mriList)
    array = np.zeros((xLen, yLen, numSlices, len_bVals_and_params))
    for slice_idx in range(len(mriList)):
        for bval_idx in range(len(mriList[slice_idx])):
            array[:,:, slice_idx, bval_idx] = mriList[slice_idx][bval_idx].pixelData



    return mriList, array
    print("Finished Calculating ADCs")    



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="IVIM analysis by Caleb Sample"
    )
    relPath = "July21/DICOM/21072609/05010002"
    mriList = Sort(relPath, addBValues=False)

    # mris, array = GetADC_ls(mriList)
    # Visuals.SliderPlot(array, ["ADC map"])

    mris, array = GetBiExp_Params(mriList, numIterations=30)
    Visuals.SliderPlot(array, ["Diffusion Coefficient (mm^2/2)", "Pseudodiffusion Coefficient (mm^2/s)", "Perfusion Fraction"])
    
