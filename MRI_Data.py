#This class will old the necessary MRI data needed. 

class MRI_Data():
    def __init__(self, pixelData, fileName=None, timeStamp=None, sliceLocation=None):
        self.pixelData = pixelData
        self.fileName = fileName
        self.timeStamp = timeStamp
        self.sliceLocation = sliceLocation
        self.bVal = None

    def SetBValue(self, b):
        self.bVal = b

