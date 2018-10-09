#!/usr/bin/env python
"""Convert a numpy array (3D) to a dicom CT series"""

import os
from datetime import datetime
from math import floor
import pydicom
import numpy as np

CTIMAGE_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.2"
def make_dicom_boilerplate(SeriesInstanceUID=None, StudyInstanceUID=None, FrameOfReferenceUID=None):
    # Populate required values for file meta information
    file_meta = pydicom.dataset.Dataset()
    file_meta.FileMetaInformationGroupLength = 204
    file_meta.FileMetaInformationVersion = b'00\01'
    file_meta.MediaStorageSOPClassUID = CTIMAGE_SOP_CLASS_UID
    file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
    file_meta.ImplementationClassUID = '2.25.229451600072090404564844894284998027179' #arbitrary specific to this library
    file_meta.ImplementationVersionName = "PyMedImage"
    file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
    ds = pydicom.dataset.Dataset()
    ds.preamble = b"\0" * 128
    ds.file_meta = file_meta
    ds.is_little_endian = True
    ds.is_implicit_VR = True

    datestr = datetime.now().strftime('%Y%m%d')
    timestr = datetime.now().strftime('%H%M%S')
    ds.ContentDate = datestr
    ds.ContentTime = timestr
    ds.StudyDate = datestr
    ds.StudyTime = timestr
    ds.PatientID = 'ANON0001'
    ds.StudyID = 'ANON0001'
    ds.SeriesNumber = '0001'
    ds.StudyDate = datestr
    ds.StudyTime = timestr
    ds.AccessionNumber = ''
    ds.ReferringPhysiciansName = ''
    ds.PatientName = 'ANON0001'
    ds.PatientSex = ''
    ds.PatientAge = ''
    ds.PatientBirthDate = ''
    ds.PatientOrientation = 'LA'
    ds.PatientPosition = 'HFS'
    ds.ImagePositionPatient = [0, 0, 0]
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.InstanceNumber = 1
    ds.StudyInstanceUID = pydicom.uid.generate_uid() if StudyInstanceUID is None else StudyInstanceUID
    ds.SeriesInstanceUID = pydicom.uid.generate_uid() if SeriesInstanceUID is None else SeriesInstanceUID
    ds.FrameOfReferenceUID = pydicom.uid.generate_uid() if FrameOfReferenceUID is None else FrameOfReferenceUID
    ds.SOPInstanceUID = pydicom.uid.generate_uid()
    ds.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL']
    ds.Modality = ''
    ds.SOPClassUID = CTIMAGE_SOP_CLASS_UID
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.RescaleIntercept = 0
    ds.RescaleSlope = 1.0
    ds.KVP = ''
    ds.AcquisitionNumber = 1
    ds.PixelRepresentation = 0
    ds.SliceLocation = 0.0
    ds.Rows = 0
    ds.Columns = 0
    ds.PixelSpacing = [1.0, 1.0]
    ds.SliceThickness = 1.0
    ds.Units = 'HU'
    ds.RescaleType = 'HU'
    return ds

class FrameOfReference:
    """Defines a dicom frame of reference to which BaseVolumes can be conformed for fusion of pre-registered
    image data
    """
    def __init__(self, start=None, spacing=None, size=None, UID=None):
        """Define a dicom frame of reference

        Args:
            start    -- (x,y,z) describing the start of the FOR (mm)
            spacing  -- (x,y,z) describing the spacing of voxels in each direction (mm)
            size     -- (x,y,z) describing the number of voxels in each direction (integer)
            UID      -- dicom FrameOfReferenceUID can be supplied to support caching in BaseVolume

        Standard Anatomical Directions Apply:
            x -> increasing from patient right to left
            y -> increasing from patient anterior to posterior
            z -> increasing from patient inferior to superior
        """
        self.start = start
        self.spacing = spacing
        self.size = size
        self.UID = UID

class BaseVolume:
    """Defines basic storage for volumetric voxel intensities within a dicom FrameOfReference
    """
    def __init__(self):
        """Entrypoint to class, initializes members
        """
        self.data = None
        self.init_object = None
        self.frameofreference = None
        self.modality = None
        self.feature_label = None

    @property
    def nslices(self):
        if len(self.frameofreference.size)>=3:
            return self.frameofreference.size[-1]
        else:
            return 1

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, v):
        self._data = v


    @classmethod
    def fromArray(cls, array, frameofreference=None):
        """Constructor: from a numpy array and FrameOfReference object

        Args:
            array             -- numpy array
            frameofreference  -- FrameOfReference object
        """
        # ensure array matches size in frameofreference
        self = cls()
        if array.ndim == 2:
            array = np.atleast_3d(array)
        if frameofreference is not None:
            self.data = array.reshape(frameofreference.size[::-1])
            self.frameofreference = frameofreference
        else:
            self.data = array
            self.frameofreference = FrameOfReference((0,0,0), (1,1,1), (*array.shape[::-1], 1))

        return self

    def toDicom(self, dname, fprefix=''):
        SeriesInstanceUID   = pydicom.uid.generate_uid()
        StudyInstanceUID    = pydicom.uid.generate_uid()
        FrameOfReferenceUID = pydicom.uid.generate_uid()
        min_val = np.min(self.data)
        for i in range(self.frameofreference.size[2]):
            ds = make_dicom_boilerplate(SeriesInstanceUID, StudyInstanceUID, FrameOfReferenceUID)
            ds.SliceThickness = self.frameofreference.spacing[2]
            ds.PixelSpacing = list(self.frameofreference.spacing[:2])
            ds.SliceLocation = self.frameofreference.start[2] + i*self.frameofreference.spacing[2]
            ds.ImagePositionPatient = [*self.frameofreference.start[:2], ds.SliceLocation]
            ds.Columns = self.frameofreference.size[0]
            ds.Rows = self.frameofreference.size[1]
            ds.AcquisitionNumber = i+1
            ds.Modality = self.modality if self.modality is not None else ''
            ds.DerivationDescription = self.feature_label if self.feature_label is not None else ''
            ds.PixelData = ((self.data[i, :, :]-min_val).flatten().astype(np.uint16)).tostring()
            ds.RescaleSlope = 1.0
            ds.RescaleIntercept = floor(min_val)
            ds.PixelRepresentation = 0 # unsigned integers
            os.makedirs(dname, exist_ok=True)
            ds.save_as(os.path.join(dname, '{}{:04d}.dcm'.format(fprefix, i)))


#############################################################################################################

# DEFINE GEOMETRY
if __name__ == "__main__":
    voxel = (1,1,1) # [unit: mm]
    size = np.ceil(np.array([101,101,300])/voxel).astype(int)
    frame = FrameOfReference(start=(0,0,0), spacing=voxel, size=size[::-1])
    arr = np.ones(size)
    vol = BaseVolume.fromArray(arr, frame)
    print('volume shape: '+str(arr.shape))
    vol.toDicom('water_phantom')
