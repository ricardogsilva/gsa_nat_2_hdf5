#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
Script's doctring goes here.
"""

import sys
from struct import unpack
import numpy as np
import tables

# TODO
#
#   - Add the relevant attributes to the HDF5 file
#   - Remove the outPath argument from GSA_nat.to_hdf5() and add an outDir.
#   The output name can be set by the class, based on the input nat file.
#   - Remove the unneeded attributes from the HDF5 file (the ones 
#   automatically added by pytables).
#   - Add this class to the processing lines and incorporate a call to it
#   in the G2GSA class. It should be called after the run() method.

class GSANat(object):

    MAIB_HEADER = {
            'start' : 0, 
            'nBytes' : 136
            }

    PRODUCT_HEADER = {
            'start' : MAIB_HEADER['start'] + MAIB_HEADER['nBytes'], 
            'nBytes' : 572
            }

    LINES_PART = {
            'start' : PRODUCT_HEADER['start'] + PRODUCT_HEADER['nBytes'], 
            'nBytes' : None
            }

    datasets = [
            ['BHRiso', None],
            ['DHR', None],
            ['QCFlag', None],
            ['NumSol', None],
            ['NumInputSlot', None],
            ['NumberSlotASM', None],
            ['IdxSurface', None],
            ['IdxOptical', None],
            ['R_0', None],
            ['Err_R_0', None],
            ['DHRErr_10D', None] ,
            ['NumDays', None] ,
            ['BestDayNbr', None],
            ['Ch2ASM', None],
            ['Ch2DCP', None],
            ['ProbIdx', None],
            ['DHRError', None],
            ['Err_K', None],
            ['Err_T', None],
            ['Err_Tau', None],
            ['TauAvg', None],
            ['TauStdErrAvg', None],
            ['RadNoise', None],
            ]
                
    def __init__(self, filePath):
        fh = open(filePath, 'rb')
        self.natHeader = self._decode_header(fh)
        self._decode_data(fh)
        fh.close()

    def _decode_header(self, fh):
        fh.seek(self.PRODUCT_HEADER['start']) # skip the MAIB header
        headerFields = [
                # Prodct Identification
                ('iYear1' ,                      'i'),
                ('iDayinYear1' ,                 'i'),
                ('iYear2' ,                      'i'),
                ('iDayinYear2' ,                 'i'),
                ('iProductTime' ,                'i'),
                ('shMeteosatID' ,                '4s'),
                ('szProductName' ,               '32s'),
                ('szMSAVersion' ,                '8s'),
                ('iWorkAreaFirstPxl' ,           'i'),
                ('iWorkAreaLastPxl' ,            'i'),
                ('iWorkAreaFirstRow' ,           'i'),
                ('iWorkAreaLastRow' ,            'i'),
                ('IIPSDataVersion' ,             'i'),
                ('iIPSAlgoVersion' ,             'i'),
                ('ICalVersion' ,                 '12s'),
                ('iNbrLine' ,                    'i'),
                ('iSpare' ,                      'i'),
                # Product set up parameters
                ('fl_DAM_WaterReflectanceTH' ,   'f'),
                ('fl_DAM_CloudForSureTH' ,       'f'),
                ('fl_DCP_PSmooth' ,              'f'),
                ('fl_NominalSSP' ,               'f'),
                ('iNbrOptIdx' ,                  'i'),
                ('iNbrSrfIdx' ,                  'i'),
                ('flOptThick' ,                  '10f'),
                ('fl_RPV_K' ,                    '10f'),
                ('fl_RPV_Theta' ,                '10f'),
                ('iNbrProbTresh' ,               'i'),
                ('fl_ProbTresh' ,                '20f'),
                ('flProbabilityPAlpha' ,         'f'),
                ('flAutoCorrCoeff' ,             'f'),
                ('flChi2AvgVal' ,                'f'),
                ('flScalingDHR' ,                'f'),
                ('flScalingBHRiso' ,             'f'),
                ('flScalingR0' ,                 'f'),
                ('flScalingAvgTau' ,             'f'),
                ('flScalingChi2DCP' ,            'f'),
                ('flScalingChi2ASM' ,            'f'),
                ('flScalingRadRelErr' ,          'f'),
                ('flScalingErr_DHR' ,            'f'),
                ('flScaling_Err_R0' ,            'f'),
                ('flScaling_Err_K' ,             'f'),
                ('flScaling_Err_Theta' ,         'f'),
                ('flScaling_Err_Tau' ,           'f'),
                ('Ispare' ,                      '2i'),
                # Product Statistics
                ('iActualNbrDay' ,               'i'),
                ('lNbrProduct' ,                 'i'),
                ('flPercentProd' ,               'f'),
                ('fl_AvgRelRadErr' ,             'f'),
                ('flAveAvailableSlot' ,          'f'),
                ('flAveProcessedSlot' ,          'f'),
                ('flAveValidPixels' ,            'f'),
                ('flAveNWeakSolutions' ,         'f'),
                ('flAveFreqDubious' ,            'f'),
                ('flAveNSolutions' ,             'f'),
                ('flAveOpt' ,                    'f'),
                ('FlAveActualSspLat' ,           'f'),
                ('flAveActualSspLon' ,           'f'),
                ('flAveProbability' ,            'f'),
                ('flAveDHR30' ,                  'f'),
                ('flAveDHR30RelErr' ,            'f'),
                ('spare' ,                       '26i'),
                ]
        headerBytes = fh.read(self.PRODUCT_HEADER['nBytes'])
        formatString = ''.join([f[1] for f in headerFields])
        values = unpack(formatString, headerBytes)
        return dict(zip([f[0] for f in headerFields], values))

    def _decode_data(self, fh):
        fh.seek(self.LINES_PART['start'])
        nLines, nCols = self._get_dimensions()
        lineDtype = np.dtype([
            ('header', 'i4', 6), # the line header is repeated twice
            ('data', 'u1', (nCols, len(self.datasets))),
            ])
        dataArr = np.fromfile(fh, dtype=lineDtype)
        for index, dsList in enumerate(self.datasets):
            self.datasets[index][1] = dataArr['data'][:,:,index]

    def _get_dimensions(self):
        '''return a tuple with number of lines and columns.'''

        nLines = self.natHeader['iWorkAreaLastPxl'] - \
                 self.natHeader['iWorkAreaFirstPxl'] + 1
        nCols = self.natHeader['iWorkAreaLastRow'] - \
                self.natHeader['iWorkAreaFirstRow'] + 1
        return nLines, nCols

    def to_hdf5(self, outPath):
        '''
        Convert the GSA_NAT file to HDF5.
        '''

        h5f = tables.openFile(outPath, mode='w', title='GSA')
        for dsList in self.datasets:
            name, arr = dsList
            h5f.createArray('/', name, arr)
        h5f.close()
