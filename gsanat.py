#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
Script's doctring goes here.
"""

import sys
from struct import unpack
import numpy as np
import tables
import re

# TODO
#
#   - Add the relevant attributes to the HDF5 file
#       - measurement units (ds)
#       - scaling factor (ds)
#       - missing value (ds)
#       - subsatellite point (root)
#       - ...
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

    DATASETS = [
            # [datasetName, empty placeholder for the data, scaling factor]
            ['BHRiso', None, 254],
            ['DHR', None, 254],
            ['QCFlag', None, 1],
            ['NumSol', None, 1],
            ['NumInputSlot', None, 1],
            ['NumberSlotASM', None, 1],
            ['IdxSurface', None, 1],
            ['IdxOptical', None, 1],
            ['R_0', None, 254],
            ['Err_R_0', None, 254],
            ['DHRErr_10D', None, 1] , # is the scaling factor correct?
            ['NumDays', None, 1] ,
            ['BestDayNbr', None, 1],
            ['Ch2ASM', None, 127],
            ['Ch2DCP', None, 254],
            ['ProbIdx', None, 1],
            ['DHRError', None, 1],
            ['Err_K', None, 254],
            ['Err_T', None, 254],
            ['Err_Tau', None, 254],
            ['TauAvg', None, 1], # is the scaling factor correct?
            ['TauStdErrAvg', None, 1], # is the scaling factor correct?
            ['RadNoise', None, 1], # is the scaling factor correct?
            ]

    __HDF5_ATTRIBUTES = {
            'SAF'                               : 'GEOLAND2',
            'CENTRE'                            : 'IM-PT',
            'ARCHIVE_FACILITY'                  : 'IM-PT',
            'PRODUCT'                           : 'GSA',
            'PARENT_PRODUCT_NAME'               : ['-', '-', '-', '-'],
            'SPECTRAL_CHANNEL_ID'               : 0,
            'PRODUCT_ALGORITHM_VERSION'         : None, # szMSAVersion value
            'CLOUD_COVERAGE'                    : '-', 
            'OVERALL_QUALITY_FLAG'              : 'OK',
            'ASSOCIATED_QUALITY_INFORMATION'    : '-',
            'REGION_NAME'                       : None, # satellite area
            'COMPRESSION'                       : 0,
            'FIELD_TYPE'                        : 'Product',
            'FORECAST_STEP'                     : 0,
            'NC'                                : None, # nCols
            'NL'                                : None, # nLines
            'NB_PARAMETERS'                     : len(DATASETS),

            #'NOMINAL_PRODUCT_TIME'              : None, # product generation time
            #'SATELLITE'  # satelites que contribuiram para o produto (ver tipo de dados esperado)
            #'INSTRUMENT_ID' # nome dos sensores usados (ver tipo de dados)
            #'INSTRUMENT_MODE'
            #'IMAGE_ACQUISITION_TIME' #igual à timeslot do nome do ficheiro (primeira do período)
            #'ORBIT_TYPE'
            #'PROJECTION_NAME' #ver num dos produtos intermédios GEOS(-075.0)
            #'NOMINAL_LONG' # SSP
            #'NOMINAL_LAT' #0
            #'CFAC' # acrescentar aos settings do sistema
            #'LFAC' # acrescentar aos settings do sistema
            #'COFF' # acrescentar aos settings do sistema
            #'LOFF' # acrescentar aos settings do sistema
            #'START_ORBIT_NUMBER'
            #'END_ORBIT_NUMBER'
            #'SUB_SATELLITE_POINT_START_LAT'
            #'SUB_SATELLITE_POINT_START_LON'
            #'SUB_SATELLITE_POINT_END_LAT'
            #'SUB_SATELLITE_POINT_END_LON'
            #'SENSING_START_TIME' # igual à timeslot do nome do ficheiro (primeira do período)
            #'SENSING_END_TIME' # timeslot do fim do periodo
            #'PIXEL_SIZE' # acrescentar aos settings do sistema (4KM)
            #'GRANULE_TYPE'
            #'PROCESSING_LEVEL' # 03?
            #'PRODUCT_TYPE' #GEOGSA
            #'PRODUCT_ACTUAL_SIZE' #tamanho (em bytes) que os dados ocupam
            #'PROCESSING_MODE'
            #'DISPOSITION_FLAG'
            #'TIME_RANGE' #10-day
            #'STATISTIC_TYPE'
            #'MEAN_SSLAT'
            #'MEAN_SSLON'
            #'PLANNED_CHAN_PROCESSING'
            #'FIRST_LAT'
            #'FIRST_LON'
            }

    OUTPUT_PATTERN = r'g2_BIOPAR_GSA_#_#_GEO_v1'

    def __init__(self, filePath):
        fh = open(filePath, 'rb')
        self.natHeader = self._decode_header(fh)
        self._decode_data(fh)
        self.hdf5Attributes = self.__HDF5_ATTRIBUTES.copy()
        self._update_dynamic_attributes()
        fh.close()

    def _update_dynamic_attributes(self):
        '''Update the HDF5 attributes based on the input file's info.'''

        pass

    def _decode_header(self, fh):
        '''Decode the nat file's header information.'''

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
        cleanValues = self._clean_values(values)
        return dict(zip([f[0] for f in headerFields], cleanValues))

    def _clean_values(self, values):
        newValues = []
        for v in values:
            if isinstance(v, str):
                newV = re.sub('\x00', '', v).strip()
            else:
                newV = v
            newValues.append(newV)
        return newValues

    def _decode_data(self, fh):
        fh.seek(self.LINES_PART['start'])
        nLines, nCols = self._get_dimensions()
        lineDtype = np.dtype([
            ('header', 'i4', 6), # the line header is repeated twice
            ('data', 'u1', (nCols, len(self.DATASETS))),
            ])
        dataArr = np.fromfile(fh, dtype=lineDtype)
        for index, dsList in enumerate(self.DATASETS):
            self.DATASETS[index][1] = dataArr['data'][:,:,index]

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
        for k, v in self.natHeader.iteritems():
            exec('h5f.root._v_attrs.%s = "%s"' % (k, v))
        for k, v in self.hdf5Attributes.iteritems():
            exec('h5f.root._v_attrs.%s = "%s"' % (k, v))
        for dsList in self.DATASETS:
            name, arr, scalingFactor = dsList
            ds = h5f.createArray('/', name, arr)
            ds._v_attrs.MISSING_VALUE = 255
            ds._v_attrs.SCALING_FACTOR = scalingFactor
        h5f.close()
