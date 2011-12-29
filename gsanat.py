#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
This module provides a class to convert outputs from EUMETSAT's Geostationary 
Surface Albedo (GSA) algorithm from their native binary format to HDF5.
"""

import sys
from struct import unpack
import numpy as np
import tables
import re
import datetime as dt

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

    def __init__(self, filePath):
        self.params = self._extract_params(filePath)
        fh = open(filePath, 'rb')
        self.natHeader = self._decode_header(fh)
        self._decode_data(fh)
        fh.close()
        self.hdf5Attributes = self._get_hdf5_attrs()

    def _extract_params(self, theString):
        '''Extract relevant file parameters from its filename.'''

        patt = re.compile(r'(?P<source>[a-zA-Z]{4})_(?P<ssp>\d{3})_.*' \
                          '_(?P<year>\d{4})_(?P<firstDoy>\d{3})' \
                          '_(?P<lastDoy>\d{3})')
        match = patt.search(theString)
        params = None
        if match is not None:
            items = match.groupdict()
            for k, v in items.iteritems():
                if k in ('ssp', 'year', 'firstDoy', 'lastDoy'):
                    items[k] = int(v)
            firstTs = dt.datetime(items['year'], 1, 1) + \
                      dt.timedelta(days=items['firstDoy'] - 1)
            lastTs = firstTs + dt.timedelta(days=items['lastDoy'] - 1)
            params = {
                    'firstDay' : firstTs, 'lastDay' : lastTs,
                    'source' : items['source'], 'ssp' : items['ssp']}
        return params

    def _get_hdf5_attrs(self):

        nLines, nCols = self._get_dimensions()
        now = dt.datetime.utcnow()
        hdf5Attrs = {
            'SAF'                               : 'GEOLAND2',
            'CENTRE'                            : 'IM-PT',
            'ARCHIVE_FACILITY'                  : 'IM-PT',
            'PRODUCT'                           : 'GSA',
            'PARENT_PRODUCT_NAME'               : ['-', '-', '-', '-'],
            'SPECTRAL_CHANNEL_ID'               : 0,
            'PRODUCT_ALGORITHM_VERSION'         : self.natHeader['szMSAVersion'],
            'CLOUD_COVERAGE'                    : '-', 
            'OVERALL_QUALITY_FLAG'              : 'OK',
            'ASSOCIATED_QUALITY_INFORMATION'    : '-',
            'REGION_NAME'                       : None, # satellite area
            'COMPRESSION'                       : 0,
            'FIELD_TYPE'                        : 'Product',
            'FORECAST_STEP'                     : 0,
            'NC'                                : nCols,
            'NL'                                : nLines,
            'NB_PARAMETERS'                     : len(self.DATASETS),
            'NOMINAL_PRODUCT_TIME'              : now.strftime('%Y%m%d%H%M'),
            'SATELLITE'                         : None, # satelites used (list), get from g2System settings
            'INSTRUMENT_ID'                     : None, # sensor name, get from g2System settings 
            'INSTRUMENT_MODE'                   : 'STATIC_VIEW',
            'IMAGE_ACQUISITION_TIME'            : self.params['firstDay'],
            'ORBIT_TYPE'                        : 'GEO',
            'PROJECTION_NAME'                   : None, # 'GEOS(-075.0)', get from g2System settings
            'NOMINAL_LONG'                      : None, # SSP, get from g2System settings
            'NOMINAL_LAT'                       : 0.0,
            'CFAC'                              : None, # get from g2System settings
            'LFAC'                              : None, # get from g2System settings
            'COFF'                              : None, # get from g2System settings
            'LOFF'                              : None, # get from g2System settings
            'START_ORBIT_NUMBER'                : 0,
            'END_ORBIT_NUMBER'                  : 0,
            'SUB_SATELLITE_POINT_START_LAT'     : 0.0,
            'SUB_SATELLITE_POINT_START_LON'     : 0.0,
            'SUB_SATELLITE_POINT_END_LAT'       : 0.0,
            'SUB_SATELLITE_POINT_END_LON'       : 0.0,
            'SENSING_START_TIME'                : self.params['firstDay'],
            'SENSING_END_TIME'                  : self.params['lastDay'],
            'PIXEL_SIZE'                        : None, # get from g2System settings (4KM)
            'GRANULE_TYPE'                      : 'DP',
            'PROCESSING_LEVEL'                  : '03', # is this correct?
            'PRODUCT_TYPE'                      : 'GEOGSA',
            'PRODUCT_ACTUAL_SIZE'               : nCols * nLines * len(self.DATASETS),
            'PROCESSING_MODE'                   : 'N', # it would be nice to have this as a dynamic attribute
            'DISPOSITION_FLAG'                  : 'O', # it would be nice to have this as a dynamic attribute
            'TIME_RANGE'                        : '10-day',
            'STATISTIC_TYPE'                    : '-',
            'MEAN_SSLAT'                        : 0.0,
            'MEAN_SSLON'                        : 0.0,
            'PLANNED_CHAN_PROCESSING'           : 0,
            'FIRST_LAT'                         : 0.0,
            'FIRST_LON'                         : 0.0,
            }
        return hdf5Attrs

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
