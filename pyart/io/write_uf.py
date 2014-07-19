"""
pyart.io.write_uf
=================

Utilities for writing Universal Format(UF) files.

.. autosummary::
    :toctree: generated/

    write_uf
    _write_mandatory_header
    _write_ray
"""

import datetime
import platform
import struct

import numpy as np

from ..config import FileMetadata
from .common import stringarray_to_chararray
from ..core.radar import Radar


def write_uf(filename, radar):
    """
    Write a Radar object to a Universal Format(UF) file.

    The files produced by this routine follow the standard as laid out in the
    Vaisala Programmers Document.


    Parameters
    ----------
    filename : str
        Filename to create.
    radar : Radar
        Radar object.
    """

    num_rays = radar.nrays

    with open(filename, 'wb') as file_id:
        for ray in range(num_rays):
            print('Writing ray',ray)
            record_size = 2*_get_ray_size(radar)
            print('Record Size',record_size)
            print('Start of Ray',file_id.tell())
            file_id.write(struct.pack('>I', record_size))
            print('After Ray Size',file_id.tell())

            _write_mandatory_header(file_id, radar, ray)
            print('After mandatory_header:',file_id.tell())
            _write_optional_header(file_id, radar)
            print('After Optional Header:',file_id.tell())

            _write_data_header(file_id, radar, ray)
            print('After Data Header', file_id.tell())

            for field_num, field in enumerate(radar.fields.keys()):
                print('Before Field Header:',file_id.tell())
                offset = 50 + 27 * len(radar.fields.keys()) + radar.ngates * (len(radar.fields.keys())-1)
                _write_field_header(file_id, radar, field, ray, offset)
                print('After Field Header',file_id.tell())
                _write_data(file_id, radar, field, ray)
                print('After Data Write',file_id.tell())
            file_id.write(struct.pack('>I', record_size))



def _write_mandatory_header(file_id, radar, ray_num):
    # We build a dictionary to keep this slightly cleaner
    uf_mhpd = {'ID': 'UF'}

    uf_mhpd['RecordSize'] = _get_ray_size(radar)  # Fix this later
    uf_mhpd['OptionalHeaderPosition'] = 46
    uf_mhpd['LocalUseHeaderPosition'] = 60
    uf_mhpd['DataHeaderPosition'] = 60

    uf_mhpd['RecordNumber'] = 1
    uf_mhpd['VolumeNumber'] = 1
    uf_mhpd['RayNumber'] = ray_num + 1
    uf_mhpd['RecordInRay'] = 1
    uf_mhpd['SweepNumber'] = 1 + radar.nsweeps -\
        sum(radar.sweep_end_ray_index['data'] > ray_num)

    uf_mhpd['RadarName'] = radar.metadata['instrument_name'][0:8]
    uf_mhpd['SiteName'] = _prep_string('UNSET', 8)  # Fix this later
    uf_mhpd['LatDegrees'] = 25  # Fix this later
    uf_mhpd['LatMinutes'] = 25  # Fix this later
    uf_mhpd['LatSeconds'] = 59  # Fix this later
    uf_mhpd['LonDegrees'] = 25  # Fix this later
    uf_mhpd['LonMinutes'] = 25  # Fix this later
    uf_mhpd['LonSeconds'] = 25  # Fix this later
    uf_mhpd['Altitude'] = 25  # Fix this later
    uf_mhpd['Year'] = 25  # Fix this later
    uf_mhpd['Month'] = 8  # Fix this later
    uf_mhpd['Day'] = 15  # Fix this later
    uf_mhpd['Hour'] = 21  # Fix this later
    uf_mhpd['Minute'] = 35  # Fix this later
    uf_mhpd['Second'] = 45  # Fix this later
    uf_mhpd['TimeZone'] = 'UT'  # Fix this later
    uf_mhpd['Azimuth'] = 35  # Fix this later
    uf_mhpd['Elevation'] = 12  # Fix this later
    uf_mhpd['SweepMode'] = 1  # Fix this later
    uf_mhpd['FixedAngle'] = 25  # Fix this later
    uf_mhpd['SweepRate'] = 25  # Fix this later
    uf_mhpd['ConvertYear'] = 25  # Fix this later
    uf_mhpd['ConvertMonth'] = 8  # Fix this later
    uf_mhpd['ConvertDay'] = 25  # Fix this later
    uf_mhpd['ConvertName'] = 'PYART   '  # Fix this later
    uf_mhpd['MissingDataValue'] = 25  # Fix this later

    format_string = '>' + ''.join([rec_type[1]
                                  for rec_type in UF_MANDATORY_HEADER2])
    data_string = [uf_mhpd[rec_type[0]] for rec_type in UF_MANDATORY_HEADER2]
    file_id.write(struct.pack(format_string, *data_string))


def _write_data_header(file_id, radar, ray):
    '''
    Write the uf_data_header2 structure
    '''

    fh_pos1 = 62+ 2 * len(radar.fields.keys()) + 1
    fhk = 25 + radar.ngates

    uf_dhd = {}
    uf_dhd['FieldsThisRay'] = len(radar.fields)
    uf_dhd['RecordsThisRay'] = 1
    uf_dhd['FieldsThisRecord'] = len(radar.fields)

    format_string = '>' + ''.join([rec_type[1]
                                  for rec_type in uf_data_header2])
    data_string = [uf_dhd[rec_type[0]] for rec_type in uf_data_header2]
    file_id.write(struct.pack(format_string, *data_string))

    uf_fhl = []
    for idx, field in enumerate(radar.fields.keys()):
        # Write Field Header Position, and Field Data Type
        file_id.write(struct.pack('>2sh', field[0:2], fh_pos1 + fhk * idx))


def _write_field_header(file_id, radar, field, ray_num, offset):
    '''
    Generates and writes uf_field_header2 structure to file.
    '''
    print('offset:%d' %offset)
    uf_fhd = {}
    uf_fhd['DataPosition'] = offset
    uf_fhd['ScaleFactor'] = 1
    uf_fhd['StartRangeKm'] = 0
    uf_fhd['StartRangeMeters'] = 0
    uf_fhd['BinSpacing'] = radar.range['meters_between_gates']
    uf_fhd['BinCount'] = radar.ngates
    uf_fhd['PulseWidth'] = 0
    uf_fhd['BeamWidthH'] = 1
    uf_fhd['BeamWidthV'] = 1
    uf_fhd['BandWidth'] = 0
    uf_fhd['Polarization'] = 0
    uf_fhd['Wavelength'] = 0
    uf_fhd['SampleSize'] = 0
    uf_fhd['ThresholdData'] = 'NC'
    uf_fhd['ThresholdValue'] = 0
    uf_fhd['Scale'] = 0
    uf_fhd['EditCode'] = 'UF'
    uf_fhd['PRT'] = 0
    uf_fhd['BitsPerBin'] = 16

    format_string = '>' + ''.join([rec_type[1]
                                  for rec_type in uf_field_header2])
    data_string = [uf_fhd[rec_type[0]] for rec_type in uf_field_header2]
    file_id.write(struct.pack(format_string, *data_string))
    file_id.write(struct.pack('>12s','NOTSETYET!!'))


def _write_data(file_id, radar, field, ray_num):
    format_string = '>' + str(radar.ngates) + 'h'
    file_id.write(struct.pack(
        format_string, *list(radar.fields[field]['data'][ray_num].data)))

def _write_optional_header(file_id, radar):
    uf_ohd = {}
    uf_ohd['ProjectName'] = 'PYART   '  # Fix this later
    uf_ohd['BaselineAzimuth']=2
    uf_ohd['BaselineElevation']= 3
    uf_ohd['VolumeScanHour'] = 4
    uf_ohd['VolumeScanMinute'] = 5
    uf_ohd['VolumeScanSecond']= 6
    uf_ohd['FieldTapeName']='FIELDTAP'
    uf_ohd['iFlag']=7

    format_string = '>' + ''.join([rec_type[1]
                                  for rec_type in uf_optional_header])
    data_string = [uf_ohd[rec_type[0]] for rec_type in uf_optional_header]
    file_id.write(struct.pack(format_string, *data_string))


def _prep_string(string, length):
    '''
    This just pads strings out and null terminates them. I'm probably
    missing an obvious python function to do this.
    '''
    length_of_string = len(string)
    if(length_of_string >= length):
        return ''.join((string[0:length], '\n'))
    else:
        return ''.join((string, '\n', ' ' * (length - length_of_string - 1)))


def _get_ray_size(radar):
    ''' Calculate the ray size for each ray in the UF file

    '''
    #ray_size_const =50+14
    #dh_size_per_field = 3 + 2 * len(radar.fields.keys())
    #fh_size = 25 * len(radar.fields.keys())
    #data_size = radar.ngates * len(radar.fields.keys())
    #record_size = ray_size_const + dh_size_per_field + fh_size + data_size
    record_size = 52 + (27 + radar.ngates)*len(radar.fields.keys())
    return record_size


UF_MANDATORY_HEADER2 = (
    ('ID', '2s'),
    ('RecordSize', 'h'),
    ('OptionalHeaderPosition', 'h'),
    ('LocalUseHeaderPosition', 'h'),
    ('DataHeaderPosition', 'h'),
    ('RecordNumber', 'h'),
    ('VolumeNumber', 'h'),
    ('RayNumber', 'h'),
    ('RecordInRay', 'h'),
    ('SweepNumber', 'h'),
    ('RadarName', '8s'),
    ('SiteName', '8s'),
    ('LatDegrees', 'h'),
    ('LatMinutes', 'h'),
    ('LatSeconds', 'h'),
    ('LonDegrees', 'h'),
    ('LonMinutes', 'h'),
    ('LonSeconds', 'h'),
    ('Altitude', 'h'),
    ('Year', 'h'),
    ('Month', 'h'),
    ('Day', 'h'),
    ('Hour', 'h'),
    ('Minute', 'h'),
    ('Second', 'h'),
    ('TimeZone', '2s'),
    ('Azimuth', 'h'),
    ('Elevation', 'h'),
    ('SweepMode', 'h'),
    ('FixedAngle', 'h'),
    ('SweepRate', 'h'),
    ('ConvertYear', 'h'),
    ('ConvertMonth', 'h'),
    ('ConvertDay', 'h'),
    ('ConvertName', '8s'),
    ('MissingDataValue', 'h')
)

uf_optional_header = (
    ('ProjectName', '8s'),
    ('BaselineAzimuth', 'h'),
    ('BaselineElevation', 'h'),
    ('VolumeScanHour', 'h'),
    ('VolumeScanMinute', 'h'),
    ('VolumeScanSecond', 'h'),
    ('FieldTapeName', '8s'),
    ('iFlag', 'h')
)


uf_dsi2 = (
    ('DataType', '2s'),
    ('FieldHeaderPosition', 'h')
)


uf_data_header2 = (
    ('FieldsThisRay', 'h'),
    ('RecordsThisRay', 'h'),
    ('FieldsThisRecord', 'h')
)

uf_field_header2 = (
    ('DataPosition', 'h'),
    ('ScaleFactor', 'h'),
    ('StartRangeKm', 'h'),
    ('StartRangeMeters', 'h'),
    ('BinSpacing', 'h'),
    ('BinCount', 'h'),
    ('PulseWidth', 'h'),
    ('BeamWidthH', 'h'),
    ('BeamWidthV', 'h'),
    ('BandWidth', 'h'),
    ('Polarization', 'h'),
    ('Wavelength', 'h'),
    ('SampleSize', 'h'),
    ('ThresholdData', '2s'),
    ('ThresholdValue', 'h'),
    ('Scale', 'h'),
    ('EditCode', '2s'),
    ('PRT', 'h'),
    ('BitsPerBin', 'h')
)

uf_field_dm = (
    ('RadarConstant', 'h'),
    ('NoisePower', 'h'),
    ('ReceiverGain', 'h'),
    ('PeakPower', 'h'),
    ('AntennaGain', 'h'),
    ('PulseDuration', 'h')
)

uf_field_ve = (
    ('NyquistVelocity', 'I')
)
