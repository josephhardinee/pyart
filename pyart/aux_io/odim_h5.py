"""
pyart.aux_io.odim_h5
====================

Routines for reading ODIM_H5 files.

.. autosummary::
    :toctree: generated/

    read_odim_h5

"""

import datetime

import numpy as np
import h5py

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, radar_coords_to_cart
from ..io.common import _test_arguments
from ..core.radar import Radar


ODIM_H5_FIELD_NAMES = {
    'TH': 'total_power',        # uncorrected reflectivity, horizontal
    'TV': 'total_power',        # uncorrected reflectivity, vertical
    'DBZH': 'reflectivity',     # corrected reflectivity, horizontal
    'DBZV': 'reflectivity',     # corrected reflectivity, vertical
    'ZDR': 'differential_reflectivity',     # differential reflectivity
    'RHOHV': 'cross_correlation_ratio',
    'LDR': 'linear_polarization_ratio',
    'PHIDP': 'differential_phase',
    'KDP': 'specific_differential_phase',
    'SQI': 'normalized_coherent_power',
    'SNR': 'signal_to_noise_ratio',
    'VRAD': 'velocity',
    'WRAD': 'spectrum_width',
    'QIND': 'quality_index',
}


def read_odim_h5(filename, field_names=None, additional_metadata=None,
                 file_field_names=False, exclude_fields=None, **kwargs):
    """
    Read a ODIM_H5 file.

    Parameters
    ----------
    filename : str
        Name of the ODIM_H5 file to read.
    field_names : dict, optional
        Dictionary mapping ODIM_H5 field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.


    Returns
    -------
    radar : Radar
        Radar object containing data from ODIM_H5 file.

    """
    # TODO before moving to pyart.io
    # * unit test
    # * add default field mapping, etc to default config
    # * auto-detect file type with pyart.io.read function
    # * instrument parameters
    # * add additional checks for HOW attributes
    # * support for other objects (SCAN, XSEC)

    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    if field_names is None:
        field_names = ODIM_H5_FIELD_NAMES
    filemetadata = FileMetadata('odim_h5', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # open the file
    hfile = h5py.File(filename, 'r')
    odim_object = hfile['what'].attrs['object']
    if odim_object != 'PVOL':
        raise NotImplementedError(
            'object: %s not implemented.' % (odim_object))

    # determine the number of sweeps by the number of groups which
    # begin with dataset
    datasets = [k for k in hfile if k.startswith('dataset')]
    datasets.sort()
    nsweeps = len(datasets)

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    h_where = hfile['where'].attrs
    latitude['data'] = np.array([h_where['lat']], dtype='float64')
    longitude['data'] = np.array([h_where['lon']], dtype='float64')
    altitude['data'] = np.array([h_where['height']], dtype='float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata['source'] = hfile['what'].attrs['source']
    metadata['original_container'] = 'odim_h5'
    metadata['odim_conventions'] = hfile.attrs['Conventions']

    h_what = hfile['what'].attrs
    metadata['version'] = h_what['version']
    metadata['source'] = h_what['source']

    try:
        ds1_how = hfile[datasets[0]]['how'].attrs
    except KeyError:
        # if no how group exists mock it with an empty dictionary
        ds1_how = {}
    if 'system' in ds1_how:
        metadata['system'] = ds1_how['system']
    if 'software' in ds1_how:
        metadata['software'] = ds1_how['software']
    if 'sw_version' in ds1_how:
        metadata['sw_version'] = ds1_how['sw_version']

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    rays_per_sweep = [hfile[d]['where'].attrs['nrays'] for d in datasets]
    total_rays = sum(rays_per_sweep)
    ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
    seri = np.cumsum(rays_per_sweep).astype('int32') - 1
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # sweep_number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    # scan_type
    scan_type = 'ppi'

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    sweep_el = [hfile[d]['where'].attrs['elangle'] for d in datasets]
    fixed_angle['data'] = np.array(sweep_el, dtype='float32')

    # elevation
    elevation = filemetadata('elevation')
    if 'elangles' in ds1_how:
        edata = np.empty(total_rays, dtype='float32')
        for d, start, stop in zip(datasets, ssri, seri):
            edata[start:stop+1] = hfile[d]['how'].attrs['elangles'][:]
        elevation['data'] = edata
    else:
        elevation['data'] = np.repeat(sweep_el, rays_per_sweep)

    # range
    _range = filemetadata('range')
    # check that the gate spacing is constant between sweeps
    rstart = [hfile[d]['where'].attrs['rstart'] for d in datasets]
    if any(rstart != rstart[0]):
        raise ValueError('range start changes between sweeps')
    rscale = [hfile[d]['where'].attrs['rscale'] for d in datasets]
    if any(rscale != rscale[0]):
        raise ValueError('range scale changes between sweeps')

    nbins = hfile['dataset1']['where'].attrs['nbins']
    _range['data'] = (np.arange(nbins, dtype='float32') * rscale[0] +
                      rstart[0] * 1000.)
    _range['meters_to_center_of_first_gate'] = rstart[0] * 1000.
    _range['meters_between_gates'] = float(rscale[0])

    # azimuth
    azimuth = filemetadata('azimuth')
    if ('startazA' in ds1_how) and ('stopazA' in ds1_how):
        # average between start and stop azimuth angles
        az_data = np.ones((total_rays, ), dtype='float32')
        for dset, start, stop in zip(datasets, ssri, seri):
            startaz = hfile[dset]['how'].attrs['startazA']
            stopaz = hfile[dset]['how'].attrs['stopazA']
            sweep_az = np.angle(
                (np.exp(1.j*np.deg2rad(startaz)) +
                 np.exp(1.j*np.deg2rad(stopaz))) / 2., deg=True)
            az_data[start:stop+1] = sweep_az
        azimuth['data'] = az_data
    else:
        # assume 1 degree per ray, starting at where/a1gate
        az_data = np.ones((total_rays, ), dtype='float32')
        for dset, start, stop in zip(datasets, ssri, seri):
            start_az = hfile[dset]['where'].attrs['a1gate']
            nrays = stop - start + 1
            sweep_az = np.fmod(start_az + np.arange(nrays), 360.)
            az_data[start:stop+1] = sweep_az
        azimuth['data'] = az_data

    # time
    _time = filemetadata('time')
    if ('startazT' in ds1_how) and ('stopazT' in ds1_how):
        # average between startazT and stopazT
        t_data = np.empty((total_rays, ), dtype='float32')
        for dset, start, stop in zip(datasets, ssri, seri):
            t_start = hfile[dset]['how'].attrs['startazT']
            t_stop = hfile[dset]['how'].attrs['stopazT']
            t_data[start:stop+1] = (t_start + t_stop) / 2
        start_epoch = t_data.min()
        start_time = datetime.datetime.utcfromtimestamp(start_epoch)
        _time['units'] = make_time_unit_str(start_time)
        _time['data'] = t_data - start_epoch
    else:
        t_data = np.empty((total_rays, ), dtype='int32')
        # interpolate between each sweep starting and ending time
        for dset, start, stop in zip(datasets, ssri, seri):
            dset_what = hfile[dset]['what'].attrs
            start_str = dset_what['startdate'] + dset_what['starttime']
            end_str = dset_what['enddate'] + dset_what['endtime']
            start_dt = datetime.datetime.strptime(start_str, '%Y%m%d%H%M%S')
            end_dt = datetime.datetime.strptime(end_str, '%Y%m%d%H%M%S')

            time_delta = end_dt - start_dt
            delta_seconds = time_delta.seconds + time_delta.days * 3600 * 24
            rays = stop - start + 1
            sweep_start_epoch = (
                start_dt - datetime.datetime(1970, 1, 1)).total_seconds()
            t_data[start:stop+1] = (sweep_start_epoch +
                                    np.linspace(0, delta_seconds, rays))
        start_epoch = t_data.min()
        start_time = datetime.datetime.utcfromtimestamp(start_epoch)
        _time['units'] = make_time_unit_str(start_time)
        _time['data'] = (t_data - start_epoch).astype('float32')

    # fields
    fields = {}
    h_field_keys = [k for k in hfile['dataset1'] if k.startswith('data')]
    odim_fields = [hfile['dataset1'][d]['what'].attrs['quantity'] for d in
                   h_field_keys]
    for odim_field, h_field_key in zip(odim_fields, h_field_keys):
        field_name = filemetadata.get_field_name(odim_field)
        if field_name is None:
            continue
        fdata = np.ma.zeros((total_rays, nbins), dtype='float32')
        start = 0
        # loop over the sweeps, copy data into correct location in data array
        for dset, rays_in_sweep in zip(datasets, rays_per_sweep):
            sweep_data = _get_odim_h5_sweep_data(hfile[dset][h_field_key])
            fdata[start:start + rays_in_sweep] = sweep_data[:]
            start += rays_in_sweep
        # create field dictionary
        field_dic = filemetadata(field_name)
        field_dic['data'] = fdata
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic

    # instrument_parameters
    instrument_parameters = None

    return Radar(
        _time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def _get_odim_h5_sweep_data(group):
    """ Get ODIM_H5 sweet data from an HDF5 group. """

    # mask raw data
    what = group['what']
    raw_data = group['data'][:]

    if 'nodata' in what.attrs:
        nodata = what.attrs.get('nodata')
        data = np.ma.masked_equal(raw_data, nodata)
    else:
        data = np.ma.masked_array(raw_data)
    if 'undetect' in what.attrs:
        undetect = what.attrs.get('undetect')
        data[data == undetect] = np.ma.masked

    offset = 0.0
    gain = 1.0
    if 'offset' in what.attrs:
        offset = what.attrs.get('offset')
    if 'gain' in what.attrs:
        gain = what.attrs.get('gain')
    return data * gain + offset
