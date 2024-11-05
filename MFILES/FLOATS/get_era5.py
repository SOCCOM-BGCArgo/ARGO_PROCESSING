#!/usr/bin/env python
import cdsapi
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': 'surface_pressure',
        'year': '2024',
        'month': ['1','2','3','4','5','6','7','8'],
        'day': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'],
        'time': ['00:00','06:00','12:00','18:00'],
        'format': 'netcdf',
    },
    '//atlas/chem/argo_processing/data/era5/sp/ERA5_sp_2024.nc')