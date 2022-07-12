import netCDF4
import numpy as np


class GEMSL2_Nc_Reader:

    @staticmethod
    def read(filepath, variables: list = None):
        ds = netCDF4.Dataset(filepath)
        # ds.set_auto_maskandscale(False)

        datas = GEMSL2_Nc_Reader._read_datafield(ds['Data Fields'], variables)
        geo_datas = GEMSL2_Nc_Reader._read_geo(ds['Geolocation Fields'])

        ds.close()
        ds = None

        return geo_datas, datas

    @staticmethod
    def _read_geo(group: netCDF4.Group):
        variables = [
            'Longitude',
            'Latitude',
            'SolarZenithAngle',
            'ViewingZenithAngle'
        ]

        data = dict()
        for variable in variables:
            data[variable] = group[variable][:].filled(np.nan)

        return data

    @staticmethod
    def _read_datafield(group: netCDF4.Group, variables: list = None):
        if variables is None:
            variables = group.variables

        if "CloudFraction" not in variables:
            variables.append("CloudFraction")

        data = dict()
        for variable in variables:
            try:
                data[variable] = group[variable][:].filled(np.nan)
            except:
                data[variable] = group[variable][:].data

        return data
