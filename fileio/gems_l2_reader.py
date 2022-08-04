#-*- encoding: utf8 -*-
import netCDF4
import numpy as np


class GEMSL2Reader:
    def read_data(self, filepath):
        ds = netCDF4.Dataset(filepath)

        result = self._read_impl(ds['Data Fields'])

        ds.close()
        ds = None

        return result

    def read_geo(self, filepath):
        ds = netCDF4.Dataset(filepath)

        result = self._read_impl(ds['Geolocation Fields'])

        ds.close()
        ds = None

        return result

    def _read_impl(self, group, result=None):
        if result is None:
            result = dict()

        for name, variable in group.variables.items():
            try:
                result[name] = variable[:].filled(np.nan)
            except:
                result[name] = variable[:].data

        for name, group in group.groups.items():
            self._read_impl(group, result)

        return result
