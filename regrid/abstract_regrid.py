from abc import ABCMeta, abstractmethod

import numpy as np


class AbstractRegrid(metaclass=ABCMeta):
    lon = None
    lat = None
    nrows = 0
    ncols = 0
    datas = None

    def __init__(self, lon: np.ndarray, lat: np.ndarray, datas: dict):
        self.lon = lon.copy()
        self.lat = lat.copy()
        self.datas = datas.copy()

        self.nrows = self.lon.shape[0]
        self.ncols = self.lon.shape[1]
        self._init_grid()

    def __call__(self, new_lon, new_lat):
        self.execute(new_lon, new_lat)

    @abstractmethod
    def _init_grid(self):
        pass

    @abstractmethod
    def execute(self, new_lon, new_lat):
        pass
