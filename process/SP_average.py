import datetime
import glob
import os

import netCDF4
import numpy as np

from fileio.gems_l2_reader import GEMSL2Reader
from process.process import Process

from regrid.sp_merge import SpatialMerge

class SptAverageProcess(Process):
    path_config = None

    reader = GEMSL2Reader()

    def __init__(self, config: dict):
        super().__init__(config)
        self.path_config = self.config['environment_directory']
        self.grid_config = self.config['grid']

        os.makedirs(self.path_config['save_dir'], exist_ok=True)

    def _init_process(self):
        pass

    def run(self, dt: datetime.datetime, mode: str, alg_name: str):
        self.logger.info(f"Start Spatio-Temporal Average Process : {dt}")

        dt_filter = self.get_dt_filter(dt, mode)

        for region_name, grid_params in self.grid_config.items():
            filelist = glob.glob(os.path.join(self.path_config['intermediate_dir'], alg_name, "**", f"GK2_GEMS_L3_*{dt_filter}*_{alg_name}_GRID-{region_name.upper()}_hourly_{int(grid_params['dx'] * 100):03d}.nc"),
                                 recursive=True)

            if len(filelist) == 0:
                self.logger.warning(f"Files Not Exist!! : {dt} {dt_filter} {region_name} {grid_params['dx']}")
                return

            self.run_spt(filelist)
            #lon, lat, avg_datas = self.compute_average(filelist)

            #output_dir = os.path.join(self.path_config['save_dir'], mode, alg_name, dt_filter)
            #os.makedirs(output_dir, exist_ok=True)
            #output_filepath = os.path.join(output_dir,
            #                               f"GK2_GEMS_L3_{dt_filter}_{alg_name}_GRID-{region_name.upper()}_{mode}_{int(grid_params['dx'] * 100):03d}.nc")

            #self.write_product(output_filepath, lon, lat, avg_datas)

    def run_spt(self, filelist: list):
        self.logger.info("Start SPT Average")
        lon = None
        lat = None

        stack_datas = dict()

        for filepath in filelist:
            self.logger.debug(filepath)

            if lon is None or lat is None:
                geos = self.reader.read_geo(filepath)
                lon = geos['Longitude']
                lat = geos['Latitude']

            datas = self.reader.read_data(filepath)

            merger = SpatialMerge(lon=lon, lat=lat, datas=datas)
            datas = merger.execute()
            #for key, data in datas.items():
            #    if key not in stack_datas.keys():
            #        stack_datas[key] = []
            #    stack_datas[key].append(data)

        #avg_datas = dict()
        #for key, stack_data in stack_datas.items():
        #    avg_datas[key] = np.nanmean(np.array(stack_data), axis=0)

        self.logger.info("Done Compute Average")

        #return lon, lat, avg_datas

    def compute_average(self, filelist: list):
        self.logger.info("Start Compute Average")
        lon = None
        lat = None

        stack_datas = dict()

        for filepath in filelist:
            self.logger.debug(filepath)

            if lon is None or lat is None:
                geos = self.reader.read_geo(filepath)
                lon = geos['Longitude']
                lat = geos['Latitude']

            datas = self.reader.read_data(filepath)

            for key, data in datas.items():
                if key not in stack_datas.keys():
                    stack_datas[key] = []
                stack_datas[key].append(data)

        avg_datas = dict()
        for key, stack_data in stack_datas.items():
            avg_datas[key] = np.nanmean(np.array(stack_data), axis=0)

        self.logger.info("Done Compute Average")

        return lon, lat, avg_datas

    def write_product(self, output_filepath: str, lon: np.ndarray, lat: np.ndarray, avg_datas: dict):
        self.logger.info("Start Write GEMS L3 Product")

        image, spatial = lon.shape

        ds = netCDF4.Dataset(output_filepath, mode="w")

        ds.createDimension("image", image)
        ds.createDimension("spatial", spatial)

        geo_group = ds.createGroup("Geolocation Fields")
        variable = geo_group.createVariable("Longitude", lon.dtype.str, dimensions=("image", "spatial"))
        variable[:] = lon
        variable = geo_group.createVariable("Latitude", lat.dtype.str, dimensions=("image", "spatial"))
        variable[:] = lat

        data_group = ds.createGroup("Data Fields")

        for name, avg_data in avg_datas.items():
            variable = data_group.createVariable(name, avg_data.dtype.str, dimensions=("image", "spatial"))
            variable[:] = avg_data

        ds.close()
        ds = None

        self.logger.info("Done Write GEMS L3 Product")

    def get_dt_filter(self, dt, mode):
        if mode == "daily":
            dt_filter = dt.strftime("%Y%m%d")
        elif mode == "monthly":
            dt_filter = dt.strftime("%Y%m")
        elif mode == "yearly":
            dt_filter = dt.strftime("%Y")
        else:
            raise Exception(f"Not Supported Mode : {mode}")

        return dt_filter
