import datetime
import glob
import os.path
import re

import netCDF4
import numpy as np
from matplotlib import pyplot as plt

from common.constant import GEMS_L2_FILENAME_PATTERN
from fileio.fileio import GEMSL2_Nc_Reader
from fileio.gems_l2_reader import GEMSL2Reader
from process.process import Process
from regrid.inversedistanceweight import InverseDistanceWeight
from regrid.tessellation import Tessellation


class L2convertProcess(Process):
    path_config = None
    grid_param_config = None
    alg_param_config = None
    input_param_config = None

    reader = GEMSL2Reader()

    # Path Variables
    target_dir = None
    intermediate_dir = None

    # latlon range Variables
    lon_min = None
    lon_max = None
    lat_min = None
    lat_max = None
    grid_size = None

    new_lon = None
    new_lat = None

    # Alg Variabls
    cf_thresh = None
    sza_thresh = None
    vza_thresh = None
    regrid_alg = None
    fill_value = None

    def __init__(self, config: dict):
        super().__init__(config)
        self.path_config = self.config['environment_directory']
        self.grid_param_config = self.config['grid_param']
        self.alg_param_config = self.config['alg_param']
        self.input_param_config = self.config['input_param']

    def _init_process(self) -> None:
        self.logger.info("Init Process")
        self.init_variables()

    def init_variables(self):
        self.init_path_variables()
        self.init_grid_params()
        self.init_alg_params()

    def init_path_variables(self):
        self.target_dir = self.path_config['target_dir']
        self.intermediate_dir = self.path_config['intermediate_dir']

        os.makedirs(self.intermediate_dir, exist_ok=True)

    def init_grid_params(self):
        self.lon_min = self.grid_param_config['lon_min']
        self.lon_max = self.grid_param_config['lon_max']
        self.lat_min = self.grid_param_config['lat_min']
        self.lat_max = self.grid_param_config['lat_max']
        self.grid_size = self.grid_param_config['grid_size']

        self.new_lon = np.arange(self.lon_min, self.lon_max, self.grid_size)
        self.new_lat = np.arange(self.lat_max, self.lat_min, -self.grid_size)

        self.new_lon, self.new_lat = np.meshgrid(self.new_lon, self.new_lat)

    def init_alg_params(self):
        self.cf_thresh = self.alg_param_config['cf_thresh']
        self.sza_thresh = self.alg_param_config['sza_thresh']
        self.vza_thresh = self.alg_param_config['vza_thresh']
        self.regrid_alg = self.alg_param_config['regrid_alg']

    def run(self, dt: datetime.datetime) -> None:
        self.logger.debug(f"Run : {dt}")
        self._init_process()

        for alg_name, input_params in self.input_param_config.items():
            base_filename = f"GK2_GEMS_L2_{dt.strftime('%Y%m%d_%H')}45_{alg_name}_*_{input_params['binning']}.nc"

            input_file_regex = os.path.join(self.target_dir, alg_name, dt.strftime('%Y%m'), dt.strftime('%d'),
                                            base_filename)

            input_filepath = glob.glob(input_file_regex)
            print(input_file_regex)
            if len(input_filepath) < 1:
                raise FileNotFoundError(input_filepath)
            elif len(input_filepath) > 1:
                raise Exception(f"Input File is so Many!! : {input_filepath}")

            input_filepath = input_filepath[0]
            input_filename = os.path.basename(input_filepath)

            matcher = re.match(GEMS_L2_FILENAME_PATTERN, input_filename)
            if matcher:
                area = matcher.group("area")

                self.logger.info(f"Read Input File : {input_filepath}")
                geo_datas, datas = GEMSL2_Nc_Reader.read(input_filepath, input_params['variables'])

                cf = datas['CloudFraction']
                del datas['CloudFraction']

                lon = geo_datas['Longitude']
                lat = geo_datas['Latitude']
                sza = geo_datas['SolarZenithAngle']
                vza = geo_datas['ViewingZenithAngle']

                nan_mask = (cf >= self.cf_thresh) | (sza >= self.sza_thresh) | (vza >= self.vza_thresh)

                for alg in datas.keys():
                    datas[alg][nan_mask] = np.nan

                # Auto-detect Half/Full PATTERNS
                if area.startswith("H"):
                    self.logger.warning("Half Product Process. Maybe Not Need Soon... (SOL On)")
                    lat = (lat[:, 0::2] + lat[:, 1::2]) / 2
                    lon = (lon[:, 0::2] + lon[:, 1::2]) / 2
                    for alg in datas.keys():
                        data = datas[alg]
                        datas[alg] = (data[:, 0::2] + data[:, 1::2]) / 2

                # Process Tessellation
                self.logger.info(f'Regrid Algorithm : {self.regrid_alg}')
                st = datetime.datetime.now()
                if self.regrid_alg == "tessellation":
                    worker = Tessellation(lon, lat, datas)
                elif self.regrid_alg == "idw":
                    worker = InverseDistanceWeight(lon, lat, datas)
                else:
                    raise Exception(f"This Algorithm Not Supported!! : {self.regrid_alg}")

                result = worker.execute(self.new_lon, self.new_lat)
                self.logger.debug(f"Regrid Process Elapsed Time : {(datetime.datetime.now() - st).total_seconds()} s")

                # fig, axs = plt.subplots(1, 2, figsize=[16, 9], sharex=True, sharey=True)
                # axs[0].pcolormesh(lon, lat, datas['ColumnAmountNO2'], cmap='jet')
                # axs[1].pcolormesh(self.new_lon, self.new_lat, result['ColumnAmountNO2'])
                # plt.show()

                output_dir = os.path.join(self.intermediate_dir, alg_name, dt.strftime("%Y%m"), dt.strftime("%d"))
                os.makedirs(output_dir, exist_ok=True)

                output_filepath = os.path.join(output_dir,
                                               f"GK2_GEMS_L3_{dt.strftime('%Y%m%d_%H')}45_{alg_name}_WGS84.nc")
                self.logger.info(f'Write GEMS L3 File : {output_filepath}')
                self.write_product(output_filepath, result)
            else:
                raise FileNotFoundError(f"This File is Not GEMS L2 Product : {input_filepath}")

    def write_product(self, output_filepath, result):
        image, spartial = self.new_lon.shape

        ds = netCDF4.Dataset(output_filepath, 'w')

        ds.createDimension('image', image)
        ds.createDimension('spatial', spartial)

        geo_group = ds.createGroup('Geolocation Fields')

        var_lon = geo_group.createVariable('Longitude', self.new_lon.dtype.str, dimensions=("image", "spatial"),
                                           zlib=True)
        var_lon[:] = self.new_lon
        var_lat = geo_group.createVariable('Latitude', self.new_lat.dtype.str, dimensions=("image", "spatial"),
                                           zlib=True)
        var_lat[:] = self.new_lat

        data_group = ds.createGroup('Data Fields')
        for key, data in result.items():
            var_data = data_group.createVariable(key, data.dtype.str, dimensions=("image", "spatial"), zlib=True)
            var_data[:] = data

        ds.close()
        ds = None
