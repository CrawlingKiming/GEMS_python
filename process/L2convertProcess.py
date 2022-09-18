import datetime
import glob
import os.path
import re

import netCDF4
import numpy as np
from fontTools.misc.intTools import bit_count

from common.constant import GEMS_L2_FILENAME_PATTERN
from fileio.fileio import GEMSL2_Nc_Reader
from process.process import Process
from regrid.inverse_distance_weight import InverseDistanceWeight
from regrid.tessellation import Tessellation
from util.grid_utils import make_grid, half_scene_avg


class L2convertProcess(Process):
    path_config = None
    grid_config = None
    alg_param_config = None
    input_param_config = None

    # Path Variables
    target_dir = None
    intermediate_dir = None

    # latlon range Variables
    grid_lon = None
    grid_lat = None
    grid_size = None

    # Alg Variabls
    cf_thresh = None
    sza_thresh = None
    vza_thresh = None
    regrid_alg = None
    fill_value = None
    qf_usage = None
    p = 2

    def __init__(self, config: dict):
        super().__init__(config)
        self.path_config = self.config['environment_directory']
        self.grid_config = self.config['grid']
        self.alg_param_config = self.config['alg_param']
        self.input_param_config = self.config['input_param']

    def _init_process(self) -> None:
        self.logger.info("Init Process")
        self.init_variables()

    def init_variables(self):
        self.init_path_variables()
        self.init_alg_params()

    def init_path_variables(self):
        self.target_dir = self.path_config['target_dir']
        self.intermediate_dir = self.path_config['intermediate_dir']
        os.makedirs(self.intermediate_dir, exist_ok=True)

    def init_alg_params(self):
        self.cf_thresh = self.alg_param_config['cf_thresh']
        self.sza_thresh = self.alg_param_config['sza_thresh']
        self.vza_thresh = self.alg_param_config['vza_thresh']
        self.regrid_alg = self.alg_param_config['regrid_alg']
        self.qf_usage = self.alg_param_config['qf_usage']
        self.p = self.alg_param_config['p']

    def init_grid_params(self, grid_params):
        self.grid_lon, self.grid_lat = make_grid(**grid_params)
        self.grid_size = abs(grid_params['dx'])

        """
        # TODO : assumptions?
        self.lon1 = np.arange(self.lon_min, self.sublon_min, self.grid_size)
        self.lon2 = np.arange(self.sublon_min, self.sublon_max, self.subgrid_size)
        self.lon3 = np.arange(self.sublon_max, self.lon_max, self.grid_size)

        self.lat1 = np.arange(self.lat_max, self.sublat_max, -self.grid_size)
        self.lat2 = np.arange(self.sublat_max, self.sublat_min, -self.subgrid_size)
        self.lat3 = np.arange(self.sublat_min, self.lat_min, -self.grid_size)

        self.pc1 = np.meshgrid(self.lon1, self.lat1); self.pc2 = np.meshgrid(self.lon2, self.lat1); self.pc3 = np.meshgrid(self.lon3, self.lat1)
        self.pc4 = np.meshgrid(self.lon1, self.lat2); self.pc5 = np.meshgrid(self.lon2, self.lat2); self.pc6 = np.meshgrid(self.lon3, self.lat2)
        self.pc7 = np.meshgrid(self.lon1, self.lat3); self.pc8 = np.meshgrid(self.lon2, self.lat3); self.pc9 = np.meshgrid(self.lon3, self.lat3)

        self.lon_row1 = np.hstack((self.pc1[0], self.pc2[0], self.pc3[0]))
        # self.lon_row1 = np.round(np.hstack((self.pc1[0], self.pc2[0], self.pc3[0])),1)
        self.lon_row2 = np.hstack((self.pc4[0], self.pc5[0], self.pc6[0]))
        self.lon_row3 = np.hstack((self.pc7[0], self.pc8[0], self.pc9[0]))
        # self.lon_row3 = np.round(np.hstack((self.pc7[0], self.pc8[0], self.pc9[0])),1)

        self.lat_row1 = np.hstack((self.pc1[1], self.pc2[1], self.pc3[1]))
        # self.lat_row2 = np.hstack((np.round(self.pc4[1],1), self.pc5[1], np.round(self.pc6[1],1)))
        self.lat_row2 = np.hstack((self.pc4[1], self.pc5[1], self.pc6[1]))
        self.lat_row3 = np.hstack((self.pc7[1], self.pc8[1], self.pc9[1]))

        self.new_lon = np.vstack((self.lon_row1, self.lon_row2, self.lon_row3))
        self.new_lat = np.vstack((self.lat_row1, self.lat_row2, self.lat_row3))
        """

    def qf_algo(self, datas):
        """

        Calculates qf weights

        :param datas:
        :return: qf matrix
        """
        #############
        #           #
        # Former Code
        #           #
        #############
        """
        binary_list = [k for k in range(2 ** 12)]
        binary_list = [np.binary_repr(k, width=12) for k in binary_list]
        algorithm_binary_list = [int(k[11 - 1]) + int(k[11 - 3]) + int(k[11 - 9]) for k in binary_list]
        amf_binary_list = [int(k[11 - 1]) for k in binary_list]

        qf1 = datas['AlgorithmQualityFlags']  # 1 3 9인 QF만 사용
        for i in range(qf1.shape[0]):
            for j in range(qf1.shape[1]):
                qf1[i, j] = algorithm_binary_list[qf1[i, j]]

        del datas['AlgorithmQualityFlags']

        qf2 = datas['AMFQualityFlags']  # 1인 QF만 사용
        for i in range(qf2.shape[0]):
            for j in range(qf2.shape[1]):
                qf2[i, j] = amf_binary_list[qf2[i, j]]

        del datas['AMFQualityFlags']

        qf = qf1 + qf2
        for i in range(qf.shape[0]):
            for j in range(qf.shape[1]):
                if qf[i, j] == 0:  # Weight=1 if All QF positive
                    qf[i, j] = 1

                else:
                    qf[i, j] = 2  # Weight=0.5 if even one negative

        return qf
        """

        to_bitcount = np.vectorize(bit_count)

        qf1 = datas['AlgorithmQualityFlags'].astype(int)
        qf1 = to_bitcount(qf1 & 0b1000001010)  # 1 3 9인 QF만 사용

        del datas['AlgorithmQualityFlags']

        qf2 = datas['AMFQualityFlags'].astype(int)
        qf2 = to_bitcount(qf2 & 0b0000000010)  # 1인 QF만 사용

        del datas['AMFQualityFlags']

        qf = np.where((qf1 + qf2) == 0, 1, 2)


        return qf

    def run(self, dt: datetime.datetime) -> None:
        self.logger.debug(f"Run : {dt}")
        self._init_process()

        for alg_name, input_params in self.input_param_config.items():
            base_filename = f"GK2_GEMS_L2_{dt.strftime('%Y%m%d_%H')}45_{alg_name}_*_{input_params['binning']}.nc"

            input_file_regex = os.path.join(self.target_dir, alg_name, dt.strftime('%Y%m'), dt.strftime('%d'),
                                            base_filename)

            input_filepath = glob.glob(input_file_regex)

            if len(input_filepath) < 1:
                raise FileNotFoundError(input_filepath)
            elif len(input_filepath) > 1:
                raise Exception(f"Input File is so Many!! : {input_filepath}")

            input_filepath = input_filepath[0]
            input_filename = os.path.basename(input_filepath)

            matcher = re.match(GEMS_L2_FILENAME_PATTERN, input_filename)
            if matcher:
                self.logger.info(f"Read Input File : {input_filepath}")
                geo_datas, datas = GEMSL2_Nc_Reader.read(input_filepath, input_params['variables'])

                # Flip LR
                for key in geo_datas.keys():
                    geo_datas[key] = np.fliplr(geo_datas[key])

                for key in datas.keys():
                    datas[key] = np.fliplr(datas[key])

                cf = datas['CloudFraction']
                del datas['CloudFraction']

                lon = geo_datas['Longitude']
                lat = geo_datas['Latitude']
                sza = geo_datas['SolarZenithAngle']
                vza = geo_datas['ViewingZenithAngle']

                nan_mask = (cf >= self.cf_thresh) | (sza >= self.sza_thresh) | (vza >= self.vza_thresh)

                for alg in datas.keys():
                    if "float" in datas[alg].dtype.name:
                        datas[alg][nan_mask] = np.nan

                # Auto-detect Half/Full PATTERNS
                if matcher.group("area").startswith("H") and lon.shape[1] > 450:
                    self.logger.warning("Half Product Process. Maybe Not Need Soon... (SOL On)")
                    lat = half_scene_avg(lat)
                    lon = half_scene_avg(lon)
                    for alg in datas.keys():
                        datas[alg] = half_scene_avg(datas[alg])

                # QF works
                if self.qf_usage:
                    self.logger.info(f"Use Quality Flag")
                    qf = self.qf_algo(datas)
                else:
                    self.logger.info(f"Not Use Quality Flag")
                    qf = np.ones(lon.shape)

                st = datetime.datetime.now()
                if self.regrid_alg == "tessellation":
                    worker = Tessellation(lon, lat, datas, qf)
                elif self.regrid_alg == "idw":
                    worker = InverseDistanceWeight(lon, lat, datas, qf, self.alg_param_config['radius'],
                                                   self.alg_param_config['p'])
                else:
                    raise Exception(f"This Algorithm Not Supported!! : {self.regrid_alg}")

                # Process Regrid
                for region_name, grid_params in self.grid_config.items():
                    self.logger.info(f"Regrid Area : {region_name}")
                    self.init_grid_params(grid_params)

                    self.logger.info(f'Regrid Algorithm : {self.regrid_alg}')
                    result = worker.execute(self.grid_lon, self.grid_lat)
                    self.logger.debug(
                        f"Regrid Process Elapsed Time : {(datetime.datetime.now() - st).total_seconds()} s")

                    output_dir = os.path.join(self.intermediate_dir, alg_name, dt.strftime("%Y%m"), dt.strftime("%d"))
                    os.makedirs(output_dir, exist_ok=True)

                    output_filename = self.get_output_filename(input_filename, matcher, region_name)
                    output_filepath = os.path.join(output_dir, output_filename)

                    self.logger.info(f'Write GEMS L3 File : {output_filepath}')
                    self.write_product(output_filepath, result)
            else:
                raise FileNotFoundError(f"This File is Not GEMS L2 Product : {input_filepath}")

    def get_output_filename(self, input_filename, matcher, region_name):
        output_filename = input_filename

        output_filename = output_filename.replace(matcher.group("level"), "L3")
        output_filename = output_filename.replace(matcher.group("area"), f"GRID-{region_name.upper()}")
        output_filename = output_filename.replace(matcher.group("prosmode"), "hourly")
        output_filename = output_filename.replace(matcher.group("binning"), f"{int(self.grid_size * 100):03d}")

        return output_filename

    def write_product(self, output_filepath, result):
        image, spartial = self.grid_lon.shape

        ds = netCDF4.Dataset(output_filepath, 'w')

        ds.createDimension('image', image)
        ds.createDimension('spatial', spartial)

        geo_group = ds.createGroup('Geolocation Fields')

        var_lon = geo_group.createVariable('Longitude', self.grid_lon.dtype.str, dimensions=("image", "spatial"),
                                           zlib=True)
        var_lon[:] = self.grid_lon
        var_lat = geo_group.createVariable('Latitude', self.grid_lat.dtype.str, dimensions=("image", "spatial"),
                                           zlib=True)
        var_lat[:] = self.grid_lat

        data_group = ds.createGroup('Data Fields')
        for key, data in result.items():
            var_data = data_group.createVariable(key, data.dtype.str, dimensions=("image", "spatial"), zlib=True)
            var_data[:] = data

        ds.close()
        ds = None
