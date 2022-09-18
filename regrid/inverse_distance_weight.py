import numpy as np

from regrid.abstract_regrid import AbstractRegrid


class InverseDistanceWeight(AbstractRegrid):

    def __init__(self, lon: np.ndarray, lat: np.ndarray, datas: dict, qf: np.ndarray = None, radius:float=0.1, p: float = 2):
        super().__init__(lon, lat, datas, qf)
        self.radius = radius
        self.p = p

    def _init_grid(self):
        self.lon = self.lon.ravel()
        self.lat = self.lat.ravel()
        self.qf = self.qf.ravel()

        nan_mask = np.isnan(self.lon) | np.isnan(self.lat)

        for key, data in self.datas.items():
            self.datas[key] = data.ravel()[~nan_mask]
            # NC datas

        # self.lon / self.lat : L2 grid
        self.lon = self.lon[~nan_mask]  # lon : (1421312,)
        self.lat = self.lat[~nan_mask]  # lat : (1421312,)
        self.qf = self.qf[~nan_mask]

    def execute(self, grid_lon: np.ndarray, grid_lat: np.ndarray):
        grid_nrows = grid_lon.shape[0]
        grid_ncols = grid_lon.shape[1]

        # 얘네 뭐하는 애들인지 모르겠음.
        # lon0 / lat0에서 얼만큼 +-로 subset할지를 이걸로 정한건가?
        # grid_dx = np.mean(np.diff(grid_lon, axis=1))
        # grid_dy = np.mean(np.diff(grid_lat, axis=0))

        result = dict()
        for key in self.datas.keys():
            result[key] = np.full((grid_nrows, grid_ncols), np.nan)

        # Step 3. IDW
        for i in range(grid_nrows):
            for j in range(grid_ncols):
                lon0 = grid_lon[i, j]
                lat0 = grid_lat[i, j]

                # idx : Desired L2 points matrix
                # weight : Desired unnormalized weights of L2
                # @TODO : clarify role of degree / dx / dy

                idx = np.where(
                    (self.lon > (lon0 - self.radius)) &
                    (self.lon < (lon0 + self.radius)) &
                    (self.lat < (lat0 + self.radius)) &
                    (self.lat > (lat0 - self.radius)))[0]

                if len(idx) == 0:  # intersection에서 걸리는게 아무것도 없으면 skip
                    continue

                x0 = np.array([lon0, lat0])
                lon_sub = self.lon[idx]
                lat_sub = self.lat[idx]
                qf_sub = self.qf[idx]

                d = np.linalg.norm((x0[0] - lon_sub, x0[1] - lat_sub), axis=0)

                d_inv = 1 / ((d ** self.p) * qf_sub)

                for key, data in self.datas.items():
                    src_data = data[idx]
                    if np.isnan(src_data).all():
                        continue

                    nan_mask = np.isnan(src_data)

                    weight_idx = d_inv[~nan_mask]
                    z = src_data[~nan_mask]

                    result[key][i, j] = np.sum(weight_idx * z) / np.sum(weight_idx)

        return result


