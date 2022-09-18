import numpy as np
from shapely.geometry import Polygon

from regrid.abstract_regrid import AbstractRegrid


class Tessellation(AbstractRegrid):
    poly = []

    def _init_grid(self):
        self.dx = np.nanmean(np.diff(self.lon, axis=1))
        self.dy = np.nanmean(np.diff(self.lat, axis=0))

        self.lon = self.lon.ravel()
        self.lat = self.lat.ravel()
        self.qf = self.qf.ravel()

        nan_mask = np.isnan(self.lon) | np.isnan(self.lat)

        for key, data in self.datas.items():
            self.datas[key] = data.ravel()[~nan_mask]

        self.lon = self.lon[~nan_mask]
        self.lat = self.lat[~nan_mask]
        self.qf = self.qf[~nan_mask]

        self.poly.clear()
        for i in range(len(self.lon)):
            poly = Polygon([
                (self.lon[i] - self.dx / 2, self.lat[i] - self.dy / 2),
                (self.lon[i] + self.dx / 2, self.lat[i] - self.dy / 2),
                (self.lon[i] + self.dx / 2, self.lat[i] + self.dy / 2),
                (self.lon[i] - self.dx / 2, self.lat[i] + self.dy / 2),
            ])

            self.poly.append(poly)

        self.poly = np.array(self.poly)

    def execute(self, grid_lon: np.ndarray, grid_lat: np.ndarray):
        grid_nrows = grid_lon.shape[0]
        grid_ncols = grid_lon.shape[1]

        grid_dx = np.mean(np.diff(grid_lon, axis=1))
        grid_dy = np.mean(np.diff(grid_lat, axis=0))

        result = dict()
        for key in self.datas.keys():
            result[key] = np.full((grid_nrows, grid_ncols), np.nan)

        # Step 3. 목표 격자 별로 Tessellation 적용
        for i in range(grid_nrows):
            for j in range(grid_ncols):
                lon0 = grid_lon[i, j]
                lat0 = grid_lat[i, j]

                # Step 3-1. 주변의 원본 데이터 파일 찾기
                idx = np.where(
                    (self.lon > (lon0 - (2 * grid_dx))) &
                    (self.lon < (lon0 + (2 * grid_dx))) &
                    (self.lat < (lat0 - (2 * grid_dy))) &
                    (self.lat > (lat0 + (2 * grid_dy)))
                )[0]

                if len(idx) == 0:
                    continue

                target_poly = Polygon([
                    (lon0 - grid_dx / 2, lat0 - grid_dy / 2),
                    (lon0 + grid_dx / 2, lat0 - grid_dy / 2),
                    (lon0 + grid_dx / 2, lat0 + grid_dy / 2),
                    (lon0 - grid_dx / 2, lat0 + grid_dy / 2),
                ])

                src_poly = self.poly[idx]
                qf_sub = self.qf[idx]

                # 0: k, 1: a, 2: s
                kas = []

                # Step 3-2. 목표 격자랑 원본데이터 Intersect & 면적 계산
                for k in range(len(idx)):
                    inter_poly = target_poly.intersection(src_poly[k])
                    if inter_poly.area == 0:
                        continue

                    kas.append([k, inter_poly.area, src_poly[k].area, qf_sub[k]])

                # Step 3-3. 목표 격자 업데이트 (공식은 PPT 등 참조)
                if len(kas) == 0:
                    continue

                kas = np.array(kas) # (4,3)
                # 3-2에서 inter_poly.area==0이면 스킵해서 0, 1, ..., n-1까지 모두 있지는 않은 듯
                # [[2.00000000e+00 5.33196912e-04 8.92362208e-04]
                #  [3.00000000e+00 5.35633950e-04 8.92362208e-04]
                #  [6.00000000e+00 8.18980276e-04 8.92362208e-04]
                #  [7.00000000e+00 8.16543237e-04 8.92362208e-04]]

                w = kas[:, 1] / (kas[:, 2] * kas[:, 3]) # (4,)


                for key, data in self.datas.items():
                    src_data = data[idx]
                    if np.isnan(src_data).all():
                        continue

                    c = src_data[kas[:, 0].astype(int)]
                    if np.isnan(c).all():
                        continue
                    nan_mask = np.isnan(c)

                    wc = w[~nan_mask]
                    cc = c[~nan_mask]

                    result[key][i, j] = np.sum(wc * cc) / np.sum(wc)

        return result
