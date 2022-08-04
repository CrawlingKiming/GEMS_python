import numpy as np
from shapely.geometry import Polygon

from regrid.abstract_regrid import AbstractRegrid


class Tessellation(AbstractRegrid):
    poly = []

    def _init_grid(self):
        diff_lon = np.nanmean(np.diff(self.lon, axis=1))
        diff_lat = np.nanmean(np.diff(self.lat, axis=0))

        self.lon = self.lon.ravel()
        self.lat = self.lat.ravel()

        nan_mask = np.isnan(self.lon) | np.isnan(self.lat)

        for key, data in self.datas.items():
            self.datas[key] = data.ravel()[~nan_mask]

        self.lon = self.lon[~nan_mask]
        self.lat = self.lat[~nan_mask]

        self.poly.clear()
        for i in range(len(self.lon)):
            poly = Polygon([
                (self.lon[i], self.lat[i]),
                (self.lon[i] + diff_lon, self.lat[i]),
                (self.lon[i] + diff_lon, self.lat[i] + diff_lat),
                (self.lon[i], self.lat[i] + diff_lat),
            ])

            self.poly.append(poly)

        self.poly = np.array(self.poly)
        # print("self.lon :", self.lon.shape)
        # print("self.lat :", self.lat.shape)

    def execute(self, new_lon: np.ndarray, new_lat: np.ndarray):
        new_nrows = new_lon.shape[0]
        new_ncols = new_lon.shape[1]
        # print("new_lon :", new_lon.shape)
        # print("new_lat :", new_lat.shape)

        new_dx = np.mean(np.diff(new_lon, axis=1))
        new_dy = np.mean(np.diff(new_lat, axis=0))
        # print("new_dx :", new_dx)
        # print("new_dy :", new_dy)

        result = dict()
        for key in self.datas.keys():
            result[key] = np.full((new_nrows, new_ncols), np.nan)

        # Step 3. 목표 격자 별로 Tessellation 적용
        for i in range(new_nrows):
            for j in range(new_ncols):
                lon0 = new_lon[i, j]
                lat0 = new_lat[i, j]

                # Step 3-1. 주변의 원본 데이터 파일 찾기
                idx = np.where(
                    (self.lon > (lon0 - new_dx)) &
                    (self.lon < (lon0 + (2 * new_dx))) &
                    (self.lat > (lat0 + new_dy)) &
                    (self.lat < (lat0 - (2 * new_dy)))
                )[0]
                # idx : (16,)
                # [ 9004  9005  9698  9699  9700  9701 10392 10393 10394 10395 11086 11087
                #  11088 11089 11782 11783]
                # if (i==20 & j==20):
                    # print("idx :", idx.shape)
                    # print(idx)

                if len(idx) == 0:
                    continue

                target_poly = Polygon([
                    (lon0, lat0),
                    (lon0 + new_dx, lat0),
                    (lon0 + new_dx, lat0 - new_dy),
                    (lon0, lat0 - new_dy),
                ])

                src_poly = self.poly[idx]

                # 0: k, 1: a, 2: s
                kas = []

                # Step 3-2. 목표 격자랑 원본데이터 Intersect & 면적 계산
                for k in range(len(idx)):
                    inter_poly = target_poly.intersection(src_poly[k])
                    if inter_poly.area == 0:
                        continue

                    kas.append([k, inter_poly.area, src_poly[k].area])

                # Step 3-3. 목표 격자 업데이트 (공식은 PPT 등 참조)
                if len(kas) == 0:
                    continue

                kas = np.array(kas) # (4,3)
                # 3-2에서 inter_poly.area==0이면 스킵해서 0, 1, ..., n-1까지 모두 있지는 않은 듯
                # [[2.00000000e+00 5.33196912e-04 8.92362208e-04]
                #  [3.00000000e+00 5.35633950e-04 8.92362208e-04]
                #  [6.00000000e+00 8.18980276e-04 8.92362208e-04]
                #  [7.00000000e+00 8.16543237e-04 8.92362208e-04]]

                w = kas[:, 1] / kas[:, 2] # (4,)
                # if (i==20 & j==20):
                    # print("kas :", kas.shape)
                    # print(kas)
                    # print("w :", w.shape)

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
            print("Regrid {} / {} done.".format(i, new_nrows), end='\r')
        return result

    # def arithmeticSequence(self, degree, min, max):
    #     start = float(min) + (degree / 2)
    #     end = float(max)  # -(degree/2)
    #     result = np.arange(start, end, degree)
    #     return result
