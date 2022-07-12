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

    def execute(self, new_lon: np.ndarray, new_lat: np.ndarray):
        new_nrows = new_lon.shape[0]
        new_ncols = new_lon.shape[1]

        new_dx = np.mean(np.diff(new_lon, axis=1))
        new_dy = np.mean(np.diff(new_lat, axis=0))

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

                kas = np.array(kas)

                w = kas[:, 1] / kas[:, 2]

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
