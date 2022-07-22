import numpy as np

from .abstract_regrid import AbstractRegrid


class InverseDistanceWeight(AbstractRegrid):

    def _init_grid(self):
        self.lon = self.lon.ravel()
        self.lat = self.lat.ravel()

        nan_mask = np.isnan(self.lon) | np.isnan(self.lat)

        for key, data in self.datas.items():
            self.datas[key] = data.ravel()[~nan_mask]
            # NC datas

        # self.lon / self.lat : L2 grid
        self.lon = self.lon[~nan_mask]  # lon : (1421312,)
        self.lat = self.lat[~nan_mask]  # lat : (1421312,)

    def execute(self, new_lon: np.ndarray, new_lat: np.ndarray):
        print("IDW START")
        # new_lon / new_lat : L3 grid
        # new_lon : (1040, 1320)
        # new_lat : (1040, 1320)
        new_nrows = new_lon.shape[0]
        new_ncols = new_lon.shape[1]

        # 얘네 뭐하는 애들인지 모르겠음.
        # lon0 / lat0에서 얼만큼 +-로 subset할지를 이걸로 정한건가?
        new_dx = np.mean(np.diff(new_lon, axis=1))
        new_dy = np.mean(np.diff(new_lat, axis=0))

        result = dict()
        for key in self.datas.keys():
            result[key] = np.full((new_nrows, new_ncols), np.nan)

        # Step 3. IDW
        for i in range(new_nrows):
            for j in range(new_ncols):
                lon0 = new_lon[i, j]
                lat0 = new_lat[i, j]

                # idx : Desired L2 points matrix
                # weight : Desired unnormalized weights of L2
                # @TODO : clarify role of degree / dx / dy
                degree = 0.05
                p = 2
                dx = 0.03#0.06
                dy = 0.015#0.03 #/ 2
                # intersection==0인 경우 continue로 넘겨야 해서 따로 함수 안 만들었음
                # idx, weight = self.l3_calculate(lon0, lat0, degree=0.05,
                #                                 p=2, dx=new_dx, dy=new_dy)
                # 주변의 원본 데이터 (L2 : lon / lat) 값 찾기
                idx = np.where(
                    (self.lon > (lon0 - dx)) &
                    (self.lon < (lon0 + degree + dx)) &
                    (self.lat > (lat0 - dy)) &
                    (self.lat < (lat0 + degree + dy))
                    )[0]


                if len(idx) == 0:  # intersection에서 걸리는게 아무것도 없으면 skip
                    continue

                # 0 : k / 1 : d_inv
                weight = []

                x0 = np.array([lon0, lat0])
                lon_sub = self.lon[idx]
                lat_sub = self.lat[idx]

                # 사실 아래 loop도 exception 없어서 matrix 연산으로 바꿔도 될듯
                for k in range(len(idx)):
                    xk = np.array([lon_sub[k], lat_sub[k]])
                    d = np.linalg.norm(x0 - xk)
                    weight.append([k, d])

                weight = np.array(weight)
                d_inv = 1 / weight[:, 1] ** p

                for key, data in self.datas.items():
                    src_data = data[idx]
                    if np.isnan(src_data).all():
                        continue
                    # 아래 weight[:,0]은 따로 생각한 weight dimension이 있어서 이렇게 한건지
                    # 아니면 tessellation에서 kas의 형식을 그대로 따온건지..!?
                    z = src_data[weight[:, 0].astype(int)]
                    if np.isnan(z).all():
                        continue
                    nan_mask = np.isnan(z)

                    weight_idx = d_inv[~nan_mask]
                    z_idx = z[~nan_mask]

                    result[key][i, j] = np.sum(weight_idx * z_idx) / np.sum(weight_idx)

            print("Regrid {} / {} done.".format(i, new_nrows), end='\r')
        return result


