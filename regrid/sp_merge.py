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

        # self.lon / self.lat : Now we have L3 grid
        self.lon = self.lon[~nan_mask]  # lon : (1421312,)
        self.lat = self.lat[~nan_mask]  # lat : (1421312,)
        self.qf = self.qf[~nan_mask]

        # for AOD
        self.delta_L = [6, 12, 18, 24]
        self.adv_idx = self.s_idx()

    def s_idx(self):
        """
        Generates (donout or marginal) grid w.r.t. center (0,0)
        """
        clist = [set((0, 0))]
        for nu in range(6):
            xv, yv = np.meshgrid(np.linspace(-nu, nu, num=nu * 2 + 1), np.linspace(-nu, nu, num=nu * 2 + 1),
                                 indexing='ij')
            X = xv.reshape((np.prod(xv.shape),))
            Y = yv.reshape((np.prod(yv.shape),))
            coords = set(zip(X, Y))
            clist.append(coords)
        res = [x - y for x, y in zip(clist[1:], clist)]
        return res[1:4]

    def execute(self, grid_lon: np.ndarray, grid_lat: np.ndarray):

        ### Ensuring NA to be zero.. ?
        # 결측치 처리 어캐 하는거? not sure..
        # fectch data
        grid_nrows = grid_lon.shape[0]
        grid_ncols = grid_lon.shape[1]

        ### Calculate sigma_(1)
        sigma_1_dL = dict() # key dL/6 : matrix I X J
        sigma_1_N = dict() # key dL/6 : matrix I X J
        for dL in self.delta_L:
            thr = int(dL/6)
            ## Is this shape correct? I am not sure
            temp_sig = np.zeros(shape=(grid_nrows, grid_ncols))#sigma_1[thr]
            num_sig = np.zeros(shape=(grid_nrows, grid_ncols))

            for i in range(4, grid_nrows-4): # what happens to the corrdinate outside the 4, grid_nrows-4
                for j in range(4, grid_ncols-4):
                        adv_idx = self.adv_idx[thr]
                        # center with i,j, generates donut grid.
                        idx = np.asarray(list(adv_idx), dtype=np.int32) + [i,j]
                        sigma = data[idx[:, 0], idx[:, 1]]
                        # 여기서 예외 처리를 하고, N 갯수를 고려해야하나?? 그러면 sum 불가..
                        temp_sig[i,j] = np.sum(np.square(sigma - data[i,j]))
                        num_sig = N(i,j) # 여기서 각 좌표별 결칙치 제거한 N의 갯수를 저장 ( 결측치 제거가 필요 없다면, N 일정 )

            sigma_1_dL[thr] = temp_sig
            sigma_1_dL[thr] = num_sig

        #calculate sigma_1
        sigma_1 = np.zeros(shape=(grid_nrows, grid_ncols))
        sigma_N= np.zeros(shape=(grid_nrows, grid_ncols))
        for dL in self.delta_L:
            thr = int(dL/6)
            sigma_1 += sigma_1_dL[thr] * dL
            sigma_N += num_sig[thr]
        sigma_1 = sigma_1 / sigma_N

        # Regression Models ###########################################################################################
        # sigma_1 oc/land dist categorize
        # linear models
        ###############################################################################################################


        # 얘네 뭐하는 애들인지 모르겠음.
        # lon0 / lat0에서 얼만큼 +-로 subset할지를 이걸로 정한건가?
        # grid_dx = np.mean(np.diff(grid_lon, axis=1))
        # grid_dy = np.mean(np.diff(grid_lat, axis=0))


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


