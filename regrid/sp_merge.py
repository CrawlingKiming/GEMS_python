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

        # For AOD
        self.delta_L = [6, 12, 18, 24]
        # For indexes.
        # adv_idx stores marginal indexes, N_idx stores total grid index inside the boundary
        self.adv_idx, self.N_idx = self.super_idx()

    def super_idx(self):
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
        res = [np.asarray(list(x-y), dtype=np.int32) for x, y in zip(clist[1:], clist)]
        return res[1:5], np.asarray(list(clist[-2]), dtype=np.int32)

    def execute(self, grid_lon: np.ndarray, grid_lat: np.ndarray):
        # Data should be numpy array.

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
            ## Is this shape (I,J) correct? I am not sure
            temp_sig = np.zeros(shape=(grid_nrows, grid_ncols))#sigma_1[thr]
            num_sig = np.zeros(shape=(grid_nrows, grid_ncols))

            for i in range(4, grid_nrows-4): # what happens to the corrdinate outside the 4, grid_nrows-4
                for j in range(4, grid_ncols-4):
                        adv_idx = self.adv_idx[thr]
                        # Center with i,j, generates donut grid. Shift Center.
                        idx = adv_idx + [i,j]
                        sigma = data[idx[:, 0], idx[:, 1]]
                        # 여기서 예외 처리를 하고, N 갯수를 고려해야하나?? 그러면 sum 불가..
                        temp_sig[i,j] = np.sum(np.square(sigma - data[i,j]))
                        num_sig = N(i,j) # 여기서 각 좌표별 결칙치 제거한 N의 갯수를 저장 ( 결측치 제거가 필요 없다면, N 일정 )
            sigma_1_dL[thr] = temp_sig
            sigma_1_dL[thr] = num_sig

        # Calculate sigma_1
        sigma_1 = np.zeros(shape=(grid_nrows, grid_ncols))
        sigma_N= np.zeros(shape=(grid_nrows, grid_ncols))
        for dL in self.delta_L:
            thr = int(dL/6)
            sigma_1 += sigma_1_dL[thr]
            sigma_N += num_sig[thr]
        sigma_1_squ = sigma_1 / sigma_N

        # Calculate AOD_est
        # 이거 계산할 떄 sigma_zero 필요 없는 것 맞지?
        weight_est = np.zeros(shape=(grid_nrows, grid_ncols))
        sig_est = np.zeros(shape=(grid_nrows, grid_ncols))
        for i in range(4, grid_nrows - 4):  # what happens to the corrdinate outside the 4, grid_nrows-4
            for j in range(4, grid_ncols - 4):
                idx = self.N_idx + [i, j]
                sig_est[i, j] = np.reciprocal(np.sum(np.reciprocal(sigma_1_squ[idx[:, 0], idx[:, 1]], -1)), -1)
                weight_est[i, j] = sig_est[i, j] / sigma_1_squ[i, j]

        AOD_est = np.zeros(shape=(grid_nrows, grid_ncols))
        for i in range(4, grid_nrows - 4):  # what happens to the corrdinate outside the 4, grid_nrows-4
            for j in range(4, grid_ncols - 4):
                idx = self.N_idx + [i, j]
                AOD_est[i,j] = data[idx[:,0], idx[:,1]]

        # Regression Models ###########################################################################################
        # sigma_1 oc/land dist categorize
        # linear models
        # get intercepts
        # return sigma_zero (the intercept) matrix, I X J NOTE THIS MATRIX SHOULD BE SQUARE, NOT SQRT
        sigma_zero = np.zeros(shape=(grid_nrows, grid_ncols))
        ###############################################################################################################

        # Calculate sigma_pure & AOD_pure
        # 여기서 원래 IDW에서 사용한 결측치 기준을 적용하면 안 될 것 같은데?
        sigma_pure = np.squrt(sigma_zero + sig_est)
        # In AOD_pure, M.A. could be treated as zero. (GUESS_NEED DOUBLE CHECK)
        AOD_pure = data
        condition_mask = (data > (AOD_est + 2.58 * sigma_pure))
        AOD_pure[condition_mask] = 0.0

        # Calculate AOD_merged
        sigma_merged = np.zeros(shape=(grid_nrows, grid_ncols))
        for i in range(4, grid_nrows - 4):  # what happens to the corrdinate outside the 4, grid_nrows-4
            for j in range(4, grid_ncols - 4):
                idx = self.N_idx + [i, j]
                sigma_merged[i, j] = np.reciprocal(np.sum(np.reciprocal(sigma_pure[idx[:, 0], idx[:, 1]], -1)), -1)
                weight_est[i, j] = sigma_merged[i, j] / sigma_pure[i, j]

        AOD_merged = np.zeros(shape=(grid_nrows, grid_ncols))
        for i in range(4, grid_nrows - 4):  # what happens to the corrdinate outside the 4, grid_nrows-4
            for j in range(4, grid_ncols - 4):
                idx = self.N_idx + [i, j]
                AOD_merged[i,j] = AOD_pure[idx[:,0], idx[:,1]]

        return AOD_merged


