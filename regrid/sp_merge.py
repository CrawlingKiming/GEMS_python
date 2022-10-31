import numpy as np

from regrid.abstract_regrid import AbstractRegrid


class SpatialMerge(AbstractRegrid):

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
        # self.qf = self.qf[~nan_mask]

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
        # NEED MORE ELABORATION ON NA -> How should we deal with it?

        # Data should be numpy array.
        grid_nrows = grid_lon.shape[0]
        grid_ncols = grid_lon.shape[1]

        result = dict()
        for key, data in self.datas.items():
            # for data, we need to ensure that nan is np.nan
            result[key] = np.full((grid_nrows, grid_ncols), np.nan)

            ### Calculate sigma_(1)
            sigma_1_dL = dict() # key dL/6 : matrix I X J
            sigma_1_N = dict() # key dL/6 : matrix I X J

            for dL in self.delta_L:
                thr = int(dL/6)
                # default zeros
                temp_sig = np.zeros(shape=(grid_nrows, grid_ncols))
                num_sig = np.zeros(shape=(grid_nrows, grid_ncols))

                for i in range(grid_nrows):
                    for j in range(grid_ncols):
                        # Center with i,j, generates donut grid. Shift Center.
                        adv_idx = self.adv_idx[thr]
                        idx = adv_idx + [i, j] # shifts center
                        # kills inappropriate indexes + considers marginal indexes
                        m1 = (idx >= 0.)
                        m2 = (np.sum(m1, axis=1) == 2)

                        src_data = data[idx[m2, 0], idx[m2, 1]]
                        if np.isnan(src_data).all():
                            continue
                        nan_mask = np.isnan(src_data)
                        num_sig[i, j] = np.sum(~nan_mask) # num. of data used, erase NA
                        temp_sig[i, j] = np.sum(np.square(src_data[~nan_mask] - data[i, j]))

                sigma_1_dL[thr] = temp_sig
                sigma_1_N[thr] = num_sig
            ###############################################################################################
            # Calculate sigma_1
            # default zero
            sigma_1 = np.zeros(shape=(grid_nrows, grid_ncols))
            sigma_N = np.zeros(shape=(grid_nrows, grid_ncols))
            for dL in self.delta_L:
                thr = int(dL/6)
                sigma_1 += sigma_1_dL[thr]
                sigma_N += sigma_1_N[thr]
            sigma_1_squ = sigma_1 / sigma_N

            # Calculate AOD_est
            # 이거 계산할 떄 sigma_zero 필요 없는 것 맞지?
            weight_est = np.zeros(shape=(grid_nrows, grid_ncols))
            sig_est = np.zeros(shape=(grid_nrows, grid_ncols))
            for i in range(grid_nrows):
                for j in range(grid_ncols):
                    idx = self.N_idx + [i, j]
                    m1 = (idx >= 0.)
                    m2 = (np.sum(m1, axis=1) == 2)
                    sig_est[i, j] = np.reciprocal(np.sum(np.reciprocal(sigma_1_squ[idx[m2, 0], idx[m2, 1]], -1)), -1)
                    weight_est[i, j] = sig_est[i, j] / sigma_1_squ[i, j]

            # Second-Order Regression Models ##############################################################################
            # sigma_1 oc/land dist categorize
            # linear models
            # get intercepts
            # return sigma_zero (the intercept) matrix, I X J NOTE THIS MATRIX SHOULD BE SQUARE, NOT SQRT
            sigma_zero = np.zeros(shape=(grid_nrows, grid_ncols))
            ###############################################################################################################

            AOD_est = np.zeros(shape=(grid_nrows, grid_ncols))
            for i in range(grid_nrows):
                for j in range(grid_ncols):
                    idx = self.N_idx + [i, j]
                    m1 = (idx >= 0.)
                    m2 = (np.sum(m1, axis=1) == 2)
                    AOD_est[i, j] = weight_est[i, j] * data[idx[m2, 0], idx[m2, 1]]

            sigma_pure_squ = sigma_zero + sig_est
            AOD_pure = np.full((grid_nrows, grid_ncols), np.nan)
            # Treated as nan (including as na Value)
            condition_mask = (data < (AOD_est + 2.58 * np.sqrt(sigma_pure_squ)))
            AOD_pure[condition_mask] = data[condition_mask]

            # Calculate AOD_merged
            sigma_merged = np.zeros(shape=(grid_nrows, grid_ncols))
            weight_final = np.zeros(shape=(grid_nrows, grid_ncols))
            for i in range(grid_nrows):
                for j in range(grid_ncols):
                    idx = self.N_idx + [i, j]
                    m1 = (idx >= 0.)
                    m2 = (np.sum(m1, axis=1) == 2)

                    sigma_merged[i, j] = np.reciprocal(np.sum(np.reciprocal(sigma_pure_squ[idx[m2, 0], idx[m2, 1]], -1)), -1)
                    weight_final[i, j]= sigma_merged[i, j] / sigma_pure_squ[i, j]

            for i in range(grid_nrows):
                for j in range(grid_ncols):
                    idx = self.N_idx + [i, j]
                    m1 = (idx >= 0.)
                    m2 = (np.sum(m1, axis=1) == 2)
                    # Erase NA
                    rs = AOD_pure[idx[m2, 0], idx[m2, 1]]
                    cc = np.isnan(rs)
                    result[key][i, j] = weight_final[idx[m2, 0], idx[m2, 1]][~cc] * rs[~cc]

        return result


