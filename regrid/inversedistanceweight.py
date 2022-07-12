from regrid.abstract_regrid import AbstractRegrid


class InverseDistanceWeight(AbstractRegrid):
    def _init_grid(self):
        raise NotImplementedError("This Method is Not Available!")

    def execute(self, new_lon, new_lat):
        raise NotImplementedError("This Method is Not Available!")
