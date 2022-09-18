import numpy as np


def make_grid(x0, x1, dx, y0, y1, dy=None, ref="c"):
    if dy is None:
        dy = -dx

    x = np.arange(x0, x1, dx)
    y = np.arange(y0, y1, dy)

    if ref.startswith("c"):
        x += dx / 2.0
        y += dy / 2.0

    return np.meshgrid(x, y)


def half_scene_avg(data: np.ndarray):
    _, ncol = data.shape

    if ncol % 2 == 0:
        return (data[:, 0::2] + data[:, 1::2]) / 2
    else:
        return (data[:, 0:-1:2] + data[:, 1:-1:2]) / 2
