# -*- encoding: utf8 -*-
import cProfile
import datetime
import sys
import traceback

import line_profiler

from common.common import init
from common.profile_utils import profile_report, line_profile_report
from process.L2convertProcess import L2convertProcess
from regrid.abstract_regrid import AbstractRegrid
from regrid.tessellation import Tessellation


def parse_args():
    """
    입력 파라미터를 파싱하는 메소드
    :return: 파싱한 입력 파라미터 구조체
    """
    if len(sys.argv) > 1:
        yr = int(sys.argv[1])
        mm = int(sys.argv[2])
        dd = int(sys.argv[3])
        hh = int(sys.argv[4])
        dt = datetime.datetime(year=yr, month=mm, day=dd, hour=hh)
    else:
        dt = datetime.datetime.now(datetime.timezone.utc)
        dt = datetime.datetime(year=dt.year, month=dt.month, day=dt.day, hour=dt.hour)

    return dt


def main():
    config, logger = init("config.yml")
    logger.info("Start GEMS L3 Tessellation Process")

    dt = parse_args()
    # 실행 예 --> python Application.py 2022 01 06 00
    # dt = datetime.datetime(2022,1,6,0,0,0) # for Test
    logger.debug(f"Input Args. : {dt}")

    startTime = datetime.datetime.now()

    try:
        worker = L2convertProcess(config)
        worker.run(dt)
    except:
        logger.error(traceback.format_exc())

    logger.debug("running time : {} s".format((datetime.datetime.now() - startTime).total_seconds()))
    logger.info("Done GEMS L3 Tessellation Process")


def profile():
    pr = cProfile.Profile(builtins=False)
    pr.enable()
    main()
    pr.disable()

    profile_report(pr.getstats())


def line_profile():
    pr = line_profiler.LineProfiler()
    pr.enable()
    pr.add_function(L2convertProcess.run)
    pr.add_function(L2convertProcess.write_product)
    pr.add_function(Tessellation.execute)
    pr.add_function(Tessellation._init_grid)
    pr.add_function(AbstractRegrid.__init__)
    # pr.add_function(BrowserProcess._l2_process)
    # pr.add_function(BrowserProcess._apply_renderer)
    # pr.add_function(BrowserProcess._add_accessory)
    # pr.add_function(add_colorbar_v2)
    # pr.add_function(compute_moving_average_exp)
    # pr.add_function(get_geos_bbox)
    # pr.add_function(MosaicProcess.Run)
    main()
    pr.disable()

    line_profile_report(pr.get_stats())


if __name__ == '__main__':
    # line_profile()
    # profile()
    main()
