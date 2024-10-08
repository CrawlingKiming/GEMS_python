import argparse
import datetime
import traceback

from common.common import init
from process.average_process import AverageProcess
from process.SP_average import SptAverageProcess


def parse_args():
    """
    입력 파라미터를 파싱하는 메소드
    :return: 파싱한 입력 파라미터 구조체
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--date", type=lambda s: datetime.datetime.strptime(s, '%Y%m%d'), required=True)
    parser.add_argument("-m", "--mode", choices=["yearly", "monthly", "daily"], default="daily", required=True)
    parser.add_argument("-t", "--target", choices=["NO2", "AERAOD"], required=True)

    options = parser.parse_args()

    return options


def main():
    config, logger = init("config.yml")
    logger.info("Start GEMS L3 Average Process")

    options = parse_args()

    startTime = datetime.datetime.now()

    try:
        if options.target == "AERAOD":
            print("AERAOD spt running")
            worker = SptAverageProcess(config)
        else:
            worker = AverageProcess(config)
        worker.run(options.date, options.mode, options.target)
    except:
       logger.error(traceback.format_exc())

    logger.debug("running time : {} s".format((datetime.datetime.now() - startTime).total_seconds()))
    logger.info("Done GEMS L3 Average Process")

if __name__ == '__main__':
    main()
