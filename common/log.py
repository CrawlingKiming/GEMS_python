import logging
import os
from logging.handlers import TimedRotatingFileHandler

backup_count = 30  # possible exceeded log file counts.


def init_logger(config):
    logger_config = config['logger']
    logger = logging.getLogger(logger_config['name'])
    logger.setLevel(logging.DEBUG)

    formatter = _get_formatter()

    sh = _create_handler(formatter, level=logging.DEBUG)
    fh = _create_handler(formatter, level=logging.DEBUG,
                         path=os.path.join(logger_config['path'], logger_config['name'].upper() + ".log"))

    logger.addHandler(sh)
    logger.addHandler(fh)

    return logger


def _create_handler(formatter, level, path=None):
    if path is None:
        handler = logging.StreamHandler()
    else:
        current_dir = os.curdir
        logger_filepath = os.path.join(current_dir, path)
        os.makedirs(os.path.dirname(logger_filepath), exist_ok=True)
        handler = TimedRotatingFileHandler(filename=logger_filepath, when='midnight', backupCount=backup_count)

    handler.setLevel(level)
    handler.setFormatter(formatter)

    return handler


def _get_formatter():
    return logging.Formatter('[%(asctime)s | %(levelname)s | %(filename)s:%(lineno)s | %(module)s] > %(message)s')


def get_stream_logger(name="default"):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = _get_formatter()

    sh = _create_handler(formatter, level=logging.DEBUG)

    logger.addHandler(sh)

    return logger

