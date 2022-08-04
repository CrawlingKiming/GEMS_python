import logging
import tempfile
from abc import ABCMeta, abstractmethod


class Process(metaclass=ABCMeta):
    config = None
    logger = None
    temp_dir = None

    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger(config['logger']['name'])
        self.temp_dir = tempfile.TemporaryDirectory()

    def __del__(self):
        self.temp_dir.cleanup()

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    @abstractmethod
    def _init_process(self, *args, **kwargs):
        pass

    @abstractmethod
    def run(self, *args, **kwargs):
        pass
