import yaml
from common.log import init_logger

def init(path):
    config = load_config(path)
    logger = init_logger(config)

    return config, logger

def load_config(path):
    with open(path, 'r', encoding="utf-8") as f:
        return yaml.load(f, Loader=yaml.FullLoader)
