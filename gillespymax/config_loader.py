import yaml

from os import PathLike


def load(config_file: str | PathLike):
    with open(config_file, "r") as stream:
        config = yaml.safe_load(stream)
    return config
