import os
from experiments.commands import run_from_config

os.environ['VERIFYTA_PATH'] = '/opt/uppaal-5.0.0-linux64/bin/verifyta'


if __name__ == '__main__':
    run_from_config('./experiments/config.json', './experiments/results/')
