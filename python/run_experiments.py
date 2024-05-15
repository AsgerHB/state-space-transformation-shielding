import os
from experiments.commands import run_from_config

os.environ['VERIFYTA_PATH'] = '/opt/uppaal-5.0.0-linux64/bin/verifyta'

config_files = [
    './experiments/configs/bouncing_ball.json',
    './experiments/configs/cartpole.json',
    './experiments/configs/spiral.json'
]

if __name__ == '__main__':
    for conf in config_files:
        run_from_config(conf, './experiments/results/')
