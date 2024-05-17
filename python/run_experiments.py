import os
from experiments.commands import run_from_config

os.environ['VERIFYTA_PATH'] = '/opt/uppaal-5.0.0-linux64/bin/verifyta'

config_files = [
    './experiments/configs/bouncing_ball.json',
    './experiments/configs/cartpole.json',
    './experiments/configs/spiral.json'
]

if __name__ == '__main__':
    for i in range(1, 11):
        print(f"================{i}================")
        results_dir = f'./experiments/results/{i}/'
        if not os.path.isdir(results_dir):
            os.mkdir(results_dir)
        for conf in config_files:
            run_from_config(conf, results_dir)
