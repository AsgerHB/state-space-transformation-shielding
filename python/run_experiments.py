import os
import numpy as np
import pandas as pd

from experiments.commands import shield_experiment, controller_experiment, \
    combined_experiment

np.set_printoptions(precision=3, suppress=True)

# this must be set appropriately in order to recompute the UPPAAL strategies
# os.environ['VERIFYTA_PATH'] = '/opt/uppaal-5.0.0-linux64/bin/verifyta'


# define data columns
shield_cols = [
    'input size',
    'mp size',
    'mp safe',
    'mp # unsafe runs',
    'mp violations',
    'viper size',
    'viper safe',
    'viper # unsafe runs',
    'viper violations'
]
control_cols = [
    'input size',
    'mp size',
    'mp # unsafe runs',
    'mp violations',
    'mp perf',
    'mp conf',
    'viper size',
    'viper # unsafe runs',
    'viper violations',
    'viper perf',
    'viper conf'
]
combined_cols = [
    'evil perf',
    'evil conf',
    'evil corrections',
    'evil violations',
    'norm perf',
    'norm conf',
    'norm corrections',
    'norm violations',
]


def op_callback(env):
    env.v = 10
    terminated = False
    return terminated


def bb_callback(env):
    env.p = 8 + np.random.uniform(0,2)
    env.v = 0
    terminated = False
    return terminated


def dc_callback(env):
    env.x1 = 0.35
    env.x2 = 15.0
    terminated = False
    return terminated


BASE_DIR = './experiments/models/'

# change to a non-existing strategy to compute a new oracle in UPPAAL
STRAT_NAME = 'unshielded_strategy.json'

# change to rebuild the shield with MaxPartitions
MP_SHIELD_NAME = 'mp_shield.json'

config = {
    'random_walk': [(True, True), None],
    'cruise': [(False, False), lambda x: False],
    'oil_pump': [(True, False), op_callback],
    'bouncing_ball': [(True, False), bb_callback],
    'dcdc_boost_converter': [(False, False), dc_callback]
}

all_shield, all_control, all_combined = [], [], []

k = 1
# models = ['random_walk', 'cruise', 'oil_pump', 'bouncing_ball', 'dcdc_boost_converter']
models = ['cruise']
for model in models:
    nstars = len(model) + 10
    print('*' * nstars)
    print(f"{' ' * 5}{model}{' ' * 5}")
    print('*' * nstars, '\n')

    # directory of the model
    model_dir = BASE_DIR + model + '/'

    # set configs
    unlucky, callback = config[model]
    env_kwargs = { 'unlucky': unlucky[0] }
    store_path = model_dir + 'shields/' + MP_SHIELD_NAME

    print('computing shield results')
    shield_res = []
    for _ in range(k):
        results, wrapped_shield, sdata = shield_experiment(
            model_dir, store_path,
            env_kwargs=env_kwargs, callback=callback
        )
        shield_res.append(results)

    # compile results
    shield_res = np.array(shield_res).T.mean(axis=1)
    shield_res[6] = int(shield_res[6] == 1)
    all_shield.append(shield_res.tolist())

    print('compute controller results')
    env_kwargs = { 'unlucky': unlucky[1] }
    cont_results, v_tree = controller_experiment(
        wrapped_shield, model_dir, sdata, STRAT_NAME,
        env_kwargs=env_kwargs, callback=callback
    )
    all_control.append(cont_results)

    print('compute combined results')
    comb_results = combined_experiment(v_tree, wrapped_shield, sdata)
    all_combined.append(comb_results)

# convert results to pandas
shield_df = pd.DataFrame(all_shield, index=models, columns=shield_cols)
control_df = pd.DataFrame(all_control, index=models, columns=control_cols)
combined_df = pd.DataFrame(all_combined, index=models, columns=combined_cols)

# store results
shield_df.to_csv('./results/shield_res.csv')
control_df.to_csv('./results/control_res.csv')
combined_df.to_csv('./results/combined_res.csv')

print('results stored in ./results/')
