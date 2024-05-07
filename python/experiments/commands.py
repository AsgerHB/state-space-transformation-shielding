import json
import pathlib

from experiments.utils import Shield, UPPAALInstance


def run_from_config(config_path, results_dir):
    with open(config_path, 'r') as f:
        conf = json.load(f)

    for model in conf:
        if model != 'cartpole':
            continue

        model_file = conf[model]['model_file']
        state_vars = conf[model]['state_vars']
        shields = conf[model]['shields']
        training_costs = conf[model]['training_costs']
        training_bound = conf[model]['training_bound']
        training_goal = conf[model]['training_goal']
        queries = conf[model]['queries']

        model_slug = model.replace(' ', '_')
        model_results = []

        instance = UPPAALInstance(model_file)

        for cost in training_costs:
            for svars in state_vars:
                features = '{} -> {' + ','.join(svars) + '}'
                instance.add_training_query(
                    training_bound, features, training_goal, cost=cost
                )
                for query in queries:
                    instance.add_query(query)

        for shield_name in shields:
            model_results.append(shield_name)

            shield_path = shields[shield_name]['path']

            # run with no shield
            if shield_path is None:
                signature =  'int shield() { return 0; }'
                instance.update_tag('ENABLE_SHIELD', 'false')
                instance.update_tag('SHIELD_PATH', '')
                instance.update_tag('SHIELD_SIGNATURE', signature)
                instance.update_tag('SHIELD_CALL', 'shield()')

            # run with shield
            else:
                shield_call = shields[shield_name]['shield_call']

                shield = Shield(shield_path, make_tree=True)
                shield_signature = shield.get_signature()
                so_path = model_slug + f'_{shield_name}.so'
                shield.compile_so(so_path)

                shield_path = pathlib.Path(so_path).resolve()
                text = f'import "{shield_path}"'
                instance.update_tag('ENABLE_SHIELD', 'true')
                instance.update_tag('SHIELD_PATH', text)
                instance.update_tag('SHIELD_SIGNATURE', shield_signature + ';')
                instance.update_tag('SHIELD_CALL', shield_call)

            # run model and add results to results list
            output = instance.run('-Wsqy')
            parsed_output = instance.parse_ouput(output)
            model_results.append(parsed_output)

        results_path = pathlib.Path(results_dir + f'/{model_slug}_results.txt')
        with open(results_path, 'w') as f:
            for row in model_results:
                if isinstance(row, str):
                    f.write(row + '\n')
                else:
                    f.write('\t' + '\n\t'.join(map(str, row)) + '\n')
        print(f"results written to '{results_path}'")
