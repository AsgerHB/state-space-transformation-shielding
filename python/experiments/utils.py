import os
import re
import json
import shutil
import pathlib
import tempfile
import subprocess
import numpy as np

from tempfile import NamedTemporaryFile

from trees.models import DecisionTree
from trees.advanced import minimize_tree


class Shield:
    def __init__(self, path, make_tree=False, **kwargs):
        """
        Initialize shield from zip found in `path`. Set `make_tree=True` for
        automatically generating a decision tree. `kwargs` are passed on to
        the `make_tree()` function if `make_tree=True`.
        """
        grid, meta = self.unpack_zip(path)

        self.grid = grid
        self.variables = meta['variables']
        self.env_kwargs = meta.get('env_kwargs', {})
        self.bounds = np.array(meta['bounds'])
        self.granularity = meta['granularity']
        self.n_actions = meta['n_actions']
        self.id_to_actionset = meta['id_to_actionset']
        self.env_id = meta['env_id']
        self.bvars = meta['bvars']
        self.amap = self.build_action_map()

        if make_tree:
            tree = self.make_tree(**kwargs)

    def build_action_map(self):
        """build action map (what shield action corresponds to which allowed
        policy actions?)"""
        amap = [0] * len(self.id_to_actionset)
        for i, acts in self.id_to_actionset.items():
            amap[int(i)] = tuple(np.argwhere(acts).T[0])

        return amap

    def make_tree(self, minimize=True, verbose=False, **kwargs):
        tree = DecisionTree.from_grid(
            self.grid,
            self.variables,
            np.arange(len(self.id_to_actionset)),
            self.granularity,
            self.bounds.T,
            bvars=self.bvars
        )
        if minimize:
            tree, _ = minimize_tree(tree, verbose=verbose, **kwargs)

        self.tree = tree

    def get_signature(self, function_name='shield'):
        params = ', '.join([f'double {v}' for v in self.variables])
        return f'int {function_name}({params})'

    def compile_so(self, so_path, c_path=None):
        with tempfile.TemporaryDirectory() as tmpdirname:
            if c_path is None:
                c_path, o_path = tmpdirname + '/s.c', tmpdirname + '/s.o'
            else:
                o_path = tmpdirname + '/s.o'

            self.tree.export_to_c_code(signature=self.get_signature(), out_fp=c_path)

            # make object
            args = ('gcc', '-c', '-fPIC', c_path, '-o', o_path)
            subprocess.run(args)

            # make shared object
            args = ('gcc', '-shared', '-o', so_path, o_path)
            subprocess.run(args)

    def unpack_zip(self, path):
        with tempfile.TemporaryDirectory() as tmpdirname:
            shutil.unpack_archive(path, tmpdirname)
            with open(tmpdirname + '/meta.json', 'r') as f:
                meta = json.load(f)

            grid = np.load(tmpdirname + '/grid.npy')

        return grid, meta


class UPPAALInstance:
    def __init__(self, model_file, cost=None):
        self.model_file = model_file
        self.queries = []
        self.cost = cost
        self.reload_model()

    def reload_model(self):
        with open(self.model_file, 'r', newline='') as f:
            self.content = f.read()

    def extract_tags(self):
        pattern = rf'// \[ (\w+)(?:\r\n|\r|\n)([\S\s]+?)(?=// ])'
        return dict(re.findall(pattern, self.content))

    def update_tag(self, tag, value):
        pattern = rf'(// \[ {tag}(?:\r\n|\r|\n))([\S\s]+?)(?=// ])'
        repl = rf'\1{value}\n'
        result, subs = re.subn(pattern, repl, self.content)

        if subs == 0:
            print(f'tag `{tag}` not found, no tags updated')

        self.content = result

    def add_query(self, query):
        """
        Add query to be executed when calling `run`
        """
        self.queries.append(query + '\n')

    def add_training_query(
        self, bound, features, goal, exp='minE', strat_name='S', cost=None,
        store_path=None
    ):
        cost = cost or self.cost
        if cost is None:
            raise ValueError('`cost` cannot be None')

        query = 'strategy {} = {} ({}) [{}] {}: <> {}'.format(
            strat_name, exp, cost, bound, features, goal
        )
        self.add_query(query)
        if store_path is not None:
            path = pathlib.Path(store_path)
            self.add_query(f'saveStrategy("{path.resolve()}", {strat_name})')

    def run(self, *args):

        # make temp file for updated model
        m_file = NamedTemporaryFile(mode='w+t', suffix='.xml')
        m_file.write(self.content)
        m_file.flush()

        # make temp file for queries
        q_file = NamedTemporaryFile(mode='w+t')
        q_file.writelines(self.queries)
        q_file.flush()

        eval_args = (
            os.environ['VERIFYTA_PATH'],
            m_file.name, q_file.name,
            *args
        )
        result = subprocess.run(
            eval_args, capture_output=True, encoding='utf-8'
        )

        m_file.close()
        q_file.close()
        return result

    def parse_ouput(self, output):
        if output.stderr:
            raise ValueError('execution failed with: {}'.format(output.stderr))

        pattern = f'(\nVerifying formula \d+ [\s\S]+?)(?=\nVerifying|$)'
        matches = re.findall(pattern, output.stdout)

        assert len(matches) == len(self.queries)

        res_pattern = r'E\((?:max|min)\) = (?:(\d+\.?\d*).{3}(\d+\.?\d*) \(95% CI\)|≈ (\d+\.?\d*))'
        parsed_output = []
        for query, match in zip(self.queries, matches):
            query = query.strip()
            if 'strategy S' in query:
                parsed_output.append((query, 'Formula is satisfied' in match))
                continue

            if 'loadStrat' in query:
                continue

            res_matches = re.search(res_pattern, match).groups()
            if res_matches[:2] == (None, None):
                parsed_output.append((query, (float(res_matches[2]), 0.0)))
            else:
                parsed_output.append((query, tuple(map(float, res_matches[:2]))))

        return parsed_output
