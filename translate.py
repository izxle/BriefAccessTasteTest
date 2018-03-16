from read_data import read_MED_data
from glob import glob
from os import path, mkdir
import shutil
import numpy as np


def write_MED_file(filename, data):
    text = ''
    data = vars(data)
    dimension_keys = [k for k, v in data.items()
                      if isinstance(v, np.ndarray)]

    for k, v in data.items():
        if k in dimension_keys or ('_' in k and k.replace('_', ' ') in data):
            continue
        elif len(k) > 1:
            value = v
        else:
            value = f'{v:11.3f}'

        text += f'{k.capitalize()}: {value}\n'

    for k in sorted(dimension_keys):
        values = data[k]
        dim = ''
        assert len(values), f'filename: {filename}\nvalues is empty for {k}'
        for i in range(int(len(values) / 5) + 1):
            line = ''.join(f'{v:13.3f}'
                           for v in values[i*5: i*5 + 5])
            dim += f'{i*5:6d}:{line}\n'
        text += f'{k.upper()}:\n{dim}'

    with open(filename, 'w') as f:
        f.write(text)


def main():
    working_directory = r'C:\Users\Mona__Luna\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\MedPC\AiresOpenLoop\OpenLoop'
    working_directory = r'C:\Users\izxle\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\MedPC\AiresOpenLoop\BaseLine\close loop'
    relative_path = '*'
    absolute_path = path.join(working_directory, relative_path)
    file_list = glob(absolute_path)
    table = {'d': 'a',
             'm': 'h',
             'a': 'm'}
    old_dir = path.join(working_directory, 'old_files')
    if not path.exists(old_dir):
        mkdir(old_dir)
    elif not path.isdir(old_dir):
        raise NameError(f'old_files already exists and is not a directory')

    for file_path in file_list:
        if not path.isfile(file_path):
            continue
        file_name = path.basename(file_path)
        data = read_MED_data(file_path)
        tmp_data = dict(vars(data))
        for k in table:
            delattr(data, k)

        for old_key, new_key in table.items():
            setattr(data, new_key, tmp_data[old_key])

        new_path = path.join(old_dir, file_name)
        shutil.move(file_path, new_path)
        try:
            write_MED_file(file_path, data)
        except Exception as e:
            shutil.move(new_path, file_path)
            raise e


if __name__ == '__main__':
    main()
