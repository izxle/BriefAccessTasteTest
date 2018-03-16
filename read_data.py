
import re
import numpy as np
from os import path
from scipy.io import loadmat
from collections import OrderedDict


class DataObject():
    def __init__(self, data, table={}):
        for key, value in data.items():
            setattr(self, key, value)
            if ' ' in key:
                setattr(self, key.replace(' ', '_'), value)
            if table and key in table:
                setattr(self, table[key], value)

    def __getitem__(self, item):
        return getattr(self, item, None)

    def __repr__(self):
        return vars(self)['start date']


def parse_data(text, tagged=False):
    data = []
    if tagged:
        tagged_data = {}
        for line in text.split('\n'):
            nums = line.strip().split()[1:]

            for n in nums:
                tag = int(n[-1])
                val = float(n[:-1])
                if tag in tagged_data:
                    tagged_data[tag].append(val)
                else:
                    tagged_data[tag] = [val]
        data = [np.array(tagged_data[i])
                for i in sorted(tagged_data.keys())]
    else:
        for line in text.split('\n'):
            data += list(map(float, line.strip().split()[1:]))
        data = np.array(data, dtype=float)
        # TODO: remove nan values
        data = data[data!=0]
    return  data


info_regex = re.compile(r'^(?P<name>\w{3,}(?: \w+)?): +(?P<data>.+)', re.MULTILINE)
dim_regex = re.compile(r'(?P<name>\w):\s+(?P<data>(?:\d+:[\d. ]+\n *)+)')

def read_MED_data(filename, table={}, tags=[]):

    with open(filename, 'r') as f:
        text = f.read()

    try:

        info = {m.group('name').lower(): m.group('data')
                for m in info_regex.finditer(text)}
        dimensions = {m.group('name').lower(): parse_data(m.group('data'), table.get(m.group('name').lower()) in tags)
                      for m in dim_regex.finditer(text)}

        for k, v in dimensions.items():
            if isinstance(v, str) or isinstance(v, list) or v.any():
                info[k] = v
        # info.update(dimensions)

        data = DataObject(data=info, table=table)
    except Exception as e:
        raise IOError(f'Did not read {filename}\n{e}')

    return data


def read_Matlab_file(file_path: str):
    plx_data = loadmat(file_path)

    spikes_plx = {channel: signal
                  for channel, signal in zip(plx_data['Spikes'].dtype.names,
                                             plx_data['Spikes'][0, 0])}

    events_plx = {event: signal
                  for event, signal in zip(plx_data['Events'].dtype.names,
                                           plx_data['Events'][0, 0])}

    events_plx['Spikes'] = spikes_plx
    filename = path.basename(file_path)
    date = re.search('\d{6}', filename).group()
    events_plx['date'] = date
    events_plx['start_date'] = date
    events_plx['stop_date'] = date
    events_plx['subject'] = filename[:3]

    data = DataObject(data=events_plx)

    return data


def get_bouts(data, dt=0.5):
    n_bouts = 0
    bout = 0

    prev = data[0]
    for s in data[1:]:
        ds = s - prev

        if ds < dt:
            bout += 1
        else:
            if bout > 2:
                n_bouts += 1
            bout = 0
        prev = s
    return n_bouts