
import re
from argparse import ArgumentParser, Namespace
from configparser import ConfigParser
from glob import glob
from os import path, listdir
import numpy as np

import matplotlib.pyplot as plt

from plotData import raster, ordered_raster_PSTH, trials_per_day, open_loop, spike_raster_PSTH, open_loop_color_task
from plotData import colors, background_color
from read_data import read_MED_data, read_Matlab_file, get_bouts


plot_colors = [colors['blue'],
               colors['blue'],
               colors['green'],
               colors['orange'],
               colors['blue'],
               colors['blue'],
               colors['blue'],]


def get_args(argv=''):
    """

    :param argv:
    :return:
    """
    parser = ArgumentParser()
    parser.add_argument('config', default='Configuration.ini', nargs='?')
    parser.add_argument('--date')

    if argv:
        if isinstance(argv, str):
            argv = argv.split()
        elif not hasattr(argv, '__iter__'):
            raise ValueError(f'argv must be string or iterable.')
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()

    return args


def read_config(filename):
    """

    :param filename:
    :return:
    """
    config = ConfigParser(allow_no_value=True)
    assert path.exists(filename), f'{filename} does not exist, must be an existing file'
    config.read(filename)

    opts = Namespace()

    # read GENERAL section
    general = config['GENERAL']
    # extract path to directory with MED data
    opts.path = general['Directory']
    assert path.isdir(opts.path), f'{opts.path} must be an existing directory'
    # extract experiment names
    opts.experiments = {k: {}
                        for k in general['Experiment'].split()}
    # extract mouse ID's
    opts.IDs = general['Mouse ID'].split()
    opts.date = general.get('date', '').split()


    if general.get('tag'):
        opts.tag = re.split(' *, *', general.get('tag', '0'))
    else:
        opts.tag = []
    # extract table as a dictionary
    opts.table = config._sections['TABLE']
    # extract info from experiments' sections
    for exp in opts.experiments:
        # TODO: maybe raise exception if experiment not a section
        section = config._sections.get(exp, {})
        assert section, f'section {exp} not found'
        for key, value in section.items():
            if value is None:
                continue
            elif '\n' in value:
                info = {}
                if re.search('[=:]', value):
                    for line in value.strip().split('\n'):
                        k, v = re.split(' *[:=] *', line)
                        if ',' in v:
                            v = re.split(' *, *', v)
                        elif k == 'plots':
                            v = [v]
                        elif k in ['plot', 'mark_end', 'average_licks']:
                            v = False if v == 'False' else True
                        elif re.match('-?\d+$', v):
                            v = int(v)
                        elif re.match('(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$', v):
                            v = float(v)

                        info[k] = v
                else:
                    info = value
                section[key] = info
            elif ',' in value:
                section[key] = re.split(' *, *', value)
        opts.experiments[exp] = section

    opts.restrictions = {'mouse': opts.IDs,
                         'date': opts.date}

    return opts


def get_file_list(directory, experiments: dict):
    """

    :param config:
    :return:
    """
    file_list = []
    for info in experiments.values():
        dir_name = path.join(directory, info['folder'])
        file_list += glob(path.join(dir_name, '*'))

    return file_list


def group_by_mouse(data_list):
    """
    Group a list of MED_DATA objects in a dictionary with MED_DATA.subject as the key
    :param data_list: List[MED_DATA]
    :return: Dict[str, List[MED_DATA]]
    """
    mouse_data = {}
    for data in data_list:
        if mouse_data.get(data.subject):
            mouse_data[data.subject].append(data)
        else:
            mouse_data[data.subject] = [data]

    for mouse_data_list in mouse_data.values():
        mouse_data_list.sort(key=lambda data: data.start_date)

    return mouse_data


def filter_data(data, restrictions):
    """

    :param data: List[MED_DATA]
    :param restrictions: Dict[str, List[str]]
    :return: List[MED_DATA]
    """
    res = [d for d in data
           if ((not restrictions['date'] or d.start_date in restrictions['date'])
           and (not restrictions['mouse'] or d.subject in restrictions['mouse']))]
    return res


def run(config, MED_DATA_list):
    bound_lick_counter = 0
    color_counter = 0
    plot_total_counter = []
    bout_counter = 0

    for name, info in config.experiments.items():
        dir_name = path.join(config.path, info['folder'])
        experiment_files = listdir(dir_name)

        MED_DATA = [data for data in MED_DATA_list
                    if path.basename(data.file) in experiment_files]
        if not MED_DATA:
            continue
        if 'raster' in info:
            signal_in_trial = []
            for data in MED_DATA:
                res = raster(data, **info['raster'])
                signal_in_trial.append(res / (len(data.Licks) if hasattr(data, 'Licks') else 1))
        if 'raster psth' in info:
            signal_in_trial = []
            for data in MED_DATA:
                res = raster(data, **info['raster psth'])
                signal_in_trial.append(res / (len(data.Licks) if hasattr(data, 'Licks') else 1))
        if 'ordered raster' in info:
            m = len(MED_DATA)
            n = max(len(info['ordered raster'].get('valve', 'Sucrose')) for data in MED_DATA)
            label_x = info['ordered raster']['labels']
            if isinstance(label_x, str):
                label_x = [label_x]
            assert len(label_x) == n, f"{len(label_x)} labels doesn't match the number of signals ({n})"

            signal_per_conc = np.zeros(n, dtype=float)

            for data in MED_DATA:
                res = ordered_raster_PSTH(data, **info['ordered raster'])
                # print(res)
                if len(res) < n:
                    tmp = np.zeros(n, dtype=float)
                    # tmp[:len(res)] = res
                    res = tmp
                signal_per_conc += res
            signal_per_conc /= m
            if info['ordered raster'].get('average_licks'):
                plt.figure('Average rewarded licks')
                plt.title('Average rewarded licks')
                plt.bar(range(n), signal_per_conc, tick_label=label_x, color=background_color)
            # plt.ylim([0, max(signal_per_conc)+5])

        if 'trials_per_day' in info:
            # sort MED_DATA by mouse in a dictionary
            mouse_data = group_by_mouse(MED_DATA)
            trials_per_day(mouse_data)

        if 'close loop' in info:
            mouse_data = group_by_mouse(MED_DATA)
            bound_licks = {}
            for mouse, data_list in mouse_data.items():
                bound_licks[mouse] = {}
                for data in data_list:
                    res = raster(data, **info['close loop'])
                    if hasattr(data, 'Licks'):
                        index = res / len(data.Licks)
                    else:
                        index = 0
                    bound_licks[mouse][data.start_date] = index

            plt.figure(f'Bound Licks')
            for mouse, licks_with_laser in bound_licks.items():
                licks = [licks for date, licks in sorted(licks_with_laser.items())]

                i = bound_lick_counter + 1
                days = len(licks)
                plt.plot(range(i, i + days), licks, marker='o', label=f'{mouse} {name}')
                bound_lick_counter += days

            plt.xlabel('Day')
            plt.ylabel('Licks in laser period / Total licks')
            plt.legend(loc=0)
            plt.title('Bound Licks')

        if 'open loop' in info:
            mouse_data = group_by_mouse(MED_DATA)
            bound_licks = {}
            for mouse, data_list in mouse_data.items():
                bound_licks[mouse] = {}
                for data in data_list:
                    res = open_loop(data, **info['open loop'])
                    if hasattr(data, 'Licks'):
                        index = res / len(data.Licks)
                    else:
                        index = 0
                    bound_licks[mouse][data.start_date] = index

            plt.figure(f'Bound Licks')
            for mouse, licks_with_laser in bound_licks.items():
                licks = [licks for date, licks in sorted(licks_with_laser.items())]
                i = bound_lick_counter + 1
                days = len(licks)
                plt.plot(range(i, i + days), licks, marker='o', label=f'{mouse} {name}')
                bound_lick_counter += days

            plt.xlabel('Day')
            plt.ylabel('Licks in laser period / Total licks')
            plt.legend(loc=0)
            plt.title('Bound Licks')

        if 'open loop color task' in info:
            mouse_data = group_by_mouse(MED_DATA)
            bound_licks = {}
            for mouse, data_list in mouse_data.items():
                bound_licks[mouse] = {}
                for data in data_list:
                    res = open_loop_color_task(data, **info['open loop color task'])
                    if hasattr(data, 'Licks'):
                        index = res / len(data.Licks)
                    else:
                        index = 0
                    bound_licks[mouse][data.start_date] = index

            plt.figure(f'Bound Licks')
            for mouse, licks_with_laser in bound_licks.items():
                licks = [licks for date, licks in sorted(licks_with_laser.items())]
                i = bound_lick_counter + 1
                days = len(licks)
                plt.plot(range(i, i + days), licks, marker='o', label=f'{mouse} {name}')
                bound_lick_counter += days

            plt.xlabel('Day')
            plt.ylabel('Licks in laser period / Total licks')
            plt.legend(loc=0)
            plt.title('Bound Licks')

        if 'spike raster psth' in info:
            for data in MED_DATA:
                spike_raster_PSTH(data, **info['spike raster psth'])

        if 'plot' in info:
            plot = info['plot']
            if isinstance(plot, str):
                plot = [plot]
            for n, signal in enumerate(plot):


                plt.figure(f'Total {signal}')
                plt.title(MED_DATA[-1].subject)
                label = name
                # color = plot_colors[color_counter]

                y = [len(getattr(data, signal)) if hasattr(data, signal) else 0
                     for data in MED_DATA]

                try:
                    i = plot_total_counter[n]
                except IndexError:
                    plot_total_counter.append(0)
                    i = 0
                i += 1
                days = len(y)
                x = range(i, i + days)
                plt.bar(x, y, label=label)
                plt.xlabel('Day')
                plt.ylabel(signal)
                plt.legend()

                # color_counter += 1
                plot_total_counter[n] += days

        if 'bouts' in info:
            signal = info['bouts']

            y = [get_bouts(getattr(data, signal)) if hasattr(data, signal) else 0
                 for data in MED_DATA]

            plt.figure('Bouts')
            plt.title(MED_DATA[-1].subject)

            i = bout_counter + 1
            days = len(y)
            x = range(i, i + days)
            plt.bar(x, y, label=name)
            plt.xlabel('Day')
            plt.ylabel(f'Bouts of {signal}')
            plt.legend()

            bout_counter += days

    plt.show()
    return


def main(argv=''):
    args = get_args(argv)
    config = read_config(args.config)
    # print(config)
    # get list of files in specified folders in each experiment
    file_list = get_file_list(config.path, config.experiments)
    # read and parse files to obtain MED_DATA
    data_files = [(read_Matlab_file(filename)
                   if filename[-4:] == '.mat' else
                   read_MED_data(filename, config.table, config.tag))
                  for filename in file_list
                  if path.isfile(filename)]

    # filter MED_DATA
    MED_DATA = filter_data(data_files, config.restrictions)

    run(config, MED_DATA)
    # MED_data = []
    # for mouse, data in mouse_data.items():
    #     if mouse in config.IDs:
    #         MED_data += data
    #
    # if args.date:
    #     for data in MED_data:
    #         if data.start_date != args.date: continue
    #         raster(data, plots=["Licks", "Sacarose"])


if __name__ == '__main__':
    main('test.ini')
