
import matplotlib.pyplot as plt
import numpy as np
from math import ceil
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

# RGB definition of colors
colors = dict(
    white = [255, 255, 255, 255],
    black = [  0,   0,   0, 255],
    red   = [255,   0,   0, 255],
    green = [  0, 128,   0, 255],
    blue  = [  0,   0, 255, 255],
    lime  = [  0, 255,   0, 255],
    yellow = [255, 255, 0, 255],
    orange = [255, 165, 0, 255],
    crimson = [220, 60, 60, 255],
    silver = [192, 192, 192, 255],
    dodgerblue = [30, 144, 255, 255],
    deeppink = [255, 20, 147, 255],
    hotpink = [255, 105, 180, 255],
    indigo = [75, 0, 130, 255],
    lavender = [230, 230, 255, 255])
colors = {k: [v / 255 for v in vals]
          for k, vals in colors.items()}

background_color = [colors['blue'],
                    colors['dodgerblue'],
                    colors['yellow'],
                    colors['orange'],
                    colors['hotpink'],
                    colors['crimson']]


def _ordered_raster_PSTH(signal, valves, start_cue, stop_cue, ax_raster, ax_PSTH,
                         resolution=0.01, trial_time=4, sub_trial_time=0, extra_time_after=3, extra_time_before=1,
                         line_width=3, reference='stop_cue', bg_color=background_color,
                         plot=True, bin_size=0.05, labels=None, **kwargs):
    if reference=='stop_cue':
        ref = stop_cue[:]
    elif reference=='start_cue':
        ref = start_cue[:]
    else:
        raise ValueError('reference must be either `stop_cue` or `start_cue`')


    # initialize raster dimensions
    window = extra_time_before + trial_time + extra_time_after
    m = len(stop_cue)
    n = int(window / resolution) + line_width
    # define an all-white image. Each cell has RGB with (1, 1, 1) as white
    img = np.ones((m, n, 4), dtype=float)
    if labels is None:
        labels = list(range(len(valves)))
    # compute indices for start and end of trial relative to each job
    j_start = int(extra_time_before / resolution)
    j_end = int((extra_time_before + trial_time) / resolution)
    trial_indices = np.arange(j_start, j_end, dtype=int)
    if sub_trial_time:
        j_sub_end = int((extra_time_before + trial_time - sub_trial_time) / resolution)
        sub_trial_indices = np.arange(j_start, j_sub_end, dtype=int)


    i = 0
    signal_per_conc = np.zeros(len(valves))
    for conc, arr in enumerate(valves):
        j = 0
        data_index = []
        trial_num = []

        for trial_index, s in enumerate(ref):
            if i + j >= m:
                break
            # range for licks inside window
            if reference == 'start_cue':
                trial_start = s
                trial_end = s + trial_time
            else:
                trial_start = s - trial_time
                trial_end = s

            if not arr[(trial_start <= arr) & (arr <= trial_end)].any():
                continue

            window_start = trial_start - extra_time_before
            window_end = trial_end + extra_time_after

            if reference == 'start_cue':
                ref_values = start_cue[(window_start <= start_cue) & (start_cue <= window_end)]
                other_values = stop_cue[(window_start <= stop_cue) & (stop_cue <= window_end)]
            else:
                ref_values = stop_cue[(window_start <= stop_cue) & (stop_cue <= window_end)]
                other_values = start_cue[(window_start <= start_cue) & (start_cue <= window_end)]

            ref_ix = np.array((ref_values - window_start) / resolution, dtype=int)
            other_cue_ix = np.array((other_values - window_start) / resolution, dtype=int)

            values = signal[(window_start <= signal) & (signal < window_end - (line_width*resolution))]
            value_indices = np.array((values - window_start) / resolution, dtype=int)
            signal_per_conc[conc] += sum((trial_start <= signal) & (signal <= trial_end))

            data_index += value_indices.tolist()
            # colors
            img[i + j, trial_indices] = bg_color[conc]
            if sub_trial_time:
                img[i + j, sub_trial_indices] = bg_color[conc] / np.array([1,1,1,2], dtype=float)
            for w in range(line_width):
                img[i + j, value_indices + w] = colors['black']
                img[i + j, ref_ix + w] = colors['silver']
                img[i + j, other_cue_ix + w] = colors['green']

            j += 1
            trial_num.append(trial_index)
        i += j

        # for t in range(j_start, psth_end):
        #     upper_index = t + d_PSTH if t + d_PSTH < psth_end - 1 else psth_end - 1
        #     PSTH_data[t] = n_data[t - d_PSTH: upper_index].sum()
        # PSTH_data = (PSTH_data * resolution) - extra_time_before
        # ax_PSTH.plot(PSTH_data, color=bg_color[conc], linewidth=line_width)
        # x = np.arange(-extra_time_before, trial_time + extra_time_after, resolution)
        if plot:
            relative_bin_size = int(bin_size / resolution)
            n_bins = int((n - line_width) / relative_bin_size)

            # edges = np.arange(0, (n - line_width), relative_bin_size)
            # hist, bin_edges = np.histogram(data_index, edges, range=(0, n - line_width))
            hist, bin_edges = np.histogram(data_index, n_bins, range=(0, n-line_width))
            psth = hist / (m * bin_size)
            w_len = relative_bin_size if relative_bin_size % 2 else relative_bin_size + 1
            if w_len <= 3:
                w_len = 5
            y = savgol_filter(psth, w_len, 3)
            y[y<0] = 0
            x = np.linspace(-extra_time_before, trial_time + extra_time_after, len(y))
            y[x<0] = 0
            # if conc == 4:
            #     ax_raster.plot(x, psth, label=labels[conc], color=bg_color[conc])
            ax_PSTH.plot(x, y, label=labels[conc], color=bg_color[conc])

    def format_func(value, tick_number):
        return int((value * resolution) - extra_time_before)

    if plot:
        ax_PSTH.set_xlim(-extra_time_before, trial_time + extra_time_after)
        ax_raster.imshow(img, aspect='auto')
        ax_raster.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        plt.legend()

    return signal_per_conc


def ordered_raster_PSTH(data, valve='Sucrose', reference='start_cue', plot=True, **kw):
    start_cue = data.Start_Trial
    stop_cue = data.Stop_Trial

    signal = data.Licks

    valves = getattr(data, valve)

    if not isinstance(valves, list):
        valves = [valves]

    name = f'{data.subject} {data.start_date}'
    if plot:
        fig = plt.figure(name)
        ax_raster = fig.add_subplot(211)
        ax_raster.set_title(name)
        ax_PSTH = fig.add_subplot(212)
    else:
        ax_raster = None
        ax_PSTH = None


    signal_per_conc =_ordered_raster_PSTH(signal=signal, valves=valves, start_cue=start_cue, stop_cue=stop_cue,
                                          ax_raster=ax_raster, ax_PSTH=ax_PSTH, reference=reference, plot=plot, **kw)
    return signal_per_conc

def gen_colors(color_names=[], default='black'):
    # print('init gen_colors')
    for color in color_names:
        # print(f'  yield color {color}')
        yield colors[color]
    while True:
        # print(f'  yield default {default}')
        yield colors[default]

def _raster2(signal, start_cue=None, stop_cue=None, resolution=0.01, trial_time=4,
             extra_time_after=3, extra_time_before=1, line_width=3, ax=None, reference='stop_cue',
             PSTH=False, ax_PSTH=None, bin_size=0.05,
             other=None, color_names=[], labels=[], end_color='green', mark_end=True, plot=True):
    signal_in_trial = 0
    if start_cue is None and stop_cue is None:
        raise ValueError('at least one of `stop_cue` or `start_cue` must be given')

    if other is None:
        other = []

    if reference=='stop_cue':
        ref = stop_cue
        # if start_cue is not None:
        #     other.append(start_cue)
    elif reference=='start_cue':
        ref = start_cue
        # if stop_cue is not None:
        #     other.append(stop_cue)
    else:
        raise ValueError('reference must be either `stop_cue` or `start_cue`')

    # initialize raster dimensions
    window = extra_time_before + trial_time + extra_time_after
    m = len(ref)
    n = int(window / resolution) + line_width
    # define an all-white image. Each cell has RGB with (1, 1, 1) as white
    img = np.ones((m, n, 4), dtype=float)

    # compute indices for start and end of trial relative to each job
    j_start = int(extra_time_before / resolution)
    j_end = int((extra_time_before + trial_time) / resolution)
    trial_indices = np.arange(j_start, j_end, dtype=int)

    signal_indices = []
    other_indices = [[] for _ in other]

    other_ix = [[] for _ in other]
    data_ix = []

    # iterate over all reference signals
    for i, s in enumerate(ref):
        # color each section
        img[i, trial_indices] = colors['lavender']
        if mark_end:
            img[i, j_start] = colors['green']
            img[i, j_end] = colors['green']
        # range for licks inside window
        if reference=='start_cue':
            start = s - extra_time_before
            end = s + trial_time + extra_time_after
        else:
            start = s - trial_time - extra_time_before
            end = s + extra_time_after

        for o, o_values in enumerate(other):
            if len(o_values) == 0:
                other_indices[o].append([])
                continue
            o_vals = o_values[(start <= o_values) & (o_values <= end)]
            o_indices = np.array((o_vals - start) / resolution, dtype=int)
            other_indices[o].append(o_indices)
            other_ix[o] += o_indices.tolist()

        if isinstance(signal, list):
            signal_ix = []
            for i, arr in enumerate(signal):
                values = arr[(start <= arr) & (arr <= end)]
                signal_ix += list(((values - start) / resolution))
            signal_ix = np.array(signal_ix, dtype=int)
            # TODO: add data_ix
        elif isinstance(signal, np.ndarray):
            values = signal[(start <= signal) & (signal <= end)]
            signal_ix = np.array(((values - start) / resolution), dtype=int)
            data_ix += signal_ix.tolist()
        else:
            raise TypeError(f'signal must be a list or ndarray, not {type(signal)}')

        signal_in_trial += ((j_start <= signal_ix) & (signal_ix < j_end)).sum()

        signal_indices.append(signal_ix)

    # set color for licks
    signal_color = next(gen_colors(color_names))
    for i in range(len(img)):
        for w in range(line_width):
            img[i, signal_indices[i] + w] = signal_color

    color = gen_colors(color_names[1:])
    other_colors = [next(color) for _ in other]
    for i in range(len(img)):
        for w in range(line_width):
            color = gen_colors(color_names[1:])
            if mark_end:
                img[i, j_start + w] = colors[end_color]
                img[i, j_end + w] = colors[end_color]
            for o, _ in enumerate(other):
                if len(other_indices[o][i]) > 0:
                    img[i, other_indices[o][i] + w] = other_colors[o]

    def format_func(value, tick_number):
        return int((value * resolution) - extra_time_before)

    if plot:
        if ax_PSTH:
            relative_bin_size = int(bin_size / resolution)
            n_bins = int((n - line_width) / relative_bin_size)

            # edges = np.arange(0, (n - line_width), relative_bin_size)
            # hist, bin_edges = np.histogram(data_index, edges, range=(0, n - line_width))
            hist, bin_edges = np.histogram(data_ix, n_bins, range=(0, n - line_width))
            psth = hist / (m * bin_size)
            w_len = relative_bin_size if relative_bin_size % 2 else relative_bin_size + 1
            if w_len <= 3:
                w_len = 5
            y = savgol_filter(psth, w_len, 3)
            y[y < 0] = 0
            x = np.linspace(-extra_time_before, trial_time + extra_time_after, len(y))
            # y[x < 0] = 0
            # if conc == 4:
            #     ax_raster.plot(x, psth, label=labels[conc], color=bg_color[conc])
            label = iter(labels)
            color = gen_colors(color_names)
            ax_PSTH.plot(x, y, label=next(label), color=next(color))

            for o, _ in enumerate(other):
                hist, bin_edges = np.histogram(other_ix[o], n_bins, range=(0, n - line_width))
                psth = hist / (m * bin_size)
                w_len = relative_bin_size if relative_bin_size % 2 else relative_bin_size + 1
                if w_len <= 3:
                    w_len = 5
                y = savgol_filter(psth, w_len, 3)
                y[y < 0] = 0
                x = np.linspace(-extra_time_before, trial_time + extra_time_after, len(y))
                ax_PSTH.plot(x, y, color=next(color), label=next(label, f'other{o}'))
                ax_PSTH.set_xlim([-extra_time_before, trial_time + extra_time_after])

        if ax:
            ax.imshow(img, aspect='auto')
            ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
            # manage labels and colors
            label = iter(labels)
            color = gen_colors(color_names)
            ax.plot([], [], color=next(color), label=next(label, 'Signal'))
            for o, _ in enumerate(other):
                ax.plot([], [], color=next(color), label=next(label, f'other{o}'))
            ax.legend()

        else:
            # TODO: integrate
            plt.imshow(img, aspect='auto')
            plt.legend()

    return signal_in_trial


def _raster(signal, start_cue=None, stop_cue=None, resolution=0.01, trial_time=4, bg_color=colors['lavender'],
            extra_time_after=3, extra_time_before=1, line_width=3, ax=None, reference='stop_cue',
            mark_end=True, plot=True):
    signal_in_trial = 0
    if start_cue is None and stop_cue is None:
        raise ValueError('at least one of `stop_cue` or `start_cue` must be given')


    if reference == 'stop_cue':
        ref = stop_cue
    elif reference == 'start_cue':
        ref = start_cue
    else:
        raise ValueError('reference must be either `stop_cue` or `start_cue`')

    # initialize raster dimensions
    window = extra_time_before + trial_time + extra_time_after
    m = len(ref)
    n = int(window / resolution) + line_width
    # define an all-white image. Each cell has RGB with (1, 1, 1) as white
    img = np.ones((m, n, 4), dtype=float)

    # compute indices for start and end of trial relative to each job
    j_start = int(extra_time_before / resolution)
    j_end = int((extra_time_before + trial_time) / resolution)
    trial_indices = np.arange(j_start, j_end, dtype=int)

    signal_indices = []
    other_cue_indices = []

    # iterate over all reference signals
    for i, s in enumerate(ref):
        # color each section
        img[i, trial_indices] = colors['lavender']
        if mark_end:
            img[i, j_start] = colors['green']
            img[i, j_end] = colors['green']
        # range for licks inside window
        if reference == 'start_cue':
            start = s - extra_time_before
            end = s + trial_time + extra_time_after
            if stop_cue is not None:
                stop_cue_values = stop_cue[(start <= stop_cue) & (stop_cue <= end)]
                other_cue_ix = np.array((stop_cue_values - start) / resolution, dtype=int)
        else:
            start = s - trial_time - extra_time_before
            end = s + extra_time_after
            if start_cue is not None:
                start_cue_values = start_cue[(start <= start_cue) & (start_cue <= end)]
                other_cue_ix = np.array((start_cue_values - start) / resolution, dtype=int)

        if isinstance(signal, list):
            signal_ix = []
            for i, arr in enumerate(signal):
                values = arr[(start <= arr) & (arr <= end)]
                signal_ix += list(((values - start) / resolution))
            signal_ix = np.array(signal_ix, dtype=int)
        elif isinstance(signal, np.ndarray):
            values = signal[(start <= signal) & (signal <= end)]
            signal_ix = np.array(((values - start) / resolution), dtype=int)
        else:
            raise TypeError(f'signal must be a list or ndarray, not {type(signal)}')

        signal_in_trial += ((j_start <= signal_ix) & (signal_ix < j_end)).sum()

        signal_indices.append(signal_ix)
        if start_cue is not None and stop_cue is not None:
            other_cue_indices.append(other_cue_ix)

    for i in range(len(img)):
        # set color for licks
        for w in range(line_width):
            img[i, signal_indices[i] + w] = colors['black']
            img[i, j_start + w] = colors['green']
            if mark_end:
                img[i, j_end + w] = colors['green']
            if start_cue is not None and stop_cue is not None:
                img[i, other_cue_indices[i] + w] = colors['red']

    def format_func(value, tick_number):
        return int((value * resolution) - extra_time_before)

    if plot:
        if ax:
            ax.imshow(img, aspect='auto')
            ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
            ax.plot([], [], color=colors['red'], label='Head entry')
            ax.plot([], [], color=colors['black'], label='Licks')
            ax.legend()

        else:
            plt.imshow(img, aspect='auto')
            plt.legend()

    return signal_in_trial


def raster(data, plots=[], reference='stop_cue', start_cue='', stop_cue='', plot=True, labels=[],
           other=[], PSTH=False, **kw):
    # extract and rename wanted data
    if isinstance(labels, str):
        labels = [labels]
    if isinstance(other, str):
        other = [other]
    if reference == 'start_cue':
        assert start_cue, 'start_cue must be provided if it is to be the reference.'
        label_reference = start_cue
        start_cue = data[start_cue]
        if stop_cue:
            if stop_cue not in other:
                other.append(stop_cue)
                labels.append(stop_cue)
            stop_cue = None
        
    elif reference == 'stop_cue':
        assert stop_cue, 'stop_cue must be provided if it is to be the reference.'
        label_reference = stop_cue
        stop_cue = data[stop_cue]
        if start_cue:
            if start_cue not in other:
                other.append(start_cue)
                labels.append(start_cue)
            start_cue = None
        
    else:
        raise ValueError('reference must be either `stop_cue` or `start_cue`')

    
    if plots:
        signals = []
        for p in plots:
            if hasattr(data, p):
                signal = getattr(data, p)
            else:
                signal = []
            signals.append(signal)

        # signals = [getattr(data, p) if hasattr(data, p) else []
        #            for p in plots]
    else:
        signals = [data.Licks]

    if not labels:
        labels = list(other)
    other = [getattr(data, o) if hasattr(data, o) else [] for o in other]

    n = len(signals)
    title = f'{data.subject} {data.start_date} \nAligned {label_reference}'
    if PSTH:
        n += 1
        title += ' PSTH'
        if isinstance(PSTH, str):
            PSTH = [PSTH]

    name = title

    if plot:
        # fig = plt.figure(name)
        # axes = fig.add_subplot(n, 1)
        fig, axes = plt.subplots(n, 1, num=name)
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        axes[0].set_title(title)
    else:
        axes = np.zeros(n)

    # print([(signal.shape, len(signal.shape)) for signal in signals])
    for i, signal in enumerate(signals):
        if PSTH:
            ax_PSTH = axes[-1]
        else:
            ax_PSTH = None
        res = _raster2(signal, start_cue, stop_cue, ax=axes[i], reference=reference,
                       plot=plot, other=other, labels=labels, PSTH=PSTH, ax_PSTH=ax_PSTH, **kw)
        if i == 0:
            signal_in_trial = res

    return signal_in_trial


def open_loop(data, trial_time, block_lenght, time_before, time_after, resolution=0.01, line_width=1,
              reference=None, total_time=0, plot=True, **kwargs):

    signal = data.Licks
    if hasattr(data, 'Blocks'):
        blocks = data.Blocks
    elif total_time:
        blocks = list(range(block_lenght, total_time+1, block_lenght))
    else:
        raise ValueError(f'Must specify blocks or total session time.')
    if reference:
        reference = getattr(data, reference)
    else:
        reference = data.LaserEnd

    window = time_before + trial_time + time_after
    n_trials = block_lenght // (trial_time + time_after)

    m = 2 * len(blocks) * n_trials
    n = int(window / resolution) + line_width

    img = np.ones((m, n, 4))
    licks_with_laser = 0

    j_start = int(time_before / resolution)
    j_end = int((time_before + trial_time) / resolution)
    trial_indices = np.arange(j_start, j_end, dtype=int)

    for i, block_start in enumerate(np.arange(0, 2*len(blocks) * block_lenght, block_lenght)):
        block_window = (block_start <= reference) & (reference < block_start + block_lenght)
        laser = block_window.any()
        if laser:
            refs = reference[block_window]
        else:
            init = block_start + time_before + trial_time
            fin = block_start + block_lenght
            step = trial_time + time_after
            refs = np.arange(init, fin+1, step)
        if not refs.any():
            print(f'block {i} has no values!')
        for j, ref in enumerate(refs):
            start = ref - time_before - trial_time
            end = ref + time_after

            values = signal[(start <= signal) & (signal < end - (line_width*resolution))]
            indices = np.array((values - start) / resolution, dtype=int)
            # indices = indices[indices < n - 2*line_width]

            # trial background color
            if laser:
                img[i*n_trials + j, trial_indices] = colors['lavender']
                licks_with_laser += sum((start + time_before <= values) & (values < end - time_after))
            # lick color
            for w in range(line_width):
                img[i*n_trials + j, indices + w] = colors['black']

    if plot:
        def xticks(value, tick_number):
            return int((value * resolution) - time_before)


        def yticks(value, tick_number):
            return int((value / n_trials) * (block_lenght / 60))

        f = plt.figure(f'{data.subject} {data.start_date} - OpenLoop')
        plt.title(f'{data.subject} {data.start_date}')
        ax = f.add_subplot(111)
        ax.imshow(img, aspect='auto')
        ax.set_xlabel('Relative time / s')
        ax.set_ylabel('Time / min')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(xticks))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(yticks))

    # print(f'licks with laser = {licks_with_laser}')

    return licks_with_laser


def open_loop_color_task(data, signal, context, trial_time, time_before, time_after, block_lenght, total_time=0, event_time=0,
                         reference=None, aligned='start', task_reference=None, task_aligned='start', task_time=0,
                         other=None, resolution=0.01, line_width=1, event_bg_color=colors['lavender'], plot=True):
    # checking arguments
    if not total_time:
        raise ValueError(f'Must specify total session time.')
    if not aligned in ['start', 'end']:
        raise ValueError(f'`aligned` must be one of ["start", "end"], not {aligned}.')

    signal = getattr(data, signal, None)
    context = getattr(data, context, None)
    if signal is None:
        raise ValueError(f'Signal "{signal}" not found in data.')
    if context is None:
        raise ValueError(f'Signal "{context}" not found in data.')
    if not isinstance(context, list):
        raise ValueError(f'Signal "{context}" not a list')
    context_size = sum(len(s) for s in context)
    # print(len(signal), context_size)
    if len(signal) < context_size:
        raise ValueError(f'signal <{len(signal)}> must be greater then or equal to context <{context_size}>')

    # # create a correspondence vector
    # correspondence = np.zeros(len(signal)) - 1
    # for conc, vals in enumerate(context):
    #     mask = np.in1d(signal, vals)
    #     correspondence[mask] = conc
    #     assert len(vals) == mask.sum(), f'not all values found for conc {i}'

    blocks = np.arange(0, total_time, block_lenght)

    if reference:
        reference = getattr(data, reference)
    else:
        reference = data.LaserEnd
        
    if task_reference:
        task_reference = getattr(data, task_reference)
    else:
        raise ValueError('`task_reference` was not provided.')

    window = time_before + trial_time + time_after
    # number of trials in each block
    n_trials = block_lenght // trial_time

    m = n_trials * total_time // block_lenght
    n = int(window / resolution) + line_width

    img = np.ones((m, n, 4))
    licks_in_event = 0

    # for event background shading
    j_start = int(time_before / resolution)
    j_end = int((time_before + event_time) / resolution)
    trial_indices = np.arange(j_start, j_end, dtype=int)

    d_task_indices = int(task_time / resolution)

    for i, block_start in enumerate(blocks):
        block_window = (block_start <= reference) & (reference < block_start + block_lenght)
        # task_refs = task_reference[(block_start <= task_reference) & (task_reference < block_start + block_lenght)]
        event = block_window.any()
        if event:
            # get references in event window
            refs = reference[block_window]
        else:
            # construct virtual references
            init = block_start
            fin = block_start + block_lenght
            step = trial_time
            if aligned == 'end':
                init += trial_time
                fin += trial_time
            refs = np.arange(init, fin, step)

        for j, ref in enumerate(refs):
            start = ref - time_before
            end = ref + time_after
            if aligned == 'start':
                end += trial_time
            else:
                start -= trial_time
            # getting values of signal inside trial window
            values = signal[(start <= signal) & (signal < end)]
            relative_values = values - start
            indices = np.array(relative_values / resolution, dtype=int)

            task_refs = task_reference[(start - task_time <= task_reference) & (task_reference < end)]
            # relative_values = task_refs - start

            # trial background color
            if event:
                img[i*n_trials + j, trial_indices] = event_bg_color
                # getting number of values inside trial
                licks_in_event += ((start + time_before <= values) & (values < end - time_after)).sum()
            # lick color
            for tr in task_refs:
                rel_ref = tr - start
                rel_start = int(rel_ref / resolution)
                if rel_start < 0:
                    rel_start = 0
                rel_end = int((rel_ref + task_time) / resolution)
                if rel_end >= n:
                    rel_end = n - 1
                task_indices = np.arange(rel_start, rel_end, dtype=int)
                # for conc, vals in enumerate(context):
                #     if np.in1d(tr, vals).any():
                #         break
                conc = next(conc for conc, vals in enumerate(context) if np.in1d(tr, vals).any())
                img[i * n_trials + j, task_indices] = background_color[conc]

            # abs_indices = np.argwhere(values == signal)
            # color_indices = context[abs_indices]
            # signal_colors = [colors[ci] for ci in color_indices]
            for w in range(line_width):
                img[i*n_trials + j, indices + w] = colors['black']

    if plot:
        # define x ticks as int relative to trial start in seconds
        def xticks(value, tick_number):
            return int((value * resolution) - time_before)
        # define y ticks as int in minutes
        def yticks(value, tick_number):
            return int((value / n_trials) * (block_lenght / 60))

        f = plt.figure(f'{data.subject} {data.start_date} - OpenLoop')
        plt.title(f'{data.subject} {data.start_date}')
        ax = f.add_subplot(111)
        ax.imshow(img, aspect='auto')
        ax.set_xlabel('Relative time / s')
        ax.set_ylabel('Time / min')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(xticks))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(yticks))

    # print(f'licks with laser = {licks_with_laser}')

    return licks_in_event


def trials_per_day(mouse_data):
    # longer = {'date': [], 'trials': []}
    plt.figure('Trials per day')

    for mouse_ID, data_list in mouse_data.items():
        date_trials = [(data.end_date, len(data.Stop_Trial)) for data in data_list]
        trials = [trials for date, trials in sorted(date_trials)]
        plt.plot(range(1, len(trials)+1), trials, marker='o', label=mouse_ID)
        # y = [len(data.Stop_Trial) for data in data_list]
        # n = len(data_list)
        # if len(longer['date']) < n:
        #     longer['date'] = [data.date for data in data_list]
        #     longer['trials'] = range(n)
    # plt.plot(y, label=mouse_ID)
    # plt.xticks(longer['trials'], longer['date'], rotation='vertical')

    plt.xlabel('Day')
    plt.ylabel('Trials')
    plt.legend(loc=0)
    plt.title('Trials per day')
    # plt.show()


def _spike_raster_PSTH(signal, reference, ax_raster, ax_PSTH, resolution=0.01, trial_time=4, bg_color=colors['lavender'],
                       extra_time_after=3, extra_time_before=1, line_width=3, bin_size=0.05):
    window = extra_time_before + trial_time + extra_time_after
    m = len(reference)
    n = int(window / resolution) + line_width
    # define an all-white image. Each cell has RGB with (1, 1, 1) as white
    img = np.ones((m, n, 4
                   ), dtype=float)

    j_start = int(extra_time_before / resolution)
    j_end = int((extra_time_before + trial_time) / resolution)

    n_data = np.zeros(n - line_width)
    spike_index = []
    for i, first_laser in enumerate(reference):
        laser_start = first_laser
        laser_end = first_laser + trial_time

        window_start = laser_start - extra_time_before
        window_end = laser_end + extra_time_after

        values = signal[(window_start <= signal) & (signal < window_end - (line_width*resolution))]
        value_indices = np.array((values - window_start) / resolution, dtype=int)
        n_data[value_indices] += 1
        spike_index += value_indices.tolist()
        img[i, j_start:j_end + 1] = bg_color
        for w in range(line_width):
            img[i, value_indices + w] = colors['black']

    def format_func(value, tick_number):
        return (value * resolution) - extra_time_before

    ax_raster.imshow(img, aspect='auto')
    ax_raster.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    # PSTH
    relative_bin_size = int(bin_size / resolution)
    n_bins = int((n - line_width) / relative_bin_size)

    hist, bin_edges = np.histogram(spike_index, n_bins, range=(0, n-line_width))
    psth = hist / (m * bin_size)
    w_len = relative_bin_size if relative_bin_size % 2 else relative_bin_size + 1
    y = savgol_filter(psth, w_len, 3)
    y[y < 0] = 0
    x = np.linspace(-extra_time_before, trial_time + extra_time_after, len(y))
    ax_PSTH.plot(x, y, linewidth=line_width)
    ax_PSTH.set_xlim(-extra_time_before, trial_time + extra_time_after)


def spike_raster_PSTH(data, trial_time, **kwargs):
    channels = data.Spikes
    if hasattr(data, 'First_Laser'):
        reference = data.First_Laser
    else:
        laser = data.Laser
        laser_difer = np.diff(laser)
        reference = laser[laser_difer > trial_time]
        np.insert(reference, 0, laser[0])


    for i, (unit, signal) in enumerate(channels.items()):
        # if i > 4: break
        fig = plt.figure(f'{data.subject} {data.start_date} {unit}')
        ax_raster = fig.add_subplot(211)
        ax_PSTH = fig.add_subplot(212)

        _spike_raster_PSTH(signal, reference, ax_raster, ax_PSTH, trial_time=trial_time, **kwargs)


def plot_total(MED_DATA, signal, color=''):
    plt.figure('Total')

    y = [getattr(data, signal) for data in MED_DATA]

    plt.plot(y, color=color)


if __name__ == '__main__':
    print("inside plotData.py")