import argparse
import sys
import collections
import os
from operator import add

_stats_path = ""
_software_flag = ""
_graph_flag = ""
_base_path = ""
_graph_title = ""
_output_path = ""
_version = 0.0
plt = None

BOX_SEQ_LEN_LABEL = "Sequence Length"

# Ensure these flags match with main EnTAP project #
# Frame Selection   = 1
# Expression        = 2
# Similarity Search = 3
# Ontology          = 4
#***************************************************#


def init_argparse():
    global _base_path
    global _stats_path
    global _software_flag
    global _graph_flag
    global _graph_title
    global _output_path
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='stats', help='Path to graphing file', type=str)
    parser.add_argument('-s', action='store', dest='soft', help='Software flag', type=int)
    parser.add_argument('-g', action='store', dest='graph', help='Graph flag', type=int)
    parser.add_argument('-t', action='store', dest='title', help='Graph title', type=str)
    parser.add_argument('-p', action='store', dest='path', help='Output path', type=str)
    args = parser.parse_args()
    _stats_path = args.stats
    _output_path = args.path
    _graph_title = args.title
    if _graph_title is not None:
        _graph_title = _graph_title.replace("_", " ")
    _graph_flag = args.graph
    _software_flag = args.soft


def init_version():
    global _version
    _version = sys.version_info[0] + 0.1*sys.version_info[1]


def verify_package():
    is_2 = 2 <= _version < 3
    has_matplot = False
    if is_2:
        import imp
        try:
            imp.find_module("matplotlib")
            has_matplot = True
        except ImportError:
            has_matplot = False
    elif _version >= 3.4:
        import importlib.util
        check = importlib.util.find_spec("matplotlib")
        has_matplot = check is not None
    elif not is_2 and _version < 3.4:
        import importlib
        spam_loader = importlib.find_loader("matplotlib")
        has_matplot = spam_loader is not None
    else:
        print("Python version not compatable.")
        exit(1)
    if has_matplot:
        global plt
        import matplotlib.pyplot as plt
    else:
        print("Matplotlib module not found. Not able to graph data.")
        exit(1)


def init_graphs():
    pass


def frame_selection_graphs():
    if _graph_flag == 1:            # Pie chart of each gene type
        InputValues = parse_input(_stats_path)
        create_pie(_graph_title, _output_path, InputValues.values, InputValues.labels)
    elif _graph_flag == 2:          # Box plot of rejected vs kept sequences and length
        input_vals = parse_input_dict(_stats_path)
        create_boxplot(_graph_title, _output_path, input_vals, BOX_SEQ_LEN_LABEL)

    pass


def expression_graphs():
    if _graph_flag == 1:
        input_vals = parse_input_dict(_stats_path)
        create_boxplot(_graph_title, _output_path, input_vals, BOX_SEQ_LEN_LABEL)
    pass


# Flag 3
def sim_search_graphs():
    if _graph_flag == 1:
        InputValues = parse_input(_stats_path)
        create_bar(_graph_title, _output_path, InputValues.values, InputValues.labels,
                   InputValues.xlabel, InputValues.ylabel)
    elif _graph_flag == 2:  # Stacked informative/uninformative bar graph
        InputValues = parse_sim_stack(_stats_path)
        create_bar_stacked(_graph_title, _output_path, InputValues.label_map,
                           InputValues.xlabel, InputValues.ylabel)

    pass


def ontology_graphs():
    if _graph_flag == 1:
        InputValues = parse_input(_stats_path)
        create_bar(_graph_title, _output_path, InputValues.values, InputValues.labels,
                   InputValues.xlabel, InputValues.ylabel)
    pass


def create_graphs(flag):
    if flag == -1:                  # Initial test case
        exit(0)
    if flag == 1:
        frame_selection_graphs()
    elif flag == 2:
        expression_graphs()
    elif flag == 3:
        sim_search_graphs()
    elif flag == 4:
        ontology_graphs()
    else:
        init_graphs()


# dict[dict] structure
def create_bar_stacked(title, file, value_map, xlab, ylab):
    indices = range(len(value_map.keys()))
    x_labels = []
    legend_labels = []
    stacked_vals = {}
    # Frame or just a temporary flag (if no frame was sent)
    for key in value_map.keys():
        x_labels.append(key)
        for val in value_map[key].keys():
            if not stacked_vals.__contains__(val):
                stacked_vals[val] = []
            stacked_vals[val].append(value_map[key][val])

    plts = []
    prev_key = ""
    totals = [0] * len(value_map.keys())
    for key in stacked_vals:
        if not prev_key == "":
            p = plt.bar(indices, stacked_vals[key], bottom=totals)
        else:
            p = plt.bar(indices, stacked_vals[key])
        plts.append(p)
        legend_labels.append(key)
        prev_key = key
        totals = map(add, totals, stacked_vals[key])

    plt.legend(plts, legend_labels)
    plt.xticks(indices, x_labels)
    plt.ylabel(ylab)
    plt.title(title)
    plt.savefig(file, bbox_inches="tight")


def create_pie(title, file, vals, labels):
    plt.pie(vals, labels=labels, autopct=autopct_vals(vals))
    plt.axis('equal')
    plt.title(title)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(file)


# label_vals = dictionary of lists containing different series. Such as all rejected seqs and kept
def create_boxplot(title, file, label_vals, y_label):
    data = []
    labels = label_vals.keys()
    for key in label_vals.keys():
        temp = []
        for val in label_vals[key]:
            temp.append(val)
        data.append(temp)
    plt.ylabel(y_label)
    plt.boxplot(data, labels=labels)
    plt.title(title)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(file)


def create_bar(title, file, vals, labels, xlabel, ylabel):
    plt.barh(range(len(labels)), vals, align='center', alpha=0.5)   # Horizontal bar graph
    plt.yticks(range(len(labels)), labels)
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(title)
    plt.savefig(file, bbox_inches="tight")


def autopct_vals(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct, v=val)
    return my_autopct


# Return collection containing values and labels for them
# Primarily used for single label - value series, NOT multiple values to one label
def parse_input(path):
    UserVals = collections.namedtuple('Values', ['values', 'labels', 'xlabel', 'ylabel'])
    file = open(path, 'r')
    labels = file.readline().split('\t')
    UserVals.xlabel = labels[0]
    UserVals.ylabel = labels[1]
    UserVals.labels = []
    UserVals.values= []
    for line in file:
        if not line:
            continue
        values = line.strip().split('\t')
        UserVals.labels.append(values[0])
        UserVals.values.append(int(values[1]))
    return UserVals


# Used for box plot with multiple labels to a value
def parse_input_dict(path):
    file = open(path, 'r')
    rows = file.readline().split('\t')  # not used currently
    series_dict = {}
    for line in file:
        if not line:
            continue
        values = line.strip().split('\t')
        if not series_dict.__contains__(values[0]):
            series_dict[values[0]] = []
        series_dict[values[0]].append(int(values[1]))
    return series_dict


#Always has 3 columns
def parse_sim_stack(path):
    UserVals = collections.namedtuple('Values', ['label_map', 'xlabel', 'ylabel'])
    file = open(path,'r')
    rows = file.readline().split('\t')
    UserVals.xlabel = rows[0]
    UserVals.ylabel = rows[1]
    UserVals.label_map = {}

    for line in file:
        if not line or line == '\n':
            continue
        temp = line.strip().split('\t')
        if not UserVals.label_map.__contains__(temp[0]):
            UserVals.label_map[temp[0]] = {}
        UserVals.label_map[temp[0]][temp[1]] = int(temp[2])
    return UserVals


def main():
    init_version()
    verify_package()
    init_argparse()
    plt.ioff()  # disable interactiveness
    create_graphs(_software_flag)

if __name__ == "__main__":
    main()