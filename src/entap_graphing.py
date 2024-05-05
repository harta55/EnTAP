#!/usr/bin/env python

# Developed by Alexander Hart
# Plant Computational Genomics Lab
# University of Connecticut
#
# For information, contact Alexander Hart at:
#     entap.dev@gmail.com
#
# Copyright 2017-2024, Alexander Hart, Dr. Jill Wegrzyn
#
# This file is part of EnTAP.
#
# EnTAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EnTAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
#

import argparse
import sys
import collections
import os
from operator import add

# ------------------- GLOBALS -------------------
gTextFilePath = ""
gGraphType = ""
gBasePath = ""
gGraphTitle = ""
gOutputPath = ""
gIniFilePath = ""
gVersion = 0.0
gPlot = None
# -----------------------------------------------

# ------------------ CONSTANTS ------------------
# WARNING Input flags must match those defined in GraphingManager.h
INPUT_FLAG_TEXT_PATH = "-i"
INPUT_FLAG_GRAPH_TYPE = "-g"
INPUT_FLAG_GRAPH_TITLE = "-t"
INPUT_FLAG_OUT_PATH = "-o"
INPUT_FLAG_INI_FILE_PATH = "-f"

# WARNING Graphing Types must match those defined in GraphingManager.h
GRAPH_TYPE_TESTING = -1
GRAPH_TYPE_BAR_HORIZONTAL = 0
GRAPH_TYPE_BAR_VERTICAL = 1
GRAPH_TYPE_LINE_HORIZONTAL = 2
GRAPH_TYPE_LINE_VERTICAL = 3
GRAPH_TYPE_PIE_CHART = 4
GRAPH_TYPE_BOX_PLOT_VERTICAL = 5
GRAPH_TYPE_BOX_PLOT_HORIZONTAL = 6
GRAPH_TYPE_BAR_STACKED = 7

ENTAP_EXIT_OK = 0
ENTAP_EXIT_UNSUPPORTED_SOFTWARE = 1
ENTAP_EXIT_UNSUPPORTED_GRAPH_TYPE = 2

BOX_SEQ_LEN_LABEL = "Sequence Length" # temporary

# --------------CONSTANTS END -------------------
# TODO graphing needs a lot of updating

def init_argparse():
    global gBasePath
    global gTextFilePath
    global gGraphType
    global gGraphTitle
    global gOutputPath
    global gIniFilePath
    parser = argparse.ArgumentParser()
    parser.add_argument(INPUT_FLAG_TEXT_PATH, action='store', dest='stats', help='Path to graphing data', type=str)
    parser.add_argument(INPUT_FLAG_GRAPH_TYPE, action='store', dest='graph', help='Graphing type', type=int)
    parser.add_argument(INPUT_FLAG_GRAPH_TITLE, action='store', dest='title', help='Graph title', type=str)
    parser.add_argument(INPUT_FLAG_OUT_PATH, action='store', dest='path', help='Output path', type=str)
    parser.add_argument(INPUT_FLAG_INI_FILE_PATH, action='store', dest='ini_file', help='Absolute path to INI file', type=str)
    args = parser.parse_args()

    gTextFilePath = args.stats
    gOutputPath = args.path
    gGraphTitle = args.title
    gIniFilePath = args.ini_file
    if gGraphTitle is not None:
        gGraphTitle = gGraphTitle.replace("_", " ")
    gGraphType = args.graph


def init_version():
    global gVersion
    gVersion = sys.version_info[0] + 0.1 * sys.version_info[1]


def verify_package():
    is_2 = 2 <= gVersion < 3
    has_matplot = False
    if is_2:
        import imp
        try:
            imp.find_module("matplotlib")
            has_matplot = True
        except ImportError:
            has_matplot = False
    elif gVersion >= 3.4:
        import importlib.util
        check = importlib.util.find_spec("matplotlib")
        has_matplot = check is not None
    elif not is_2 and gVersion < 3.4:
        import importlib
        spam_loader = importlib.find_loader("matplotlib")
        has_matplot = spam_loader is not None
    else:
        print("Python version not compatable.")
        exit(ENTAP_EXIT_UNSUPPORTED_SOFTWARE)
    if has_matplot:
        global gPlot
        import matplotlib
        import matplotlib.pyplot as gPlot
        gPlot.switch_backend('agg')
    else:
        print("Matplotlib module not found. Not able to graph data.")
        exit(ENTAP_EXIT_UNSUPPORTED_SOFTWARE)


def init_graphs():
    pass


def create_graphs(graph_type):
    if graph_type == GRAPH_TYPE_TESTING:                  # Initial test case
        exit(ENTAP_EXIT_OK)
    elif graph_type == GRAPH_TYPE_BAR_HORIZONTAL:
        InputValues = parse_input(gTextFilePath)
        create_bar(gGraphTitle, gOutputPath, InputValues.values, InputValues.labels,
                   InputValues.xlabel, InputValues.ylabel)
    elif graph_type == GRAPH_TYPE_BAR_VERTICAL:
        InputValues = parse_input(gTextFilePath)
        create_bar(gGraphTitle, gOutputPath, InputValues.values, InputValues.labels,
                   InputValues.xlabel, InputValues.ylabel)
    elif graph_type == GRAPH_TYPE_LINE_HORIZONTAL:
        pass
    elif graph_type == GRAPH_TYPE_LINE_VERTICAL:
        pass
    elif graph_type == GRAPH_TYPE_PIE_CHART:
        InputValues = parse_input(gTextFilePath)
        create_pie(gGraphTitle, gOutputPath, InputValues.values, InputValues.labels)
        pass
    elif graph_type == GRAPH_TYPE_BOX_PLOT_VERTICAL:
        input_vals = parse_input_dict(gTextFilePath)
        create_boxplot(gGraphTitle, gOutputPath, input_vals, BOX_SEQ_LEN_LABEL)
        pass
    elif graph_type == GRAPH_TYPE_BOX_PLOT_HORIZONTAL:
        input_vals = parse_input_dict(gTextFilePath)
        create_boxplot(gGraphTitle, gOutputPath, input_vals, BOX_SEQ_LEN_LABEL)
        pass
    elif graph_type == GRAPH_TYPE_BAR_STACKED:
        InputValues = parse_sim_stack(gTextFilePath)
        create_bar_stacked(gGraphTitle, gOutputPath, InputValues.label_map,
                           InputValues.xlabel, InputValues.ylabel)
        pass
    else:
        exit(ENTAP_EXIT_UNSUPPORTED_GRAPH_TYPE)



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
            p = gPlot.bar(indices, stacked_vals[key], bottom=totals)
        else:
            p = gPlot.bar(indices, stacked_vals[key])
        plts.append(p)
        legend_labels.append(key)
        prev_key = key
        totals = map(add, totals, stacked_vals[key])

    gPlot.legend(plts, legend_labels)
    gPlot.xticks(indices, x_labels)
    gPlot.ylabel(ylab)
    gPlot.title(title)
    gPlot.savefig(file, bbox_inches="tight")


def create_pie(title, file, vals, labels):
    gPlot.pie(vals, labels=labels, autopct=autopct_vals(vals))
    gPlot.axis('equal')
    gPlot.title(title)
    gPlot.gcf().subplots_adjust(bottom=0.15)
    gPlot.savefig(file)


# label_vals = dictionary of lists containing different series. Such as all rejected seqs and kept
def create_boxplot(title, file, label_vals, y_label):
    data = []
    labels = label_vals.keys()
    for key in label_vals.keys():
        temp = []
        for val in label_vals[key]:
            temp.append(val)
        data.append(temp)
    gPlot.ylabel(y_label)
    gPlot.boxplot(data, labels=labels)
    gPlot.title(title)
    gPlot.gcf().subplots_adjust(bottom=0.15)
    gPlot.savefig(file)


def create_bar(title, file, vals, labels, xlabel, ylabel):
    gPlot.barh(range(len(labels)), vals, align='center', alpha=0.5)   # Horizontal bar graph
    gPlot.yticks(range(len(labels)), labels)
    gPlot.ylabel(xlabel)
    gPlot.xlabel(ylabel)
    gPlot.title(title)
    gPlot.savefig(file, bbox_inches="tight")


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
    gPlot.ioff()  # disable interactiveness
    create_graphs(gGraphType)

if __name__ == "__main__":
    main()
