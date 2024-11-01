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
import json

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
GRAPH_TYPE_LOG_VISUAL = 8

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
    elif graph_type == GRAPH_TYPE_LOG_VISUAL:
        input_vals = parse_inputs(gTextFilePath)
        create_gen_statistics(gOutputPath, input_vals)
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



def parse_inputs(file_path):
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split(':')
            if len(parts) >= 2:
                key, value = parts[0], ':'.join(parts[1:]).strip()
                data_dict[key.strip()] = value.strip()

    return data_dict


def create_gen_statistics(files, label_vals):
    with open(files, 'w') as outfile:
        input_no_filter = label_vals["No Filter"]
        no_filter = [int(item.strip().strip("'")) for item in input_no_filter.replace('[', '').replace(']', '').split(',')]

        input_frame_selection = label_vals["Frame Selection"]
        frame_select = [int(item.strip().strip("'")) for item in input_frame_selection.replace('[', '').replace(']', '').split(',')]

        input_expression_filter = label_vals["Expression Filter"]
        expression_filter = [int(item.strip().strip("'")) for item in input_expression_filter.replace('[', '').replace(']', '').split(',')]

        input_sim_annotated_no_contam = label_vals["Sim Search Annotated with no contam"]
        sim_annotated_no_contam = [int(item.strip().strip("'")) for item in input_sim_annotated_no_contam.replace('[', '').replace(']', '').split(',')]

        input_sim_annotated_contam = label_vals["Sim Search Annotated with contam"]
        sim_annotated_contam = [int(item.strip().strip("'")) for item in input_sim_annotated_contam.replace('[', '').replace(']', '').split(',')]

        input_sim_gene_annotated_contam = label_vals["Sim Search or Gene Families Annotated with contam"]
        sim_gene_annotated_contam = [int(item.strip().strip("'")) for item in input_sim_gene_annotated_contam.replace('[', '').replace(']', '').split(',')]

        input_sim_gene_annotated_no_contam = label_vals["Sim Search or Gene Families Annotated with no contam"]
        sim_gene_annotated_no_contam = [int(item.strip().strip("'")) for item in input_sim_gene_annotated_no_contam.replace('[', '').replace(']', '').split(',')]

        input_gene_annotated_contam = label_vals["Gene Families Annotated With Contam"]
        gene_annotated_contam = [int(item.strip().strip("'")) for item in input_gene_annotated_contam.replace('[', '').replace(']', '').split(',')]

        input_gene_annotated_no_contam = label_vals["Gene Families Annotated with no contam"]
        gene_annotated_no_contam = [int(item.strip().strip("'")) for item in input_gene_annotated_no_contam.replace('[', '').replace(']', '').split(',')]

        input_sim_unannotated = label_vals["Unannotated By Sim Search"]
        sim_unannotated = [int(item.strip().strip("'")) for item in input_sim_unannotated.replace('[', '').replace(']', '').split(',')]

        input_gene_unannotated = label_vals["Unannotated By Gene Families"]
        gene_unannotated = [int(item.strip().strip("'")) for item in input_gene_unannotated.replace('[', '').replace(']', '').split(',')]

        avg_seq = label_vals["Average Sequence Length(BP)"]
        n_num = label_vals["n50"]
        long_seq = label_vals["Longest Sequence(BP)"]
        short_seq = label_vals["Shortest Sequence(BP)"]
        align_seq = label_vals["align_seq"]
        align_percent = "{:.1f}".format(float(label_vals["align_percent"]))
        contam_seq = label_vals["contam_seq"]
        contam_percent = "{:.1f}".format(float(label_vals["contam_percent"]))
        unique_seq = label_vals["unique_seq"]
        unique_percent = "{:.1f}".format(float(label_vals["unique_percent"]))
        go = label_vals["go"]
        go_percent = "{:.1f}".format(float(label_vals["go_percent"]))
        unique_seq_ann_sim = label_vals["unique_seq_ann_sim"]
        unique_seq_ann_sim_percent = "{:.1f}".format(float(label_vals["unique_seq_ann_sim_percent"]))
        unique_seq_ann_gene = label_vals["unique_seq_ann_gene"]
        unique_seq_ann_gene_percent = "{:.1f}".format(float(label_vals["unique_seq_ann_gene_percent"]))
        unique_seq_total = label_vals["unique_seq_total"]
        unique_seq_total_percent = "{:.1f}".format(float(label_vals["unique_seq_total_percent"]))
        total_sequences = label_vals["Total sequences"]

        exp_sequences = label_vals["Expression Retained"]
        exp_average_sequences = label_vals["Expression Average Sequence"]
        exp_longest_sequence = label_vals["Expression Longest Sequence"]
        exp_shortest_sequence = label_vals["Expression Shortest Sequence"]
        exp_n50 = label_vals["Expression N50"]

        frame_sequences = label_vals["Frame Retained"]
        frame_average_sequences = label_vals["Frame Average Sequence"]
        frame_longest_sequence = label_vals["Frame Longest Sequence"]
        frame_shortest_sequence = label_vals["Frame Shortest Sequence"]
        frame_n50 = label_vals["Frame N50"]

        config_file = label_vals["config file"]
        fasta_file = label_vals["transcriptome"]
        databases = label_vals["databases"].split(",")
        databases_js_array = "[" + ", ".join("'{}'".format(db.strip()) for db in databases) + "]"
        out_dir = label_vals["out dir"]
        e_value = label_vals["e-value"]
        fpkm_input = "{:.1f}".format(float(label_vals["fpkm"]))
        t_coverage = label_vals["t-coverage"]
        q_coverage = label_vals["q-coverage"]
        if(label_vals["runN"] != "0"):
            runN = "runN"
            runP = "false"
        else:
            runN = "runP"
            runP = "true"
        diamond_database = label_vals["diamond path"]
        interpro = label_vals["interpro"]
        if(label_vals["no-trim"] != "0"):
            no_trim = "true"
        else:
            no_trim = "false"
        hgt_gff = label_vals["hgt-gff"]
        hgt_donor = label_vals["hgt-donor"]
        hgt_recipient = label_vals["hgt-recipient"]
	outfile.write('''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnTAP General Statistics Page</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
    <style>

.workflow-wrapper {{
            position: relative;
            width: 100%;
            max-width: 2048px;
            margin: 0 auto;
            padding: 20px;
            background-color: #ffffff;
            border-radius: 8px;
        }}
        .workflow-container {{
            position: relative;
            width: 100%;
            height: 900px;  /* Adjusted height for better vertical spacing */
        }}
        .workflow-box {{
            position: absolute;
            background-color: #fff;
            border: 2px solid #000;
            border-radius: 8px;
            padding: 10px;
            width: 16%;  /* Adjusted width for better spacing */
            text-align: center;
            font-size: 14px;
            box-sizing: border-box;
        }}
        .workflow-box.orange {{
            background-color: #ffcc80;
        }}
        .workflow-box.blue {{
            background-color: #81d4fa;
        }}
        .workflow-box.gray {{
            background-color: #c9bcbc;
        }}
        .arrow {{
            position: absolute;
            width: 0;
            height: 0;
            border-style: solid;
        }}
        .arrow-right {{
            border-width: 10px 0 10px 20px;
            border-color: transparent transparent transparent black;
        }}
        .arrow-down {{
            border-width: 20px 10px 0 10px;
            border-color: black transparent transparent transparent;
        }}

        /* Positioning elements for grouped appearance and spacing */
        #box1 {{ top: 5%; left: 5%; }}
        #box2 {{ top: 15%; left: 5%; }}  
        #box3 {{ top: 10%; left: 26%; }}
        #box5 {{ top: 10%; left: 47%; }}
        #box6 {{ top: 22%; left: 47.25%; }}
        #box7 {{ top: 35%; left: 38.25%; }}
        #box8 {{ top: 36%; left: 56.25%; }}
        #box9 {{ top: 47%; left: 75.75%; }}
        #box10 {{ top: 35%; left: 75.5%; }}
        #box11 {{ top: 40%; left: 5%; }}

        #arrow1 {{ top: 13%; left: 23%; width: 5%; }}
        #arrow2 {{ top: 13%; left: 28%; width: 5%; }}
        #arrow3 {{ top: 13%; left: 44%; width: 5%; }}
        #arrow4 {{ top: 36%; left: 65%; width: 5%; }}
        #arrow5 {{ top: 19%; left: 54%; height: 15%; }}
        #arrow6 {{ top: 58.5%; left: 50%; height: 15%; }}
        #arrow7 {{ top: 37.5%; left: 73%; width: 10%; }}
        #arrow9 {{ top: 44%; left: 82.5%; height: 10%; }}
        #arrow10 {{ top: 32.5%; left: 54%; height: 5%; }}

        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #f8f9fa;
            color: #333;
            margin: 0;
            padding: 0;
        }}
        body.custom-body {{
            font-family: 'Roboto', sans-serif;
            background-color: #f4f4f9;
            color: #333;
            margin: 0;
            padding: 20px;
        }}
        .custom-container {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
        }}
        .custom-box {{
            background-color: #fff;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            padding: 15px;
            width: calc(33% - 20px);
            box-sizing: border-box;
            overflow: hidden; /* This will keep the text within the box */
        }}
        .custom-box h2 {{
            margin-top: 0;
            margin-bottom: 10px;
            font-size: 1.3em;
            color: #4a90e2;
        }}
        .custom-box p, .custom-box ul {{
            margin-bottom: 5px;
            font-size: 0.9em;
            line-height: 1.4;
        }}
        .custom-box ul {{
            list-style-type: none;
            padding: 0;
        }}
        .custom-box ul li {{
            margin: 3px 0;
            padding: 3px;
            background-color: #f9f9f9;
            border-radius: 4px;
        }}
        @media (max-width: 768px) {{
            .custom-box {{
                width: calc(50% - 20px);
            }}
        }}
        @media (max-width: 480px) {{
            .custom-box {{
                width: 100%;
            }}
        }}
        header {{
            background-color: #007bff;
            color: #fff;
            padding: 20px 0;
        }}
        
        .config-container {{
            font-family: 'Arial', sans-serif;
            margin: 20px;
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        }}
        .config-wrapper {{
            max-width: 800px;
            margin: 0 auto;
            background-color: #fff;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        }}
        .config-header h2 {{
            border-bottom: 3px solid #007BFF;
            padding-bottom: 10px;
            color: #007BFF;
            font-size: 24px;
            margin-bottom: 20px;
        }}
        .subsection-header h3 {{
            border-bottom: 2px solid #007BFF;
            padding-bottom: 5px;
            color: #007BFF;
            font-size: 20px;
            margin-top: 30px;
        }}
        .config-section p {{
            font-size: 16px;
            line-height: 1.8;
            color: #333;
        }}
        .section-columns {{
    display: flex;
    justify-content: space-between;
    margin-top: 20px;
    }}

        .config-section p strong {{
            display: block;
            margin-bottom: 10px;
            font-size: 15px;
            color: #555;
        }}
        .config-container {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}
        .config-wrapper {{
            max-width: 800px;
            margin: 0 auto;
            padding: 20px;
            border: 1px solid #ddd;
            border-radius: 8px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
        }}
        .config-header h2 {{
            border-bottom: 2px solid #007BFF;
            padding-bottom: 10px;
            color: #007BFF;
        }}
        .subsection-header h3 {{
            border-bottom: 1px solid #007BFF;
            padding-bottom: 5px;
            color: #007BFF;
            margin-top: 20px;
        }}
        .config-section p {{
            font-size: 14px;
            line-height: 1.6;
        }}
        .config-section p strong {{
            display: block;
            margin-bottom: 8px;
            font-size: 13px;
        }}
        .entap-banner, .header-content {{
            display: flex;
            justify-content: center;
            align-items: center;
        }}
        .entap-logo {{
            max-width: 250px;
            margin-right: 20px;
        }}
        .entap-info {{
            font-size: 24px;
        }}
        .entap-title {{
            margin-bottom: 5px;
            font-weight: bold;
        }}
        .lab-info {{
            font-size: 18px;
        }}
        .section-container {{
            background-color: #fff;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            margin-bottom: 20px;
            padding: 20px;
        }}
        .section-header {{
            background-color: #343a40;
            border-radius: 8px 8px 0 0;
            padding: 10px;
            color: #fff;
            text-align: center;
            margin-bottom: 10px;
        }}
        .new-config-container {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}
        .new-config-wrapper {{
            border: 1px solid #ccc;
            margin-bottom: 20px;
            border-radius: 5px;
            overflow: hidden;
        }}
        .new-config-header {{
            background-color: #f4f4f4;
            padding: 10px;
            border-bottom: 1px solid #ccc;
        }}
        .new-config-header h2 {{
            margin: 0;
        }}
        .new-config-section {{
            width: 100%;
            border-collapse: collapse;
        }}
        .new-config-section th, .new-config-section td {{
            border: 1px solid #ccc;
            padding: 8px;
            text-align: left;
        }}
        .new-config-section th {{
            background-color: #f4f4f4;
        }}
        .circle-container {{
    text-align: center;
    margin-bottom: 20px;
}}

.circle {{
    position: relative;
    width: 100px; /* Circle size */
    height: 100px; /* Circle size */
    border-radius: 50%; /* Make it round */
    border: 2px solid #007bff; /* Circle border color */
    display: flex; /* Enables flexbox layout */
    justify-content: center; /* Center children horizontally */
    align-items: center; /* Center children vertically */
    overflow: hidden;
    background-color: #f8f9fa; /* Background color */
    margin: 0 auto; /* Center the circle */
}}

.circle::after {{
            content: '';
            position: absolute;
            width: 100%; /* Fill width */
            background-color: #90ee90; /* Fill color */
            bottom: 0;
            left: 0;
            height: 0; /* Initial fill height */
            transition: height 0.3s ease; /* Smooth transition */
        }}
        .circle[data-percent="{unique_seq_ann_sim_percent}"]::after {{
            height: {unique_seq_ann_sim_percent}%;
        }}
        .circle[data-percent="{unique_seq_ann_gene_percent}"]::after {{
            height: {unique_seq_ann_gene_percent}%;
        }}
        .circle[data-percent="{unique_seq_total_percent}"]::after {{
            height: {unique_seq_total_percent}%;
        }}
        .circle[data-percent="{align_percent}"]::after {{
            height: {align_percent}%;
        }}
        .circle[data-percent="{contam_percent}"]::after {{
            height: {contam_percent}%;
        }}
        .circle[data-percent="{go_percent}"]::after {{
            height: {go_percent}%;
        }}
        .circle[data-percent="{unique_percent}"]::after {{
            height: {unique_percent}%;
        }}
        .circle-content {{
            position: relative;
            z-index: 1;
            color: #333;
            font-size: 14px;
            text-align: center;
        }}
        .export-button {{
            margin-top: 20px;
            width: 100%;
        }}
        .chart-container {{
            max-width: 1200px;
            margin: 0 auto;
            border-radius: 8px;
            background-color: #fff;
            padding: 20px;
        }}
        .flex-container {{
            display: flex;
            justify-content: space-around;
            width: 100%;
        }}
        @media (max-width: 768px) {{
            .entap-banner, .header-content {{
                flex-direction: column;
                text-align: center;
            }}
            .entap-logo {{
                margin-bottom: 15px;
            }}
            .section-container, .chart-container {{
                padding: 15px;
            }}
        }}
    </style>
</head>
<body>

<header>
    <div class="entap-banner">
        <div class="header-content">
            <div class="entap-info">
                <h1 class="entap-title">EnTAP 1.3</h1>
                <p class="lab-info">Functional Annotation Report</p>
            </div>
        </div>
    </div>
</header>


<div class="new-config-container">
    <!-- Directory Of Run Section -->
    <div class="new-config-wrapper">
        <div class="new-config-header">
            <h2>Directory of Run</h2>
        </div>
        <table id="databaseTable" class="new-config-section">
            <tr>
                <th>Parameter</th>
                <th>Description</th>
            </tr>
            <tr>
                <td>Path To EnTAP Configuration File</td>
                <td>{config_file}</td>
            </tr>
            <tr>
                <td>FASTA (input)</td>
                <td>{fasta_file}</td>
            </tr>
            <tr>
                <td>Output Directory (out-dir)</td>
                <td>{out_dir}</td>
            </tr>
        </table>
    </div>

<script>
document.addEventListener("DOMContentLoaded", function() {{
    const databases = {databases_js_array};

    const table = document.getElementById("databaseTable");

    databases.forEach((database, index) => {{
        const row = table.insertRow(table.rows.length - 1);
        
        const paramCell = row.insertCell(0);
        const descCell = row.insertCell(1);

        paramCell.textContent = `Database ${{index + 1}}`;
        descCell.textContent = database;
    }});
}});
</script>


    <!-- Run Configuration Section -->
    <div class="new-config-wrapper">
        <div id = "hgt" class="new-config-header">
            <h2>Run Configuration</h2>
        </div>
        <table class="new-config-section">
            <tr>
                <th>Parameter</th>
                <th>Description</th>
            </tr>
            <tr>
                <td>Diamond Database (database)</td>
                <td>{diamond_database}</td>
            </tr>
            <tr>
                <td>FASTA Header Format (no-trim)</td>
                <td>{no_trim}</td>
            </tr>
            <tr>
                <td>Target Sequence Coverage (tcoverage)</td>
                <td>{t_coverage}</td>
            </tr>
            <tr>
                <td>Query Sequence Coverage (qcoverage)</td>
                <td>{q_coverage}</td>
            </tr>
            <tr>
                <td>E-Value (e-value)</td>
                <td>{e_value}</td>
            </tr>
            <tr>
                <td>Sequence Mode (RunP or RunN)</td>
                <td>{runN}</td>
            </tr>
            <tr>
                <td>FPKM</td>
                <td>{fpkm_input}</td>
            </tr>
            <tr>
                <td>Protein (interproscan-db)</td>
                <td>{interpro}</td>
            </tr>
            <tr>
                <td>HGT Donor Databases (hgt-donor)</td>
                <td>{hgt_donor}</td>
            </tr>

            <tr>
                <td>HGT Recipient Databases (hgt-recipient)</td>
                <td>{hgt_recipient}</td>
            </tr>
            <tr>
                <td>HGT GFF Path (hgt-gff)</td>
                <td>{hgt_gff}</td>
            </tr>
        </table>
    </div>
</div>


<div class="container">
<div class="section-container">

<div class="workflow-wrapper">
    <div class="workflow-container">
            <img src="https://github.com/harta55/EnTAP/raw/master/docs/source/ENTAP_white_50.jpg?raw=true" alt="EnTAP Logo" class="entap-logo" style="float: right; margin-left: 20px;">
        </div>
        <div id="box1" class="workflow-box gray">Nucleotide Sequences (FASTA)</div>
        <div id="box2" class="workflow-box gray">EnTAP Configuration Files (ini)</div>
        <div id="box3" class="workflow-box orange">Filtered Transcripts (FASTA)</div>
        <div id="arrow1" class="arrow arrow-right"></div>
        <div id="box5" class="workflow-box blue">Database Alignments (tsv)</div>
        <div id="arrow3" class="arrow arrow-right"></div>
        <div id="box6" class="workflow-box gray">Optional Alignment Selection / Contaminant Filtering</div>
        <div id="arrow4" class="arrow arrow-right"></div>
        <div id="box7" class="workflow-box blue">Optimal Database Alignments (tsv)</div>
        <div id="arrow5" class="arrow arrow-down"></div>
        <div id="box8" class="workflow-box orange">Transcriptome (FASTA)</div>
        <div id="arrow6" the-arrow arrow-down"></div>
        <div id="box9" class="workflow-box blue">EnTAP Summary</div>
        <div id="box10" class="workflow-box gray">EggNOG Alignments (tsv)</div>
        <div id="arrow7" class="arrow arrow-right"></div>
        <div id="arrow9" class="arrow arrow-down"></div>
        <div id="arrow10" class="arrow arrow-down"></div>
    </div>
</div>
</div>
</div>

<script>
    var runP = {runP}; 

    if (runP) {{
        var proteinBox = document.createElement("div");
        proteinBox.id = "box11";
        proteinBox.className = "workflow-box orange";
        proteinBox.innerText = "Protein Sequences (FASTA)";
        document.getElementById("box3").insertAdjacentElement("afterend", proteinBox);

        var arrowToProtein = document.createElement("div");
        arrowToProtein.className = "arrow arrow-right";
        proteinBox.insertAdjacentElement("beforebegin", arrowToProtein);

        proteinBox.style.position = "absolute";
        proteinBox.style.top = "10.125%"; 
        proteinBox.style.left = "46%";
        arrowToProtein.style.position = "absolute";
        arrowToProtein.style.top = "12.125%"; 
        arrowToProtein.style.left = "43.5%"; 

        var downArrowAfterProtein = document.createElement("div");
        downArrowAfterProtein.className = "arrow arrow-down";
        proteinBox.insertAdjacentElement("afterend", downArrowAfterProtein);

        downArrowAfterProtein.style.position = "absolute";
        downArrowAfterProtein.style.top = "19%"; 
        downArrowAfterProtein.style.left = "53%";

        document.getElementById("box5").style.top = "22%"; 
        document.getElementById("box5").style.left = "46%";

        document.getElementById("box6").style.top = "35%"; 
        document.getElementById("box6").style.left = "46%";

        document.getElementById("box7").style.top = "36%"; 
        document.getElementById("box7").style.left = "66%";

        document.getElementById("box8").style.top = "50.5%"; 
        document.getElementById("box8").style.left = "46.25%";

        document.getElementById("box9").style.top = "63%"; 
        document.getElementById("box9").style.left = "67.5%";

        document.getElementById("box10").style.top = "62%"; 
        document.getElementById("box10").style.left = "46.5%"; 

        document.getElementById("arrow3").style.top = "38.5%";
        document.getElementById("arrow3").style.left = "63%";

        document.getElementById("arrow4").style.display = "none";

        document.getElementById("arrow5").style.top = "46.5%";
        document.getElementById("arrow5").style.left = "53.5%";

        document.getElementById("arrow7").style.top = "64.5%";
        document.getElementById("arrow7").style.left = "64%";

        document.getElementById("arrow9").style.top = "58%";
        document.getElementById("arrow9").style.left = "53.5%";

        document.getElementById("arrow10").style.top = "31%";
        document.getElementById("arrow10").style.left = "53%"; 

    }}
</script>



<div class="container">
    <!-- Totals Section -->
    <div class="section-container">
        <div class="section-header">
            <h2>Final Results</h2>
        </div>
        <div class="section">
            <p><strong>Total retained sequences (after filtering and/or frame selection): </strong>38,640</p>
            <div class="row">
                <div class="col-md-4">
                    <div class="circle-container">
                        <div class="circle" data-percent="{unique_seq_ann_sim_percent}">
                            <div class="circle-content">{unique_seq_ann_sim_percent}%</div>
                        </div>
                        <p><strong>Total unique sequences annotated (similarity search alignments only): </strong>{unique_seq_ann_sim}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;"> 
                                <small style="margin-left: 5px;">entap_results.tsv</small>
                            </span>
                        </p>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="circle-container">
                        <div class="circle" data-percent="{unique_seq_ann_gene_percent}">
                            <div class="circle-content">{unique_seq_ann_gene_percent}%</div>
                        </div>
                        <p><strong>Total unique sequences annotated (gene family assignment only): </strong>{unique_seq_ann_gene}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;"> 
                                <small style="margin-left: 5px;">entap_results.tsv</small>
                            </span>
                        </p>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="circle-container">
                        <div class="circle" data-percent="{unique_seq_total_percent}">
                            <div class="circle-content">{unique_seq_total_percent}%</div>
                        </div>
                        <p><strong>Total unique sequences annotated (gene family and/or similarity search): </strong>{unique_seq_total}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;"> 
                                <small style="margin-left: 5px;">entap_results.tsv</small>
                            </span>
                        </p>
                    </div>
                </div>
            </div>
        </div>
    </div>
</div>

<div class="section-container">
    <div class="section-header">
        <h2 style="font-size: 28px;">Transcriptomes</h2>
    </div>
    <div class="section-columns" style="display: flex; justify-content: space-between; align-items: center;">
        <!-- First Section -->
        <div class="section" style="flex: 1; padding: 10px;">
            <p style="font-size: 18px;"><strong>Total Input Sequences:</strong> {total_sequences}</p>
            <p style="font-size: 18px;"><strong>Average sequence length (bp):</strong> {avg_seq}</p>
            <p style="font-size: 18px;"><strong>Longest sequence (bp):</strong> {long_seq}</p>
            <p style="font-size: 18px;"><strong>Shortest sequence (bp):</strong> {short_seq}</p>
            <p style="font-size: 18px;"><strong>N50 (bp):</strong> {n_num}</p>
        </div>

        <!-- Text above arrow and clean SVG arrow -->
        <div class="arrow-section" style="flex: 0; padding: 0 20px; text-align: center;">
            <!-- Label above the arrow -->
            <p style="font-size: 24px; font-weight: bold; margin-bottom: 10px;">Expression Analysis</p>
            <!-- Clean Arrow SVG -->
            <svg width="300" height="150" xmlns="http://www.w3.org/2000/svg">
                <polygon points="10,75 220,75 220,35 280,75 220,115 220,75 10,75" style="fill: white; stroke: black; stroke-width: 4;" />
            </svg>
        </div>

        <!-- Second Section -->
        <div class="section" style="flex: 1; padding: 10px;">
            <p style="font-size: 18px;"><strong>Total Input Sequences:</strong> {exp_sequences}</p>
            <p style="font-size: 18px;"><strong>Average sequence length (bp):</strong> {exp_average_sequences}</p>
            <p style="font-size: 18px;"><strong>Longest sequence (bp):</strong> {exp_longest_sequence}</p>
            <p style="font-size: 18px;"><strong>Shortest sequence (bp):</strong> {exp_shortest_sequence}</p>
            <p style="font-size: 18px;"><strong>N50 (bp):</strong> {exp_n50}</p>
        </div>

        <!-- Text above arrow and clean SVG arrow -->
        <div class="arrow-section" style="flex: 0; padding: 0 20px; text-align: center;">
            <!-- Label above the arrow -->
            <p style="font-size: 24px; font-weight: bold; margin-bottom: 10px;">Frame Selection</p>
            <!-- Clean Arrow SVG -->
            <svg width="300" height="150" xmlns="http://www.w3.org/2000/svg">
                <polygon points="10,75 220,75 220,35 280,75 220,115 220,75 10,75" style="fill: white; stroke: black; stroke-width: 4;" />
            </svg>
        </div>

        <!-- Third Section -->
        <div class="section" style="flex: 1; padding: 10px;">
            <p style="font-size: 18px;"><strong>Total Input Sequences:</strong> {frame_sequences}</p>
            <p style="font-size: 18px;"><strong>Average sequence length (bp):</strong> {frame_average_sequences}</p>
            <p style="font-size: 18px;"><strong>Longest sequence (bp):</strong> {frame_longest_sequence}</p>
            <p style="font-size: 18px;"><strong>Shortest sequence (bp):</strong> {frame_shortest_sequence}</p>
            <p style="font-size: 18px;"><strong>N50 (bp):</strong> {frame_n50}</p>
        </div>
    </div>
</div>

    <!-- Similarity Search Section -->
    <div class="row">
        <div class="col-md-6">
            <div class="section-container">
                <div class="section-header" style="background-color: #007bff;">
                    <h2>Similarity Search</h2>
                </div>
                <div class="section">
                    <div class="circle-container">
                        <div class="circle" data-percent="{align_percent}">
                            <div class="circle-content">{align_percent}%</div>
                        </div>
                        <p>
                            <strong>Total unique sequences with an alignment: </strong>{align_seq}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;">
                                <small style="margin-left: 5px;">diamond_annotated.tsv</small>
                            </span>
                        </p>
                    </div>
                    <div class="circle-container">
                        <div class="circle" data-percent="{contam_percent}">
                            <div class="circle-content">{contam_percent}%</div>
                        </div>
                        <p>
                            <strong>Total alignments flagged as a contaminant: </strong>{contam_seq}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;">
                                <small style="margin-left: 5px;">diamond_annotated.tsv</small>
                            </span>
                        </p>
                    </div>
                </div>
            </div>
        </div>
    

        <!-- Gene Families Section -->
        <div class="col-md-6">
            <div class="section-container">
                <div class="section-header" style="background-color: #e67e22;">
                    <h2>Gene Families<h2></h2>
                </div>
                <div class="section">
                    <div class="circle-container">
                        <div class="circle" data-percent="{unique_percent}">
                            <div class="circle-content">{unique_percent}%</div>
                        </div>
                        <p>
                            <strong>Total unique sequences with family assignment: </strong>{unique_seq}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;">
                                <small style="margin-left: 5px;">eggnog_annotated.tsv</small>
                            </span>
                        </p>
                    </div>
                    <div class="circle-container">
                        <div class="circle" data-percent="{go_percent}">
                            <div class="circle-content">{go_percent}%</div>
                        </div>
                        <p>
                            <strong>Total unique sequences with at least one GO term:  </strong>{go}
                            <span style="display: inline-flex; align-items: center; margin-left: 10px;">
                                <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJoAAACUCAMAAABcK8BVAAAAYFBMVEX///8AAABnZ2cqKiqoqKjZ2dkjIyOkpKS/v789PT0QEBAWFhbn5+fIyMji4uL09PRKSkp2dnZeXl5sbGyvr69YWFiZmZkICAjQ0NC3t7eEhIR8fHw0NDTt7e2Li4tSUlL531sDAAADMklEQVR4nO2c2ZqiMBBGWYQWRLYADQj6/m85Kgi2hsXKn29mmjq3kXCoJFgkBMNgGOa/JXKqnSKVE+kwy2wTgL3Dy50QYjdStNkZZWaaVgs1c3ycGjhu96CJMlQiKB5xQ/a3/Fqhp1xhNrQpwqnHhtSXjW2KixtOzc3vbgIg1YFT88Ma299wat+O0bkJkBtUrb99g8YCVg0aN7AaMm5oNWDc4Gq4uOHVYHHToGYkmHuvDjVQ3LSo9f1NMW561B5uSnHTpBZ1bZqoVKpJDTEWtKn1Y+FEr1SfmnLcgFnuq5pq3DBqe1/a509KY+GnWnC2xBxWHUpraZubQ/zy60R0T5JngFpkmUvk8sfg2Qftw1tTf67mfC+qvfeojnzuGFrYfqi18aJZE8jraecCflFXM8IknUdMt01Zp7kE93oCG6B27W7zzNYlPcKCqcFhNQqsRoHVKLAaBVajwGoUUGqhaLxXmnzflwan91IpjRjzdJBae5Gngt2JVqTnY/I45OkgtakEPOvE15s9PZSC1KKJqHVJd5R+ELUh6UT1tUB48Rvp2NckpTI8MT5BbGCEaoDVKLAaBVajwGoUtqAWHiE8Taui1Hau7x+U8X13h1YLIa+L3bCHuKFSyQNK7YBOJY1i6ZRrKYYqYcOg3EMoxxq3cPPAw2oUWI0Cq1FgNQobUGsLkaxCFI8JtLKWFBfjMjhqEitZn1vU/cXI86gEPYnlrDczzS4w2UQpOl9bsfY+0E/uTVxNjJ4wNZx03dxeHFuPDDu7SErTcYkeN0IX1t5li/DzpRu4eWiA1SiwGgVWo8BqFGBqUQDh6Z8KpVb5ICq0WtCsT4rmGd9b3MDUX41Sq4cqUWrRUXVTcsdxHAdbuHngYTUKrEaB1SiwGgVWkxF+Fe9k/8J/aCRf3T3B31/7nGoiK4Lna59Tys08dJZLoYglr5hejkP53xyhreQDH0/7wPjmQYHVKLAaBVajwGoUfqWaZ2r4Otkzt00KHvlIW75jHML9RULatd9Xgd3zlybOt63I/QYU0lVphtoqe/1q+2ULORVsAlfOoVp2mCJMGtvVhN0kamMsCh1NhFo+N8kwv4s/NxVdOGdiz2kAAAAASUVORK5CYII=" 
                                     alt="File Icon" 
                                     style="height: 16px; width: auto;">
                                <small style="margin-left: 5px;">eggnog_annotated.tsv</small>
                            </span>
                        </p>
                    </div>
                </div>
            </div>
        </div>
        <div class="section-container" style="display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 40px; margin: 0 auto;">

            <!-- Transcriptomes Section with 3 Violin Plots -->
            <div class="chart-container" style="width: 1200px; height: 600px; margin: 0 auto;">
                <div id="transcriptomes" style="width:100%; height:100%;"></div>
            </div>
        
            <!-- Similarity Search Section with 4 Violin Plots -->
            <div class="chart-container" style="width: 1200px; height: 600px; margin: 0 auto;">
                <div id="similarity-search" style="width:100%; height:100%;"></div>
            </div>
        
            <!-- Gene Family Section with 3 Violin Plots -->
            <div class="chart-container" style="width: 1200px; height: 600px; margin: 0 auto;">
                <div id="gene-family" style="width:100%; height:100%;"></div>
            </div>
        
            <!-- Combined Section: Similarity Search and/or Gene Family -->
            <div class="chart-container" style="width: 1200px; height: 600px; margin: 0 auto;">
                <div id="combined" style="width:100%; height:100%;"></div>
            </div>
        
        </div>
<div style = "height: 50px;"></div>




    <!-- Export PDF Button -->
    <div class="container">
        <button type="button" class="btn btn-primary btn-lg export-button" onclick="exportPDF()">Export as PDF</button>
    </div>
</div>

<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>

<script>
    // Function to create violin plots
    function createViolinPlot(id, traces, layout) {{
        Plotly.newPlot(id, traces, layout);
    }}

    // Transcriptomes Violin Plot (3 groups)
    var trace1 = {{ type: "violin", y: {no_filter}, name: "No Filtering", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace2 = {{ type: "violin", y: {frame_select}, name: "Frame Selected", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace3 = {{ type: "violin", y: {expression_filter}, name: "Expression (FPKM > 0.5)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    
    var transcriptomesLayout = {{ title: "Transcriptome", yaxis: {{ title: "Sequence Length" }}, xaxis: {{ showticklabels: false }},showlegend: true }};
    createViolinPlot("transcriptomes", [trace1, trace2, trace3], transcriptomesLayout);

    // Similarity Search Violin Plot (3 groups)
    var trace4 = {{ type: "violin", y: {sim_annotated_no_contam}, name: "Annotated (no contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace5 = {{ type: "violin", y: {sim_annotated_contam}, name: "Annotated (contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace7 = {{ type: "violin", y: {sim_unannotated}, name: "Unannotated", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    
    var similaritySearchLayout = {{ title: "Similarity Search", yaxis: {{ title: "Sequence Length" }}, xaxis: {{ showticklabels: false }}, showlegend: true }};
    createViolinPlot("similarity-search", [trace4, trace5, trace7], similaritySearchLayout);

    // Gene Family Violin Plot (3 groups)
    var trace8 = {{ type: "violin", y: {gene_annotated_no_contam}, name: "Annotated (no contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace9 = {{ type: "violin", y: {gene_annotated_contam}, name: "Annotated (contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace10 = {{ type: "violin", y: {gene_unannotated}, name: "Unannotated", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    
    var geneFamilyLayout = {{ title: "Gene Families", yaxis: {{ title: "Sequence Length" }}, xaxis: {{ showticklabels: false }},showlegend: true }};
    createViolinPlot("gene-family", [trace8, trace9, trace10], geneFamilyLayout);

    // Combined Violin Plot (2 groups)
    var trace11 = {{ type: "violin", y: {sim_gene_annotated_contam}, name: "Annotated (no contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    var trace12 = {{ type: "violin", y: {sim_gene_annotated_no_contam}, name: "Annotated (contam.)", box: {{ visible: true }}, meanline: {{ visible: true }} }};
    
    var combinedLayout = {{ title: "Similarity Search and Gene Families", yaxis: {{ title: "Sequence Length" }}, xaxis: {{ showticklabels: false }},showlegend: true }};
    createViolinPlot("combined", [trace11, trace12], combinedLayout);
</script>




'''.format(diamond_database = diamond_database, frame_sequences = frame_sequences, frame_average_sequences = frame_average_sequences, frame_longest_sequence = frame_longest_sequence, frame_shortest_sequence = frame_shortest_sequence, frame_n50 = frame_n50, hgt_gff = hgt_gff, hgt_donor = hgt_donor, hgt_recipient = hgt_recipient, no_trim = no_trim, interpro = interpro, t_coverage = t_coverage, q_coverage = q_coverage, runP = runP, e_value = e_value, runN = runN, fpkm_input = fpkm_input, out_dir = out_dir, databases_js_array = databases_js_array, fasta_file = fasta_file, config_file = config_file, exp_sequences = exp_sequences, exp_average_sequences = exp_average_sequences, exp_longest_sequence = exp_longest_sequence, exp_shortest_sequence = exp_shortest_sequence, exp_n50 = exp_n50, total_sequences = total_sequences, sim_unannotated = sim_unannotated, gene_unannotated = gene_unannotated, frame_select = frame_select, expression_filter = expression_filter, sim_annotated_contam = sim_annotated_contam, sim_annotated_no_contam = sim_annotated_no_contam, gene_annotated_contam = gene_annotated_contam, gene_annotated_no_contam = gene_annotated_no_contam, sim_gene_annotated_contam = sim_gene_annotated_contam, sim_gene_annotated_no_contam = sim_gene_annotated_no_contam, no_filter = no_filter, avg_seq = avg_seq, short_seq=short_seq,long_seq=long_seq,n_num=n_num,align_seq=align_seq,align_percent=align_percent,contam_seq=contam_seq,contam_percent=contam_percent,unique_seq=unique_seq,unique_percent=unique_percent,go=go,go_percent=go_percent,unique_seq_ann_sim=unique_seq_ann_sim,unique_seq_ann_sim_percent=unique_seq_ann_sim_percent,unique_seq_ann_gene=unique_seq_ann_gene,unique_seq_ann_gene_percent=unique_seq_ann_gene_percent,unique_seq_total=unique_seq_total,unique_seq_total_percent=unique_seq_total_percent))



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
