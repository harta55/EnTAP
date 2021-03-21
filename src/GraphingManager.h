/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2021, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENTAP_GRAPHINGMANAGER_H
#define ENTAP_GRAPHINGMANAGER_H

#include "common.h"
#include "FileSystem.h"

#define ENT_GRAPH_NULL "NULL"

class GraphingManager {

public:

    // WARNING when changing/adding values here, make sure entap_graphing.py
    //  matches ENUM
    typedef enum {
        ENT_GRAPH_TEST=-1,
        // Reserved as a TEST that graphing is working on this system
        ENT_GRAPH_BAR_HORIZONTAL=0,
        ENT_GRAPH_BAR_VERTICAL=1,
        /*
         * Bar Graph Formatting
         *
         */
        ENT_GRAPH_LINE_HORIZONTAL=2,
        ENT_GRAPH_LINE_VERTICAL=3,
        ENT_GRAPH_PIE_CHART=4,
        /*
         * Pie Chart Formatting
         *
         * Each line of the file should contain the Label for that dataset,
         *     followed by a tab character ('\t') and then the number
         *     of entries corresponding to that Label
         *
         * Example:
         *
         * Removed          34570
         * Partial 5 Prime	11554
         * Partial 3 Prime	7666
         * Complete	        37819
         * Internal	        7837
         */
        ENT_GRAPH_BOX_PLOT_VERTICAL=5,
        ENT_GRAPH_BOX_PLOT_HORIZONTAL=6,
        /*
         * Box Plot Formatting
         *
         * The first line of the file should contain the overall X Axis Label,
         *    followed by a tab character ('\t') and then the Y Axis Label
         * Note: Currently, only the Labels in subsequent lines will be shown
         *    on the X-Axis
         * Each subsequent line should contain the Label for that datapoint
         *    (used for the X-Axis labelling), followed by a tab character ('\t')
         *    and then the value associated with the datapoint
         *
         * Example:
         *
         * flag	    Sequence Length
         * Removed	225
         * Removed	1087
         * Removed	1110
         * Removed	441
         * Selected	930
         * Removed	675
         *
         * This will give us a Y-Axis Label of "Sequence Length", and an X-Axis
         * label of "Removed" and "Selected"
         */
        ENT_GRAPH_BAR_STACKED=7
        /*
         * Bar Graph Stacked Formatting
         *
         *
         */

    } ENT_GRAPHING_TYPES;

    struct GraphingData {
        std::string        text_file_path;  // Absolute path to text file containing data to graph
        std::string        graph_title;     // Title of gray
        std::string        x_axis_label;
        std::string        y_axis_label;
        std::string        fig_out_path;    // Absolute path to output file
        ENT_GRAPHING_TYPES graph_type;      // Type of graph we want to create from data
    };

    GraphingManager(std::string, FileSystem *fileSystem);
    ~GraphingManager();
    void graph(GraphingData*);
    bool is_graphing_enabled() const;
    bool initialize_graph_data(GraphingData &graphingData);
    bool add_datapoint(std::string &path, std::list<std::string> list);
    void graph_data(std::string &path);
    void close_graphing_file(std::string &path);

    class EntapGraphBase {
    public:
        EntapGraphBase(GraphingData &graphingData);

        virtual ~EntapGraphBase();

        GraphingData &getMGraphingData();
        bool close_graphing_file();

        virtual bool add_datapoint(std::list<std::string> &list)=0;

    protected:

        GraphingData mGraphingData;
        std::ofstream *mpOutputTextStream;
    };

    class EntapGraphBarVertical : public EntapGraphBase {
    public:
        EntapGraphBarVertical(GraphingData &graphingData);
        virtual ~EntapGraphBarVertical()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr uint16 MINIMUM_ARGS = 2;
    };

    class EntapGraphBarHorizontal : public EntapGraphBase {
    public:
        EntapGraphBarHorizontal(GraphingData &graphingData);
        virtual ~EntapGraphBarHorizontal()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr uint16 MINIMUM_ARGS = 2;
    };

    class EntapGraphLineHorizontal : public EntapGraphBase {
    public:
        EntapGraphLineHorizontal(GraphingData &graphingData);
        virtual ~EntapGraphLineHorizontal()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    };

    class EntapGraphLineVertical : public EntapGraphBase {
    public:
        EntapGraphLineVertical(GraphingData &graphingData);
        virtual ~EntapGraphLineVertical()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    };

    class EntapGraphPieChart : public EntapGraphBase {
    public:
        EntapGraphPieChart(GraphingData &graphingData);
        virtual ~EntapGraphPieChart()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr  uint16 MINIMUM_ARGS = 2;
        static constexpr  uint16 MAXIMUM_ARGS = 2;
    };

    class EntapGraphBoxPlotVertical : public EntapGraphBase {
    public:
        EntapGraphBoxPlotVertical(GraphingData &graphingData);
        virtual ~EntapGraphBoxPlotVertical()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr  uint16 MINIMUM_ARGS = 2;
        static constexpr  uint16 MAXIMUM_ARGS = 2;
    };

    class EntapGraphBoxPlotHorizontal : public EntapGraphBase {
    public:
        EntapGraphBoxPlotHorizontal(GraphingData &graphingData);
        virtual ~EntapGraphBoxPlotHorizontal()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr  uint16 MINIMUM_ARGS = 2;
        static constexpr  uint16 MAXIMUM_ARGS = 2;
    };

    class EntapGraphBarStacked : public EntapGraphBase {
    public:
        EntapGraphBarStacked(GraphingData &graphingData);
        virtual ~EntapGraphBarStacked()= default;
        bool add_datapoint(std::list<std::string> &list) override;
    private:
        static constexpr  uint16 MINIMUM_ARGS = 3;
        static constexpr  uint16 MAXIMUM_ARGS = 3;
    };

private:

    const std::string FLAG_GRAPH_DATA_PATH = "-i";
    const std::string FLAG_GRAPH_TYPE      = "-g";
    const std::string FLAG_TITLE           = "-t";
    const std::string FLAG_OUT_PATH        = "-o";
    const std::string FLAG_INI_FILE_PATH   = "-f";

    std::string mGraphingPath;  // Absolute path to graphing script
    bool mGraphingEnabled;      // TRUE if graphing is enabled on current system
    std::unordered_map<std::string, EntapGraphBase*> mGraphData;
    FileSystem *mpFileSystem;
};


#endif //ENTAP_GRAPHING_H
