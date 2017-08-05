//
// Created by harta on 8/4/17.
//

#include <chrono>
#include <string>
#include <ios>
#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <pstream.h>
#include <thread>
#include <boost/program_options/variables_map.hpp>
#include <unordered_map>
#include <sstream>
#include "EntapGlobals.h"


/**
 * ======================================================================
 * Function void print_debug(std::string    msg)
 *
 * Description          - Handles printing to EnTAP debug file
 *                      - Adds timestamp to each entry
 *
 * Notes                - None
 *
 * @param msg           - Message to be sent to debug file
 * @return              - None
 *
 * =====================================================================
 */
void print_debug(std::string msg) {

    std::chrono::time_point<std::chrono::system_clock> current;
    std::time_t time;

    current = std::chrono::system_clock::now();
    time = std::chrono::system_clock::to_time_t(current);
    std::string out_time(std::ctime(&time));
    std::ofstream debug_file(DEBUG_FILE_PATH, std::ios::out | std::ios::app);
    debug_file << out_time.substr(0,out_time.length()-1) << ": " + msg << std::endl;
    debug_file.close();
}


/**
 * ======================================================================
 * Function void print_statistics(std::string    &msg)
 *
 * Description          - Handles printing to EnTAP statistics/log file
 *
 * Notes                - None
 *
 * @param msg           - Message to be sent to log file
 * @return              - None
 *
 * =====================================================================
 */
void print_statistics(std::string &msg) {
    std::ofstream log_file(LOG_FILE_PATH, std::ios::out | std::ios::app);
    log_file << msg << std::endl;
    log_file.close();
}


bool file_exists(std::string path) {
    /* Non-boost implementation
    struct stat buff;
    return (stat(path.c_str(), &buff) == 0);
    */
    return boost::filesystem::exists(path);
}


int execute_cmd(std::string cmd, std::string out_path) {
    std::ofstream out_file(out_path+".out", std::ios::out | std::ios::app);
    std::ofstream err_file(out_path+".err", std::ios::out | std::ios::app);
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                err_file.write(buf, n);
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                out_file.write(buf, n).flush();
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    out_file.close();
    err_file.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}
// todo, may want to handle differently
// TODO change to sending map of flags as command
int execute_cmd(std::string cmd) {
    const redi::pstreams::pmode mode = redi::pstreams::pstdout|redi::pstreams::pstderr;
    redi::ipstream child(cmd, mode);
    char buf[1024];
    std::streamsize n;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1]) {
        if (!finished[0]) {
            while ((n = child.err().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[0] = true;
                if (!finished[1])
                    child.clear();
            }
        }
        if (!finished[1]) {
            while ((n = child.out().readsome(buf, sizeof(buf))) > 0)
                continue;
            if (child.eof()) {
                finished[1] = true;
                if (!finished[0])
                    child.clear();
            }
        }
    }
    child.close();
    if (child.rdbuf()->exited())
        return child.rdbuf()->status();
    return 1;
}


int get_supported_threads(boost::program_options::variables_map &user_map) {

    unsigned int supported_threads;
    int          threads;

    supported_threads = std::thread::hardware_concurrency();
    if (user_map["threads"].as<int>() > supported_threads) {
        print_debug("Specified thread number is larger than available threads,"
                                       "setting threads to " + std::to_string(supported_threads));
        threads = supported_threads;
    } else {
        threads = user_map["threads"].as<int>();
    }
    return threads;
}


std::string generate_command(std::unordered_map<std::string,std::string> &map,std::string exe_path) {
    std::stringstream ss;
    ss << exe_path << " ";
    for (auto &pair : map)ss << pair.first << " " << pair.second << " ";
    std::string out = ss.str();
    return out;
}
