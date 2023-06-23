//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_UTILS_HPP
#define OPENVDBBRIDGE_UTILS_HPP

#include <string>

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
    return std::find(begin, end, option) != end;
}

char *getCmdOption(char **begin, char **end, const std::string &option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) return *itr;
    return nullptr;
}

#endif //OPENVDBBRIDGE_UTILS_HPP
