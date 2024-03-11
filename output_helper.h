#pragma once

#include "field.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#    if __cplusplus >= 201703L && __has_include(<filesystem>)
#        include <filesystem>
namespace fs = std::filesystem;
#    elif __has_include(<experimental/filesystem>)
#        include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#    elif __has_include(<boost/filesystem.hpp>)
#        include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#    endif
#endif

void field_output(field& _field, std::string root_dir, std::string label, int ind);