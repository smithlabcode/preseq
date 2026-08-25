/* Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef SRC_DNMT_ERROR_HPP_
#define SRC_DNMT_ERROR_HPP_

#include <cstdint>  // for int64_t
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

struct dnmt_error : public std::exception {
  std::int64_t err{};    // error possibly from HTSlib
  int the_errno{};       // ERRNO at time of construction
  std::string msg;       // the message
  std::string the_what;  // to report
  dnmt_error(const std::int64_t err, const std::string &msg) :
    err{err}, the_errno{errno}, msg{msg} {
    std::ostringstream oss;
    // clang-format off
    oss << "[error: " << err << "][" << "ERRNO: " << the_errno << "]"
        << "[" << strerror(the_errno) << "][" << msg << "]";
    // clang-format on
    the_what = oss.str();
  }
  explicit dnmt_error(const std::string &_msg) : dnmt_error(0, _msg) {}
  const char *what() const noexcept override { return the_what.data(); }
};

#endif  // SRC_DNMT_ERROR_HPP_
