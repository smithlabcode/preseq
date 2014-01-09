/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPTION_PARSER_HPP
#define OPTION_PARSER_HPP

#include <string>
#include <vector>
#include "smithlab_utils.hpp"

class Option {
public:
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, int &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, unsigned &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, long &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, unsigned long &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, float &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, double &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, std::string &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, bool &val);
  Option(const std::string l_name, const char s_name, 
	 const std::string descr,  const bool reqd, char &val);
  
  bool parse(std::vector<std::string> &command_line);
  void parse_config_file(std::vector<std::string> &options);

  std::string format_option_name() const;
  std::string format_option_description(const size_t offset) const;
  
private:
  
  unsigned arg_type;
  std::string long_name;
  char short_name;
  std::string description;
  bool required;
  bool specified;
  
  // the values of the options: ugly but clean
  int *int_value;
  unsigned int *uint_value;
  long *long_value;
  unsigned long *ulong_value;
  float *float_value;
  double *double_value;
  std::string *string_value;
  bool *bool_value;    
  char *char_value;    
  
  void format_option(const std::string &argument);
  static void set_max_length(size_t num);
  static size_t get_max_length();
  bool option_match(const std::string &other);
};

struct SMITHLABOptionException : public SMITHLABException {
  SMITHLABOptionException(std::string m) : SMITHLABException(m) {}
};

class OptionParser {
public:

  OptionParser(const std::string nm, const std::string descr,
	       std::string noflag_msg = "");
  
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, int &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, unsigned &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, long &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, unsigned long &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, float &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, double &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, std::string &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, bool &val);
  void add_opt(const std::string l_name, const char s_name, 
	       const std::string descr,  const bool reqd, char &val);

  void parse(const int argc, const char **argv,
             std::vector<std::string> &arguments);

  void parse(const int argc, const char **argv, 
	     std::vector<std::string> &arguments,
	     std::string config_filename);

  bool help_requested() const {return help_request;}
  std::string help_message() const;
  
  bool about_requested() const {return about_request;}
  std::string about_message() const;

  bool option_missing() const {
    return !first_missing_option_name.empty();
  }
  std::string option_missing_message() const;

  static const bool OPTIONAL = false;
  static const bool REQUIRED = true;
    
private:
  std::string prog_name;
  std::string prog_descr;
  std::string noflag_message;
  std::vector<Option> options;

  bool help_request;
  bool about_request;
  std::string first_missing_option_name;
};

#endif
