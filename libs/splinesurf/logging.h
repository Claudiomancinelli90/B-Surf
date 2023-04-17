#pragma once

#include <yocto/yocto_math.h>
using namespace yocto;

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace logging {
using namespace std;

#define FORMAT_TOKEN '$'
#define FORMAT_TOKEN_LENGTH 1

// clang-format off
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */

#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
// clang-format on

// Iostream utilities for basic types
inline ostream& operator<<(ostream& os, const vec2f& v) {
  return os << "vec2f"
            << "[" << v.x << ", " << v.y << "]";
}
inline ostream& operator<<(ostream& os, const vec3f& v) {
  return os << "vec3f"
            << "[" << v.x << ", " << v.y << ", " << v.z << "]";
}
inline ostream& operator<<(ostream& os, const vec4f& v) {
  return os << "vec4f"
            << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
}
inline ostream& operator<<(ostream& os, const vec2i& v) {
  return os << "vec2i"
            << "[" << v.x << ", " << v.y << "]";
}
inline ostream& operator<<(ostream& os, const vec3i& v) {
  return os << "vec3i"
            << "[" << v.x << ", " << v.y << ", " << v.z << "]";
}
inline ostream& operator<<(ostream& os, const vec4i& v) {
  return os << "vec4i"
            << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
}
inline ostream& operator<<(ostream& os, const mat2f& v) {
  return os << "mat2f"
            << "[" << v.x << ", " << v.y << "]";
}
inline ostream& operator<<(ostream& os, const mat3f& v) {
  return os << "mat3f"
            << "[" << v.x << ", " << v.y << ", " << v.z << "]";
}
inline ostream& operator<<(ostream& os, const mat4f& v) {
  return os << "mat4f"
            << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
}
inline ostream& operator<<(ostream& os, const frame2f& v) {
  return os << "frame2f"
            << "[" << v.x << ", " << v.y << ", " << v.o << "]";
}
inline ostream& operator<<(ostream& os, const frame3f& v) {
  return os << "frame3f"
            << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.o << "]";
}

enum struct log_unit {
  nanosecond,
  millisecond,
  microsecond,
  second,
  minute,
  hour,
  full
};

inline int64_t get_time() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Print a value
template <typename T>
inline bool print_value(stringstream& stream, const T& v) {
  stream << v;
  return (bool)stream;
}

// Print a value
template <typename T>
inline void write(const T& v) {
  cout << v;
}

// Prints a string.
inline bool format_next(stringstream& stream, const string& fmt) {
  return print_value(stream, fmt);
}
template <typename Arg, typename... Args>
inline bool format_next(stringstream& stream, const string& fmt, const Arg& arg,
    const Args&... args) {
  auto pos = fmt.find(FORMAT_TOKEN);
  if (pos == string::npos) return print_value(stream, fmt);
  if (!print_value(stream, fmt.substr(0, pos))) return false;
  if (!print_value(stream, arg)) return false;
  return format_next(stream, fmt.substr(pos + FORMAT_TOKEN_LENGTH), args...);
}

// Formats a string `fmt` with values taken from `args`. Uses FORMAT_TOKEN as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args) {
  auto stream = stringstream();
  format_next(stream, fmt, args...);
  return stream.str();
}

// log levels
enum struct log_level { error = 0, check = 1, info = 2, time = 3, extra = 4 };

// Logging configutation
inline bool& _log_console() {
  static auto _log_console = true;
  return _log_console;
}

inline bool& _log_colors() {
  static auto _log_colors = true;
  return _log_colors;
}
inline ofstream& _log_filestream() {
  static auto _log_filestream = ofstream();
  return _log_filestream;
}
// inline log_level& _log_level() {
//     static auto _log_level = log_level::time;
//     return _log_level;
// }
inline bool& _log_time() {
  static bool log_time = true;
  return log_time;
}
inline bool& _log_info() {
  static bool log_info = true;
  return log_info;
}
inline bool& _log_check() {
  static bool log_check = true;
  return log_check;
}
inline bool& _log_error() {
  static bool log_error = true;
  return log_error;
}
inline bool& _log_quiet() {
  static bool log_quiet = false;
  return log_quiet;
}
inline void set_log_time(bool enabled) { _log_time() = enabled; }
inline void set_log_info(bool enabled) { _log_info() = enabled; }
inline void set_log_check(bool enabled) { _log_check() = enabled; }
inline void set_log_error(bool enabled) { _log_error() = enabled; }
inline void set_log_quiet(bool enabled) { _log_quiet() = enabled; }
inline bool is_log_level_skipped(log_level level) {
  if (_log_quiet()) return true;
  switch (level) {
    case log_level::info: return !_log_info();
    case log_level::time: return !_log_time();
    case log_level::check: return !_log_check();
    case log_level::error: return !_log_error();
    case log_level::extra: return false;
  }
}
inline int& _log_nesting() {
  static int _log_nesting = 0;
  return _log_nesting;
}
inline int& _log_last_nesting() {
  static int _log_last_nesting = 0;
  return _log_last_nesting;
}

// struct StringBuilder {
//     struct slice {
//         int count, next;
//     };
//     vector<char>  buffer;
//     vector<slice> slices;

//     StringBuilder(size_t capacity = 1024) {
//         buffer.reserve(capacity);
//         slices = {0, 0};
//     }
// };

// inline void operator+=(StringBuilder& buffer, const char* s) {
//     auto& a   = buffer.buffer;
//     auto  len = strlen(s);
//     if (buffer.strings.size() == 1) {
//         buffer.strings[0].count += len;
//     } else {
//         buffer.string.push_back({a.size(), len});
//     }
//     a.insert(a.end(), s, s + len);
// }

// inline void write(const StringBuilder& buffer) {
//     for (int i = buffer.strings.size() - 1; i >= 0; i--) {
//         auto& string = buffer.strings[i];
//         printf("%.*s\n", string.count, buffer.buffer[string.start]);
//     }
// }

template <int capacity>
struct Text {
  char data[capacity];
  int  count = 0;

  void operator+=(const char* s) {
    auto len = strlen(s);
    assert(count + len < capacity);
    memcpy(data + count, s, len);
    count += len;
  }
};

template <int N>
inline void write(const Text<N>& text) {
  printf("%.*s", text.count, text.data);
}

// Logs a message
inline void log_message(
    log_level level, const char* msg, bool new_line = true, int nesting = -1) {
  static const char* log_tags[]   = {"error", "check", "info", "time"};
  static const char* log_colors[] = {MAGENTA, RED, YELLOW, BLUE};
  static bool        ever_logged  = false;
  if (nesting == -1) nesting = _log_nesting();
  int i = int(level);
  if (_log_console()) {
    auto out = Text<255>{};
    if (new_line) {
      if (ever_logged) out += "\n";
      if (i > 1) out += " ";
      if (i != 4) {
        out += "[";
        if (_log_colors()) out += log_colors[i];
        out += log_tags[i];
        if (_log_colors()) out += RESET;
        out += "] ";
      } else {
        out += "       ";
      }
      if (nesting) {
        for (int i = 0; i < nesting; i++) out += "| ";
      }
    }
    out += msg;
    write(out);
    fflush(stdout);
  }
  if (_log_filestream()) {
    _log_filestream() << log_tags[i] << ", " << msg << "\n";
  }
  ever_logged         = true;
  _log_last_nesting() = nesting;
}

// Log info/check/error/time
template <typename... Args>
inline void info(const string& fmt, const Args&... args) {
  if (is_log_level_skipped(log_level::info)) return;
  log_message(log_level::info, format(fmt, args...).c_str());
}

struct Info {
  int nesting;
  inline ~Info() {
    if (nesting != _log_last_nesting()) {
      log_message(log_level::info, "*", true, nesting);
    } else {
    }
    _log_nesting() -= 1;
  }
};

template <typename... Args>
inline Info info_scope(const string& fmt, const Args&... args) {
  info(fmt, args...);
  _log_nesting() += 1;
  return {_log_nesting() - 1};
}

// template <typename... Args>
// inline void check(const string& fmt, const Args&... args) {
//     if (is_log_level_skipped(log_level::check)) return;
//     log_message(log_level::check, format(fmt, args...).c_str());
//     assert(0);
// }

template <typename Type>
inline void _check(const string& name, const Type& flag, const string& msg) {
  // #ifndef NDEBUG
  if (is_log_level_skipped(log_level::check)) return;
  // if (!flag) log_message(log_level::check, format(msg, args...).c_str());
  if (!flag) {
    log_message(log_level::check, format("$ == $", name, msg).c_str());
    // assert(0);
  }
  // #endif
}

template <typename A, typename B>
inline void _check_eq(const string& lname, const string& rname, const A& left,
    const B& right, const string& msg, const char* function) {
  // #ifndef NDEBUG
  if (is_log_level_skipped(log_level::check)) return;
  // if (!flag) log_message(log_level::check, format(msg, args...).c_str());
  if (left != right) {
    log_message(log_level::check, format("check failed $()", function).c_str());
    log_message(
        log_level::extra, format("  $ != $, $", lname, rname, msg).c_str());
    log_message(log_level::extra, format("  $ != $", left, right).c_str());
    // assert(0);
  }
}
template <typename Type>
inline void _check_neq(const string& lname, const string& rname,
    const Type& left, const Type& right, const string& msg,
    const char* function) {
  // #ifndef NDEBUG
  if (is_log_level_skipped(log_level::check)) return;
  // if (!flag) log_message(log_level::check, format(msg, args...).c_str());
  if (left == right) {
    log_message(log_level::check, format("check failed $()", function).c_str());
    log_message(
        log_level::extra, format("  $ == $, $", lname, rname, msg).c_str());
    log_message(log_level::extra, format("  $ == $", left, right).c_str());
    // assert(0);
  }
}

#define check(EX, msg) _check(#EX, EX, msg);
#define check_eq(L, R, msg) _check_eq(#L, #R, L, R, msg, __FUNCTION__);
#define check_neq(L, R, msg) _check_neq(#L, #R, L, R, msg, __FUNCTION__);

template <typename... Args>
inline void error(const string& fmt, const Args&... args) {
  if (is_log_level_skipped(log_level::error)) return;
  log_message(log_level::error, format(fmt, args...).c_str());
  exit(1);
}

// Log traces for time and program debugging
struct Time {
  string   message    = "";
  int64_t  start_time = -1;
  bool     scoped     = false;
  log_unit unit       = log_unit::second;
  int      nesting    = 0;  // TODO: unused?

  inline ~Time();
};

inline string convert_nanoseconds(int64_t t, log_unit unit) {
  auto next_digit = [](int64_t& n) {
    int digit = n % 10;
    n /= 10;
    return digit;
  };
  int64_t ns = t;
  int64_t µs = t / 1000;
  int64_t ms = µs / 1000;
  int64_t s  = ms / 1000;
  int64_t m  = s / 60;
  int64_t h  = m / 60;
  if (unit == log_unit::nanosecond) {
    return format("$ns", ns);
  }
  if (unit == log_unit::microsecond) {
    auto n0 = next_digit(ns);
    auto n1 = next_digit(ns);
    auto n2 = next_digit(ns);
    return format("$.$$$ms", µs, n2, n1, n0);
  }
  if (unit == log_unit::millisecond) {
    auto µ0 = next_digit(µs);
    auto µ1 = next_digit(µs);
    auto µ2 = next_digit(µs);
    return format("$.$$$ms", ms, µ2, µ1, µ0);
  }
  if (unit == log_unit::second) {
    auto ms0 = next_digit(ms);
    auto ms1 = next_digit(ms);
    auto ms2 = next_digit(ms);
    return format("$.$$$s", s, ms2, ms1, ms0);
  }
  if (unit == log_unit::minute) {
    auto s0 = next_digit(s);
    auto s1 = next_digit(s);
    return format("$m:$$s", m, s1, s0);
  }
  if (unit == log_unit::hour) {
    auto m1 = next_digit(m);
    auto m0 = next_digit(m);
    return format("$h:$$m", h, m1, m0);
  }
  if (unit == log_unit::full) {
    assert(0);
  }
  return "ERROR: unknown conversion unit, cannot convert nanosecond";
}

inline void log_time(const string& message, int nesting = 0) {
  if (is_log_level_skipped(log_level::time)) return;
  log_message(log_level::time, message.c_str(), true, nesting);
}

inline Time time_begin(const string& message, bool show_start = true,
    log_unit unit = log_unit::second, bool scoped = false) {
  if (is_log_level_skipped(log_level::time)) return Time();
  // if (show_start)
  log_message(log_level::time, message.c_str(), true, _log_nesting());
  _log_nesting() += 1;
  return Time{message, get_time(), scoped, unit, _log_nesting() - 1};
}

inline void time_end(Time& scope) {
  if (is_log_level_skipped(log_level::time)) return;
  if (scope.start_time >= 0) {
    auto time = convert_nanoseconds(get_time() - scope.start_time, scope.unit);
    if (scope.nesting != _log_last_nesting()) {
      log_message(log_level::time, (scope.message + ": " + time).c_str(), true,
          scope.nesting);
    } else {
      log_message(log_level::time, (": " + time).c_str(), false, scope.nesting);
    }
    // log_time(time, scope.nesting);
  } else {
    log_message(log_level::time, (scope.message + " ended").c_str(), true,
        scope.nesting);
  }
  _log_nesting() -= 1;
}

inline Time time_scope(
    const string& message, log_unit unit = log_unit::second) {
  return time_begin(message, true, unit, true);
}

inline Time::~Time() {
  if (scoped) time_end(*this);
}

// Configure the logging
inline void set_log_colors(bool enabled) { _log_colors() = enabled; }
inline void set_log_console(bool enabled) { _log_console() = enabled; }
inline void set_log_file(const string& filename, bool append) {
  if (_log_filestream()) {
    _log_filestream().close();
    _log_filestream() = {};
  }
  if (empty(filename)) return;
  _log_filestream().open(filename, append ? std::ios::app : std::ios::out);
}
}  // namespace logging

#define info_function() \
  auto _info = logging::info_scope("$$", __FUNCTION__, "()");
#define time_function() \
  auto _time = logging::time_scope(string(__FUNCTION__) + "()");
#define here() logging::info("here: $ - $", __FUNCTION__, __LINE__);
#define show(V) logging::info(string(#V) + ": $", V);
