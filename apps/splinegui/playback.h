#pragma once
#include <realtime/window.h>

#include "serialize/serialize.h"
using namespace window;

inline void serialize_input(Serializer& serializer, gui_input& input) {
  serialize(serializer, input);
}

#if 0
inline void save_input_record(
    vector<gui_input>& record, const std::string& filename) {
  auto serializer = make_writer(filename, 1024);
  serialize_vector(serializer, record);
  close_serializer(serializer);
}

inline void load_input_record(
    vector<gui_input>& record, const std::string& filename) {
  auto serializer = make_reader(filename, 1024);
  serialize_vector(serializer, record);
  close_serializer(serializer);
}
#else

struct Compressed_Input {
  struct Input_Diff {
    unsigned int position;
    unsigned int diff;
  };
  vector<Input_Diff> diffs = {};
};

inline Compressed_Input compress_input(const gui_input& input) {
  static const auto def_input = gui_input{};

  auto result = Compressed_Input{};

  auto ptr_def = (unsigned int*)&def_input;
  auto ptr_in  = (unsigned int*)&input;
  for (unsigned int i = 0; i < sizeof(def_input) / sizeof(unsigned int); i++) {
    if (ptr_in[i] != ptr_def[i]) {
      result.diffs.push_back({i, ptr_in[i] ^ ptr_def[i]});
    }
  }
  return result;
}
inline gui_input uncompress_input(const Compressed_Input& cinput) {
  auto result = gui_input{};

  auto ptr = (unsigned int*)&result;
  for (auto& diff : cinput.diffs) {
    ptr[diff.position] ^= diff.diff;
  }
  return result;
}

inline void serialize_compressed_input(
    Serializer& serializer, Compressed_Input& cinput) {
  serialize_vector(serializer, cinput.diffs);
}

inline void save_input_record(
    vector<gui_input>& record, const std::string& filename) {
  auto serializer = make_writer(filename, 1024);
  auto xxx        = vector<Compressed_Input>(record.size());
  for (int i = 0; i < record.size(); i++) {
    xxx[i] = compress_input(record[i]);
  }
  serialize_vector_custom(serializer, xxx, serialize_compressed_input);
  close_serializer(serializer);
}

inline void load_input_record(
    vector<gui_input>& record, const std::string& filename) {
  auto serializer = make_reader(filename, 1024);
  auto xxx        = vector<Compressed_Input>();
  serialize_vector_custom(serializer, xxx, serialize_compressed_input);
  record.resize(xxx.size());
  for (int i = 0; i < record.size(); i++) {
    record[i] = uncompress_input(xxx[i]);
  }
  close_serializer(serializer);
}
#endif
