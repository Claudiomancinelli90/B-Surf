#pragma once

#include <yocto/yocto_common.h>

#include <unordered_map>

using ProfileEntry = string;

inline std::unordered_map<ProfileEntry, float>& global_times() {
  static std::unordered_map<ProfileEntry, float> _times;
  return _times;
}

struct ProfileNode {
  string name = "void";
  float  time = 0;
  // float                delay    = 0;
  ProfileNode*         parent   = nullptr;
  vector<ProfileNode*> children = {};
};

inline ProfileNode* root_node() {
  static ProfileNode* _root = new ProfileNode();
  return _root;
}

inline ProfileNode*& current_node() {
  static ProfileNode* _parent = root_node();
  return _parent;
}

static uint64_t clock_delay = 0;
static uint64_t last_clock  = 0;

inline void stop_clock() { last_clock = get_time(); }

inline void start_clock() { clock_delay += get_time() - last_clock; }

inline uint64_t get_clock() { return get_time() - clock_delay; }

struct ScopeProfiler {
  string key;
  float  start;

  ScopeProfiler(const string& name) {
    // auto delay = get_time();

    auto parent = current_node();
    int  it     = -1;
    for (int i = 0; i < parent->children.size(); ++i) {
      if (parent->children[i]->name == name) {
        it = i;
        break;
      }
    }
    if (it == -1) {
      it            = parent->children.size();
      auto child    = new ProfileNode();
      child->name   = name;
      child->parent = parent;
      parent->children.push_back(child);
    }
    current_node() = parent->children[it];

    auto time = get_time();
    // current_node()->delay += (time - delay) / 1e6;
    start = time;
  }

  ~ScopeProfiler() {
    float t = get_time() - start;
    current_node()->time += t / 1e6;
    current_node() = current_node()->parent;
  }
};

#if 0
#define profile_scope(name) auto _timer = ScopeProfiler(name);
#define profile_function() auto _timer = ScopeProfiler(__FUNCTION__);
#else
#define profile_scope(name) ;
#define profile_function() ;
#endif

inline void print_node(const ProfileNode* node, int depth) {
  for (int i = 0; i < node->children.size(); ++i) {
    printf("%8d", int(node->children[i]->time));
    printf(": %.*s", 2 * depth, "| | | | | | | | | | | | ");

    printf("%s", node->children[i]->name.c_str());
    float mis = 0;
    for (auto& c : node->children[i]->children) mis += c->time;
    if (mis)
      printf(" [%d%%]", int(100 * (1 - (node->children[i]->time - mis) /
                                           node->children[i]->time)));
    printf("\n");
    print_node(node->children[i], depth + 1);
  }
}

// inline void compute_delay(ProfileNode* node) {
//   for (auto child : node->children) {
//     compute_delay(child);
//   }
//   for (auto child : node->children) {
//     node->delay += child->delay;
//   }
//   for (auto child : node->children) {
//     node->time -= child->delay;
//   }
// }

inline void print_profiling_result() {
  if (current_node() != root_node()) return;
  printf("\n");
  printf("%8s: function\n", "time(ms)");
  // compute_delay(root_node());
  print_node(root_node(), 0);
}
