#ifndef _REALTIME_GUI_
#define _REALTIME_GUI_

#include "window.h"

namespace gui {
using namespace yocto;
using namespace window;

// Tools for building a user interface. Maybe put this in another file?
void init_gui(Window& win, int width, bool left = true);
void gui_begin(const Window& win, const char* name);
void gui_end(const Window& win);

bool is_gui_active();
bool begin_header(Window& win, const char* title);
void end_header(Window& win);
void gui_label(Window& win, const char* lbl, const std::string& text);

void gui_text(Window& win, const string& text);
void gui_text(Window& win, const char* format, ...);
void gui_separator(Window& win);
void continue_line(Window& win);

bool gui_button(Window& win, const char* lbl, bool enabled = true);

bool gui_textinput(Window& win, const char* lbl, std::string& value);

bool gui_slider(
    Window& win, const char* lbl, float& value, float min, float max);
bool gui_slider(
    Window& win, const char* lbl, vec2f& value, float min, float max);
bool gui_slider(
    Window& win, const char* lbl, vec3f& value, float min, float max);
bool gui_slider(
    Window& win, const char* lbl, vec4f& value, float min, float max);

bool gui_slider(Window& win, const char* lbl, int& value, int min, int max);
bool gui_slider(Window& win, const char* lbl, vec2i& value, int min, int max);
bool gui_slider(Window& win, const char* lbl, vec3i& value, int min, int max);
bool gui_slider(Window& win, const char* lbl, vec4i& value, int min, int max);

bool gui_dragger(Window& win, const char* lbl, float& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool gui_dragger(Window& win, const char* lbl, vec2f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool gui_dragger(Window& win, const char* lbl, vec3f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool gui_dragger(Window& win, const char* lbl, vec4f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);

bool gui_dragger(Window& win, const char* lbl, int& value, float speed = 1,
    int min = 0, int max = 0);
bool gui_dragger(Window& win, const char* lbl, vec2i& value, float speed = 1,
    int min = 0, int max = 0);
bool gui_dragger(Window& win, const char* lbl, vec3i& value, float speed = 1,
    int min = 0, int max = 0);
bool gui_dragger(Window& win, const char* lbl, vec4i& value, float speed = 1,
    int min = 0, int max = 0);

bool gui_checkbox(Window& win, const char* lbl, bool& value);
bool gui_checkbox(Window& win, const string& lbl);

bool gui_coloredit(Window& win, const char* lbl, vec3f& value);
bool gui_coloredit(Window& win, const char* lbl, vec4f& value);

bool gui_hdrcoloredit(Window& win, const char* lbl, vec3f& value);
bool gui_hdrcoloredit(Window& win, const char* lbl, vec4f& value);

bool gui_combobox(Window& win, const char* lbl, int& idx,
    const std::vector<std::string>& labels);
bool gui_combobox(Window& win, const char* lbl, std::string& value,
    const std::vector<std::string>& labels);
bool gui_combobox(Window& win, const char* lbl, int& idx, int num,
    const std::function<std::string(int)>& labels, bool include_null = false);

template <typename T>
inline bool gui_combobox(Window& win, const char* lbl, T*& value,
    const std::vector<T*>& vals, bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = gui_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return vals[idx]->name; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

template <typename T>
inline bool gui_combobox(Window& win, const char* lbl, T*& value,
    const std::vector<T*>& vals, const std::vector<std::string>& labels,
    bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = gui_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return labels[idx]; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

void gui_progressbar(Window& win, const char* lbl, float fraction);

void gui_histogram(
    Window& win, const char* lbl, const std::vector<float>& values);
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec2f>& values);
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec3f>& values);
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec4f>& values);

bool gui_messages(Window& win);
void push_message(Window& win, const std::string& message);
bool gui_filedialog(Window& win, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter);
bool gui_filedialog_button(Window& win, const char* button_lbl,
    bool button_active, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter);

void log_info(Window& win, const std::string& msg);
void log_error(Window& win, const std::string& msg);
void clear_log(Window& win);
void gui_log(Window& win);
}  // namespace gui

#endif