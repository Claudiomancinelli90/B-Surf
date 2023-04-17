#include "gui.h"
#include <yocto/yocto_parallel.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"
#define CUTE_FILES_IMPLEMENTATION
#include "ext/cute_files.h"

using namespace std::string_literals;

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace gui {
using namespace yocto;
using namespace window;

struct Panel {
  vec2f       position;
  vec2f       size;
  const char* name;
};

void init_gui(Window& win, int width, bool left) {
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(win.glfw, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();
  win.gui_width = width;
  win.gui_left  = left;
}

void gui_begin(const Window& win, const char* name) {
  if (win.gui_width <= 0) return;
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  if (win.gui_left) {
    ImGui::SetNextWindowPos({0, 0});
    ImGui::SetNextWindowSize({float(win.gui_width), float(win.size.y)});
  } else {
    ImGui::SetNextWindowPos({float(win.size.x - win.gui_width), 0});
    ImGui::SetNextWindowSize({float(win.gui_width), float(win.size.y)});
  }

  ImGui::SetNextWindowBgAlpha(1);
  ImGui::Begin(
      name, nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
}

void gui_end(const Window& win) {
  ImGui::End();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool is_gui_active() {
  auto io = &ImGui::GetIO();
  return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

bool begin_header(Window& win, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_header(Window& win) { ImGui::PopID(); }

void open_glmodal(Window& win, const char* lbl) { ImGui::OpenPopup(lbl); }
void clear_glmodal(Window& win) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(Window& win, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(Window& win) { ImGui::EndPopup(); }
bool is_glmodal_open(Window& win, const char* lbl) {
  return ImGui::IsPopupOpen(lbl);
}

bool gui_message(Window& win, const char* lbl, const std::string& message) {
  if (ImGui::BeginPopupModal(lbl)) {
    auto open = true;
    ImGui::Text("%s", message.c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      open = false;
    }
    ImGui::EndPopup();
    return open;
  } else {
    return false;
  }
}

std::deque<std::string> _message_queue = {};
std::mutex              _message_mutex;
void                    push_message(Window& win, const std::string& message) {
  std::lock_guard lock(_message_mutex);
  _message_queue.push_back(message);
}
bool gui_messages(Window& win) {
  std::lock_guard lock(_message_mutex);
  if (_message_queue.empty()) return false;
  if (!is_glmodal_open(win, "<message>")) {
    open_glmodal(win, "<message>");
    return true;
  } else if (ImGui::BeginPopupModal("<message>")) {
    ImGui::Text("%s", _message_queue.front().c_str());
    if (ImGui::Button("Ok")) {
      ImGui::CloseCurrentPopup();
      _message_queue.pop_front();
    }
    ImGui::EndPopup();
    return true;
  } else {
    return false;
  }
}

// Utility to normalize a path
static inline std::string normalize_path(const std::string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)

    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get extension (not including '.').
static std::string get_extension(const std::string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == std::string::npos) return "";
  return filename.substr(pos);
}

struct filedialog_state {
  std::string                               dirname       = "";
  std::string                               filename      = "";
  std::vector<std::pair<std::string, bool>> entries       = {};
  bool                                      save          = false;
  bool                                      remove_hidden = true;
  std::string                               filter        = "";
  std::vector<std::string>                  extensions    = {};

  filedialog_state() {}
  filedialog_state(const std::string& dirname, const std::string& filename,
      bool save, const std::string& filter) {
    this->save = save;
    set_filter(filter);
    set_dirname(dirname);
    set_filename(filename);
  }
  void set_dirname(const std::string& name) {
    dirname = name;
    dirname = normalize_path(dirname);
    if (dirname == "") dirname = "./";
    if (dirname.back() != '/') dirname += '/';
    refresh();
  }
  void set_filename(const std::string& name) {
    filename = name;
    check_filename();
  }
  void set_filter(const std::string& flt) {
    auto globs = std::vector<std::string>{""};
    for (auto i = 0; i < flt.size(); i++) {
      if (flt[i] == ';') {
        globs.push_back("");
      } else {
        globs.back() += flt[i];
      }
    }
    filter = "";
    extensions.clear();
    for (auto pattern : globs) {
      if (pattern == "") continue;
      auto ext = get_extension(pattern);
      if (ext != "") {
        extensions.push_back(ext);
        filter += (filter == "") ? ("*." + ext) : (";*." + ext);
      }
    }
  }
  void check_filename() {
    if (filename.empty()) return;
    auto ext = get_extension(filename);
    if (std::find(extensions.begin(), extensions.end(), ext) ==
        extensions.end()) {
      filename = "";
      return;
    }
    if (!save && !exists_file(dirname + filename)) {
      filename = "";
      return;
    }
  }
  void select_entry(int idx) {
    if (entries[idx].second) {
      set_dirname(dirname + entries[idx].first);
    } else {
      set_filename(entries[idx].first);
    }
  }

  void refresh() {
    entries.clear();
    cf_dir_t dir;
    cf_dir_open(&dir, dirname.c_str());
    while (dir.has_next) {
      cf_file_t file;
      cf_read_file(&dir, &file);
      cf_dir_next(&dir);
      if (remove_hidden && file.name[0] == '.') continue;
      if (file.is_dir) {
        entries.push_back({file.name + "/"s, true});
      } else {
        entries.push_back({file.name, false});
      }
    }
    cf_dir_close(&dir);
    std::sort(entries.begin(), entries.end(), [](auto& a, auto& b) {
      if (a.second == b.second) return a.first < b.first;
      return a.second;
    });
  }

  std::string get_path() const { return dirname + filename; }
  bool        exists_file(const std::string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
  }
};
bool gui_filedialog(Window& win, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter) {
  static auto states = std::unordered_map<std::string, filedialog_state>{};
  ImGui::SetNextWindowSize({500, 300}, ImGuiCond_FirstUseEver);
  if (ImGui::BeginPopupModal(lbl)) {
    if (states.find(lbl) == states.end()) {
      states[lbl] = filedialog_state{dirname, filename, save, filter};
    }
    auto& state = states.at(lbl);
    char  dir_buffer[1024];
    strcpy(dir_buffer, state.dirname.c_str());
    if (ImGui::InputText("dir", dir_buffer, sizeof(dir_buffer))) {
      state.set_dirname(dir_buffer);
    }
    auto current_item = -1;
    if (ImGui::ListBox(
            "entries", &current_item,
            [](void* data, int idx, const char** out_text) -> bool {
              auto& state = *(filedialog_state*)data;
              *out_text   = state.entries[idx].first.c_str();
              return true;
            },
            &state, (int)state.entries.size())) {
      state.select_entry(current_item);
    }
    char file_buffer[1024];
    strcpy(file_buffer, state.filename.c_str());
    if (ImGui::InputText("file", file_buffer, sizeof(file_buffer))) {
      state.set_filename(file_buffer);
    }
    char filter_buffer[1024];
    strcpy(filter_buffer, state.filter.c_str());
    if (ImGui::InputText("filter", filter_buffer, sizeof(filter_buffer))) {
      state.set_filter(filter_buffer);
    }
    auto ok = false, exit = false;
    if (ImGui::Button("Ok")) {
      path = state.dirname + state.filename;
      ok   = true;
      exit = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      exit = true;
    }
    if (exit) {
      ImGui::CloseCurrentPopup();
      states.erase(lbl);
    }
    ImGui::EndPopup();
    return ok;
  } else {
    return false;
  }
}
bool gui_filedialog_button(Window& win, const char* button_lbl,
    bool button_active, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter) {
  if (is_glmodal_open(win, lbl)) {
    return gui_filedialog(win, lbl, path, save, dirname, filename, filter);
  } else {
    if (gui_button(win, button_lbl, button_active)) {
      open_glmodal(win, lbl);
    }
    return false;
  }
}

bool gui_button(Window& win, const char* lbl, bool enabled) {
  if (enabled) {
    return ImGui::Button(lbl);
  } else {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    auto ok = ImGui::Button(lbl);
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
    win.gui_buttons[std::string(lbl)] = ok;
    return ok;
  }
}

void gui_label(Window& win, const char* lbl, const std::string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}

void gui_text(Window& win, const string& text) {
  ImGui::Text("%s", text.c_str());
}

void gui_text(Window& win, const char* format, ...) {
  va_list arglist;
  va_start(arglist, format);
  ImGui::TextV(format, arglist);
  va_end(arglist);
}

void gui_separator(Window& win) { ImGui::Separator(); }

void continue_line(Window& win) { ImGui::SameLine(); }

bool gui_textinput(Window& win, const char* lbl, std::string& value) {
  char buffer[4096];
  auto num = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer, sizeof(buffer));
  if (edited) value = buffer;
  return edited;
}

bool gui_slider(
    Window& win, const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool gui_slider(
    Window& win, const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool gui_slider(
    Window& win, const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool gui_slider(
    Window& win, const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool gui_slider(Window& win, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool gui_slider(Window& win, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool gui_slider(Window& win, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool gui_slider(Window& win, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool gui_dragger(Window& win, const char* lbl, float& value, float speed,
    float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool gui_dragger(Window& win, const char* lbl, vec2f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool gui_dragger(Window& win, const char* lbl, vec3f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool gui_dragger(Window& win, const char* lbl, vec4f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool gui_dragger(
    Window& win, const char* lbl, int& value, float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool gui_dragger(
    Window& win, const char* lbl, vec2i& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool gui_dragger(
    Window& win, const char* lbl, vec3i& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool gui_dragger(
    Window& win, const char* lbl, vec4i& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool gui_checkbox(Window& win, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}
bool gui_checkbox(Window& win, const string& lbl) {
  auto& value = win.gui_checkboxes[lbl];
  return ImGui::Checkbox(lbl.c_str(), &value);
}

bool gui_coloredit(Window& win, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool gui_coloredit(Window& win, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool gui_hdrcoloredit(Window& win, const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(max(color.x, color.y), color.z);
  if (scale > 1) {
    color /= scale;
    exposure = yocto::log2(scale);
  }
  auto edit_exposure = gui_slider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = gui_coloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * yocto::exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool gui_hdrcoloredit(Window& win, const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(max(color.x, color.y), color.z);
  if (scale > 1) {
    color.x /= scale;
    color.y /= scale;
    color.z /= scale;
    exposure = yocto::log2(scale);
  }
  auto edit_exposure = gui_slider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = gui_coloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value.x = color.x * yocto::exp2(exposure);
    value.y = color.y * yocto::exp2(exposure);
    value.z = color.z * yocto::exp2(exposure);
    value.w = color.w;
    return true;
  } else {
    return false;
  }
}

bool gui_combobox(Window& win, const char* lbl, int& value,
    const std::vector<std::string>& labels) {
  if (!ImGui::BeginCombo(lbl, labels[value].c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool gui_combobox(Window& win, const char* lbl, std::string& value,
    const std::vector<std::string>& labels) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool gui_combobox(Window& win, const char* lbl, int& idx, int num,
    const std::function<std::string(int)>& labels, bool include_null) {
  if (num <= 0) idx = -1;
  if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx).c_str() : "<none>"))
    return false;
  auto old_idx = idx;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", idx < 0)) idx = -1;
    if (idx < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < num; i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels(i).c_str(), idx == i)) idx = i;
    if (idx == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return idx != old_idx;
}

void gui_progressbar(Window& win, const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void gui_histogram(
    Window& win, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void gui_histogram(
    Window& win, const char* lbl, const std::vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void gui_histogram(
    Window& win, const char* lbl, const std::vector<vec4f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " w"s).c_str(), (const float*)values.data() + 3,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
}

// https://github.com/ocornut/imgui/issues/300
struct ImGuiAppLog {
  ImGuiTextBuffer Buf;
  ImGuiTextFilter Filter;
  ImVector<int>   LineOffsets;  // Index to lines offset
  bool            ScrollToBottom;

  void Clear() {
    Buf.clear();
    LineOffsets.clear();
  }

  void AddLog(const char* msg, const char* lbl) {
    int old_size = Buf.size();
    Buf.appendf("[%s] %s\n", lbl, msg);
    for (int new_size = Buf.size(); old_size < new_size; old_size++)
      if (Buf[old_size] == '\n') LineOffsets.push_back(old_size);
    ScrollToBottom = true;
  }

  void Draw() {
    if (ImGui::Button("Clear")) Clear();
    ImGui::SameLine();
    bool copy = ImGui::Button("Copy");
    ImGui::SameLine();
    Filter.Draw("Filter", -100.0f);
    ImGui::Separator();
    ImGui::BeginChild("scrolling");
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
    if (copy) ImGui::LogToClipboard();

    if (Filter.IsActive()) {
      const char* buf_begin = Buf.begin();
      const char* line      = buf_begin;
      for (int line_no = 0; line != NULL; line_no++) {
        const char* line_end = (line_no < LineOffsets.Size)
                                   ? buf_begin + LineOffsets[line_no]
                                   : NULL;
        if (Filter.PassFilter(line, line_end))
          ImGui::TextUnformatted(line, line_end);
        line = line_end && line_end[1] ? line_end + 1 : NULL;
      }
    } else {
      ImGui::TextUnformatted(Buf.begin());
    }

    if (ScrollToBottom) ImGui::SetScrollHere(1.0f);
    ScrollToBottom = false;
    ImGui::PopStyleVar();
    ImGui::EndChild();
  }
  void Draw(const char* title, bool* p_opened = NULL) {
    ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(title, p_opened);
    Draw();
    ImGui::End();
  }
};

std::mutex  _log_mutex;
ImGuiAppLog _log_widget;
void        log_info(Window& win, const std::string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_error(Window& win, const std::string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_log(Window& win) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void gui_log(Window& win) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

}  // namespace gui
