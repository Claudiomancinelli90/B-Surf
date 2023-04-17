#ifndef _REALTIME_WINDOW_
#define _REALTIME_WINDOW_

#include <functional>
using std::function;

// #include <graphics/common.h>
// #include <graphics/math.h>

#include <yocto/yocto_common.h>
#include <yocto/yocto_math.h>

#include <array>
#include <unordered_map>

// Forward declaration
struct GLFWwindow;

namespace window {
using namespace yocto;

enum struct Key : int {
  // For printable keys, just use the constructor, like Key('*').
  // For letters, always use upper case, like Key('C').

  escape        = 256,
  enter         = 257,
  tab           = 258,
  backspace     = 259,
  insert        = 260,
  _delete       = 261,
  right         = 262,
  left          = 263,
  down          = 264,
  up            = 265,
  page_up       = 266,
  page_down     = 267,
  home          = 268,
  end           = 269,
  caps_lock     = 280,
  scroll_lock   = 281,
  num_lock      = 282,
  print_screen  = 283,
  pause         = 284,
  f1            = 290,
  f2            = 291,
  f3            = 292,
  f4            = 293,
  f5            = 294,
  f6            = 295,
  f7            = 296,
  f8            = 297,
  f9            = 298,
  f10           = 299,
  f11           = 300,
  f12           = 301,
  f13           = 302,
  f14           = 303,
  f15           = 304,
  f16           = 305,
  f17           = 306,
  f18           = 307,
  f19           = 308,
  f20           = 309,
  f21           = 310,
  f22           = 311,
  f23           = 312,
  f24           = 313,
  f25           = 314,
  kp_0          = 320,
  kp_1          = 321,
  kp_2          = 322,
  kp_3          = 323,
  kp_4          = 324,
  kp_5          = 325,
  kp_6          = 326,
  kp_7          = 327,
  kp_8          = 328,
  kp_9          = 329,
  kp_decimal    = 330,
  kp_divide     = 331,
  kp_multiply   = 332,
  kp_subtract   = 333,
  kp_add        = 334,
  kp_enter      = 335,
  kp_equal      = 336,
  left_shift    = 340,
  left_control  = 341,
  left_alt      = 342,
  left_super    = 343,
  right_shift   = 344,
  right_control = 345,
  right_alt     = 346,
  right_super   = 347,
  menu          = 348,
  world_1       = 161,  //  non-us #1
  world_2       = 162   //  non-us #2
};

struct Button {
  enum struct State : unsigned char {
    nothing = 0,  // button is not pressed (!holding...)
    pressing,     // button is being pressed during this frame.
    holding,      // button is down, mmmhhh
    releasing,    // button is being pressed released
  };
  enum struct Modifier : unsigned char {
    none      = 0,
    shift     = 0x0001,
    control   = 0x0002,
    alt       = 0x0004,
    super     = 0x0008,
    caps_lock = 0x0010,
    num_lock  = 0x0020,
  };
  State    state    = State::nothing;
  Modifier modifier = Modifier::none;

  operator bool() const { return state == Button::State::holding; }
};

// Data structure where all the inputs are stored between frames.
struct Input {
  // Ever-changing data
  vec2f    mouse_pos  = {0, 0};  // position excluding gui
  vec2f    mouse_last = {0, 0};  // last mouse position excluding gui
  uint64_t clock_now  = 0;       // clock now
  uint64_t clock_last = 0;       // clock last
  double   time_now   = 0;       // time now
  double   time_delta = 0;       // time delta
  int      frame      = 0;

  Button click_left   = {};
  Button click_middle = {};
  Button click_right  = {};
  vec2f  scroll       = {0, 0};  // scroll input

  bool modifier_alt   = false;  // alt modifier
  bool modifier_ctrl  = false;  // ctrl modifier
  bool modifier_shift = false;  // shift modifier

  vec2i window_size          = {0, 0};
  vec4i framebuffer_viewport = {0, 0, 0, 0};
  vec2i framebuffer_size     = {0, 0};
  bool  is_window_focused    = false;  // window is focused
  bool  is_gui_active        = false;  // gui is active

  std::array<Button, 512> key_buttons = {};
};

// Data structure where the input of a joystick is stored.
struct Joystick {
  vec2f left_stick, right_stick;
  float left_trigger, right_trigger;
  bool  buttons[15];
  int   id = -1;

  // Handy interface for buttons. The mapping may be broken on some platform.
  // Using the recommanded mapping from GLFW gave wrong results with a PS4
  // Dualshock on a Mac. Change the implementation of these methods if needed.
  bool cross() const;
  bool circle() const;
  bool triangle() const;
  bool square() const;
  bool left_bumper() const;
  bool right_bumper() const;
  bool start() const;
  bool back() const;
  bool guide() const;
  bool A() const;
  bool B() const;
  bool X() const;
  bool Y() const;
};

// Info to open and handle a new window within the OS.
struct Window {
  string      title = "";
  GLFWwindow* glfw  = nullptr;
  Input       input = {};
  // vector<Joystick> joysticks = {};

  // TODO(giacomo): cleanup
  vec2i                                 size           = {0, 0};
  int                                   gui_width      = 0;
  bool                                  gui_left       = true;
  std::unordered_map<std::string, bool> gui_buttons    = {};
  std::unordered_map<std::string, bool> gui_checkboxes = {};
};

void init_window(Window& win, const vec2i& size, const string& title,
    bool visible = true, int msaa = 1);
void delete_window(Window& win);

bool should_window_close(const Window& win);

// TODO(giaocomo): remove
vec2f get_mouse_pos_normalized(const Window& win, bool isometric = false);
vec2f get_mouse_pos_normalized(const Input& input, bool isometric = false);

bool is_pressing(const Input& input, Key key);
bool is_down(const Input& input, Key key);
bool is_releasing(const Input& input, Key key);
bool is_pressing(Button button);
bool is_down(Button button);
bool is_releasing(Button button);

void update_window_size(Window& win);
void update_input(Input& input, const Window& win);
void update_input(Window& win);
// void init_joysticks(vector<Joystick>& joysticks);
// void update_joystick_input(vector<Joystick>& joysticks);
void init_callbacks(Window& win);

void poll_events(const Window& win, bool wait);
void swap_buffers(const Window& win);
void run_draw_loop(
    Window& win, function<void(const Input&)> draw, bool wait = false);

void update_camera_frame(frame3f& frame, float& focus, const Input& win,
    bool rotating, bool panning, bool zooming);
}  // namespace window

#endif
