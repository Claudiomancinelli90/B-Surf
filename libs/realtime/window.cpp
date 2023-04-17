#include "window.h"

#include <GLFW/glfw3.h>

#include <cassert>
#include <unordered_map>

namespace window {

void update_camera_frame(frame3f& frame, float& focus, const Input& input,
    bool rotating, bool panning, bool zooming) {
  auto last_pos    = input.mouse_last;
  auto mouse_pos   = input.mouse_pos;
  auto mouse_left  = is_down(input.click_left);
  auto mouse_right = is_down(input.click_right);

  // handle mouse and keyboard for navigation
  if (mouse_left) {
    auto dolly  = 0.0f;
    auto pan    = zero2f;
    auto rotate = zero2f;
    if (rotating) {
      if (mouse_left) rotate = (mouse_pos - last_pos) / 100.0f;
    }
    if (zooming) {
      if (mouse_right) dolly = (mouse_pos.y - last_pos.y) / 100.0f;
    }
    if (panning) {
      if (mouse_left) pan = (mouse_pos - last_pos) * focus / 200.0f;
    }
    pan.x    = -pan.x;
    rotate.y = -rotate.y;
    update_turntable(frame, focus, rotate, dolly, pan);
  }
}

void init_joysticks(vector<Joystick>& joysticks) {
  joysticks.clear();
  const int max_joysticks = 8;
  for (int i = 0; i < max_joysticks; i++) {
    auto present    = glfwJoystickPresent(i);
    auto is_gamepad = glfwJoystickIsGamepad(i);

    if (present && is_gamepad) {
      auto& joystick = joysticks.emplace_back();
      joystick.id    = i;
    }
  }
}

// clang-format off
bool Joystick::cross() const { return buttons[GLFW_GAMEPAD_BUTTON_CIRCLE]; }
bool Joystick::circle() const { return buttons[GLFW_GAMEPAD_BUTTON_SQUARE]; }
bool Joystick::triangle() const { return buttons[GLFW_GAMEPAD_BUTTON_TRIANGLE]; }
bool Joystick::square() const { return buttons[GLFW_GAMEPAD_BUTTON_CROSS]; }
bool Joystick::left_bumper() const { return buttons[GLFW_GAMEPAD_BUTTON_LEFT_BUMPER]; }
bool Joystick::right_bumper() const { return buttons[GLFW_GAMEPAD_BUTTON_RIGHT_BUMPER]; }
bool Joystick::start() const { return buttons[GLFW_GAMEPAD_BUTTON_START]; }
bool Joystick::back() const { return buttons[GLFW_GAMEPAD_BUTTON_BACK]; }
bool Joystick::guide() const { return buttons[GLFW_GAMEPAD_BUTTON_GUIDE]; }
bool Joystick::A() const { return buttons[GLFW_GAMEPAD_BUTTON_A]; }
bool Joystick::B() const { return buttons[GLFW_GAMEPAD_BUTTON_B]; }
bool Joystick::X() const { return buttons[GLFW_GAMEPAD_BUTTON_X]; }
bool Joystick::Y() const { return buttons[GLFW_GAMEPAD_BUTTON_Y]; }
// clang-format on

inline void update_button_from_input(
    Button& button, bool pressed, Button::Modifier modifier) {
  if (pressed) {
    assert(button.state != Button::State::holding);
    button.state    = Button::State::pressing;
    button.modifier = modifier;
  } else {
    button.state    = Button::State::releasing;
    button.modifier = Button::Modifier::none;
  }
}

inline void update_button_state_for_next_frame(Button& button) {
  if (button.state == Button::State::pressing) {
    button.state = Button::State::holding;
  } else if (button.state == Button::State::releasing) {
    button.state = Button::State::nothing;
  }
}

void init_window(Window& win, const vec2i& size, const string& title,
    bool visible, int msaa) {
  // init glfw
  if (!glfwInit())
    throw std::runtime_error(
        "Cannot initialize GLFW context. Make sure the OpenGL context is intialized.");
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  if (!visible) {
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
  }

  // Set multisample anti-aliasing
  if (msaa > 1) {
    glfwWindowHint(GLFW_SAMPLES, msaa);
  }

  // create window
  win.glfw = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
  if (!win.glfw) throw std::runtime_error("Cannot initialize GLFW window.");
  glfwMakeContextCurrent(win.glfw);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowUserPointer(win.glfw, &win);

  init_callbacks(win);

  // TODO
  // glfwSetJoystickCallback([](int id, int event) { update_joysticks(win);
  // });

  //  init_joysticks(win);

  glfwSetWindowSizeCallback(
      win.glfw, [](GLFWwindow* glfw, int width, int height) {
        auto win = (Window*)glfwGetWindowUserPointer(glfw);
        glfwGetWindowSize(win->glfw, &win->size.x, &win->size.y);
        glfwGetFramebufferSize(win->glfw, &win->input.framebuffer_viewport.z,
            &win->input.framebuffer_viewport.w);
        win->input.framebuffer_viewport.x = 0;
        win->input.framebuffer_viewport.y = 0;
      });

  update_window_size(win);
}

void init_callbacks(Window& win) {
  // set callbacks
  glfwSetKeyCallback(win.glfw,
      [](GLFWwindow* glfw, int key, int scancode, int action, int mods) {
        auto win   = (Window*)glfwGetWindowUserPointer(glfw);
        auto press = (action == GLFW_PRESS);
        update_button_from_input(
            win->input.key_buttons[(int)key], press, (Button::Modifier)mods);
      });

  glfwSetMouseButtonCallback(
      win.glfw, [](GLFWwindow* glfw, int button, int action, int mods) {
        auto win      = (Window*)glfwGetWindowUserPointer(glfw);
        auto modifier = (Button::Modifier)mods;
        auto press    = (action == GLFW_PRESS);
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
          update_button_from_input(win->input.click_left, press, modifier);
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
          update_button_from_input(win->input.click_right, press, modifier);
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
          update_button_from_input(win->input.click_middle, press, modifier);
        }
      });

  glfwSetScrollCallback(
      win.glfw, [](GLFWwindow* glfw, double xoffset, double yoffset) {
        auto win          = (Window*)glfwGetWindowUserPointer(glfw);
        win->input.scroll = {(float)xoffset, (float)yoffset};
      });

  // TODO
  // if (win.callbacks.drop) {
  //   glfwSetDropCallback(
  //       win.glfw, [](GLFWwindow* glfw, int num, const char** paths) {
  //         auto win   = (Window*)glfwGetWindowUserPointer(glfw);
  //         auto pathv = vector<string>();
  //         for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
  //         win->drop(pathv);
  //       });
  // }

  // if (win.callbacks.focus) {
  //   glfwSetWindowFocusCallback(win.glfw, [](GLFWwindow* glfw, int focus) {
  //     auto win = (Window*)glfwGetWindowUserPointer(glfw);
  //     win->focus(focus);
  //   });
  // }
}

void delete_window(Window& win) {
  glfwDestroyWindow(win.glfw);
  glfwTerminate();
  win.glfw = nullptr;
}

bool should_window_close(const Window& win) {
  return glfwWindowShouldClose(win.glfw);
}

void poll_events(const Window& win, bool wait) {
  if (wait)
    glfwWaitEvents();
  else
    glfwPollEvents();
}

void run_draw_loop(Window& win, function<void(const Input&)> draw, bool wait) {
  while (!should_window_close(win)) {
    update_window_size(win);
    // reset input
    {
      update_button_state_for_next_frame(win.input.click_left);
      update_button_state_for_next_frame(win.input.click_right);
      for (auto& key : win.input.key_buttons) {
        update_button_state_for_next_frame(key);
      }
      // TODO
      //      for (auto& [key, value] : win.input.gui_buttons) {
      //        value = false;
      //      }
      win.input.scroll = zero2f;
    }
    poll_events(win, wait);
    update_input(win.input, win);
    //    update_joystick_input(win);
    draw(win.input);

    swap_buffers(win);
  }
}

// TODO(giacomo): remove
vec2f get_mouse_pos_normalized(const Window& win, bool isometric) {
  auto& pos    = win.input.mouse_pos;
  auto  size   = vec2f{(float)win.size.x, (float)win.size.y};
  auto  result = vec2f{2 * (pos.x / size.x) - 1, 1 - 2 * (pos.y / size.y)};
  if (isometric) {
    result.x *= float(win.size.x) / float(win.size.y);
  }
  return result;
}
vec2f get_mouse_pos_normalized(const Input& input, bool isometric) {
  auto& pos    = input.mouse_pos;
  auto  size   = vec2f{(float)input.window_size.x, (float)input.window_size.y};
  auto  result = vec2f{2 * (pos.x / size.x) - 1, 1 - 2 * (pos.y / size.y)};
  if (isometric) {
    result.x *= float(input.window_size.x) / float(input.window_size.y);
  }
  return result;
}

bool is_pressing(Button button) {
  return button.state == Button::State::pressing;
}
bool is_down(Button button) {
  return (button.state == Button::State::pressing ||
          button.state == Button::State::holding);
}
bool is_releasing(Button button) {
  return button.state == Button::State::releasing;
}
bool is_pressing(const Input& input, Key key) {
  return is_pressing(input.key_buttons[(int)key]);
}
bool is_down(const Input& input, Key key) {
  return is_down(input.key_buttons[(int)key]);
}
bool is_releasing(const Input& input, Key key) {
  return is_releasing(input.key_buttons[(int)key]);
}

void update_window_size(Window& win) {
  auto& glfw = win.glfw;
  glfwGetWindowSize(glfw, &win.size.x, &win.size.y);
  glfwGetFramebufferSize(glfw, &win.input.framebuffer_viewport.z,
      &win.input.framebuffer_viewport.w);
  win.input.framebuffer_viewport.x = 0;
  win.input.framebuffer_viewport.y = 0;
  win.input.framebuffer_size       = {
      win.input.framebuffer_viewport.z - win.input.framebuffer_viewport.x,
      win.input.framebuffer_viewport.w - win.input.framebuffer_viewport.y};

  // TODO(giacomo): cleanup
  win.input.window_size = win.size;
}

void update_input(Input& input, const Window& win) {
  // update input
  input.mouse_last = input.mouse_pos;
  auto  mouse_posx = 0.0, mouse_posy = 0.0;
  auto& glfw = win.glfw;
  glfwGetCursorPos(glfw, &mouse_posx, &mouse_posy);
  input.mouse_pos    = vec2f{(float)mouse_posx, (float)mouse_posy};
  input.modifier_alt = glfwGetKey(glfw, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
                       glfwGetKey(glfw, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
  input.modifier_shift = glfwGetKey(glfw, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                         glfwGetKey(glfw, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
  input.modifier_ctrl = glfwGetKey(glfw, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
                        glfwGetKey(glfw, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;

  // TODO(giacomo): separate gui
  // if (win.gui_width) {
  //   auto io             = &ImGui::GetIO();
  //   input.is_gui_active = io->WantTextInput || io->WantCaptureMouse ||
  //                         io->WantCaptureKeyboard;
  // }

  input.is_window_focused = glfwGetWindowAttrib(win.glfw, GLFW_FOCUSED) != 0;

  // time
  input.clock_last = input.clock_now;
  input.clock_now =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  input.time_now   = (double)input.clock_now / 1e9;
  input.time_delta = (double)(input.clock_now - input.clock_last) / 1e9;
  input.frame += 1;
}

void update_input(Window& win) { update_input(win.input, win); }

void update_joystick_input(vector<Joystick>& joysticks) {
  for (auto& joystick : joysticks) {
    GLFWgamepadstate state;
    if (glfwGetGamepadState(joystick.id, &state)) {
      int  count;
      auto axes              = glfwGetJoystickAxes(joystick.id, &count);
      joystick.left_stick    = {axes[0], -axes[1]};
      joystick.right_stick   = {axes[2], -axes[5]};
      joystick.left_trigger  = axes[4];
      joystick.right_trigger = axes[3];

      auto buttons = glfwGetJoystickButtons(joystick.id, &count);
      for (int i = 0; i < count; i++) {
        joystick.buttons[i] = buttons[i];
      }
    }
  }
}

void swap_buffers(const Window& win) { glfwSwapBuffers(win.glfw); }

}  // namespace window
