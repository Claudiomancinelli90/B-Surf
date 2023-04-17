# Graphics
Small library for computer graphics. This is simply a subset of [yocto-gl](https://github.com/xelatihy/yocto-gl) with some modifications. Here are the main differences:
- Gpu commands and window management are split into two indipendent libraries.
- No RAII. Manual cleanup of gpu resources.
- Unified concept of gpu buffer into `Arraybuffer`.
- You are given low level access to shapes (i.e. vertex-array-objects) and draw calls. There are no high-level concepts like "scene", "material", "instance".
- Render buffers, useful for deferred shading and multipass techniques.
- Joystick input.