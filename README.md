# b/Surf: Interactive Bézier Splines on Surface Meshes

This repository contains the implementation of the algorithms described in "b/Surf: Interactive Bézier Splines on Surface Meshes", C. Mancinelli, G.Nazzaro, F. Pellacini and E. Puppo. The code consists of one GUI (splinegui) supporting all the algorithms for curve tracing described in the paper.

## Compilation
You can compile the project using the standard CMake routine

`mkdir build`<br/>
`cd build`<br/>
`cmake ../`<br/>
`make`<br/>

## Dependencies
Mac users only need to have C++ 17 or later installed. Linux users may have to install some additional packages, as

`libxrandr-dev`<br/>
`libxinerama-dev`<br/>
`libxcursor-dev`<br/>
`libxi-dev.`<br/>

Similarly, Windows users may need to install some OS-related packages, 
but there are no dependencies that prevent the use of the code on such operating system. 

## Run
Once the projects are built, you can run the app from within the project repository by issuing

`./bin/splinegui <path_to_a_model>`


## Using the GUI

                                 ---Basics---
Rotate: SHIFT + dragging.
Panning: ALT + dragging.                                 
                                 ---Control Points Selection/Editing---
Control points can be picked with left click and moved on the mesh by left-clicking on them + dragging

                                ---Curve Tracing---
Once the third control point P3 is picked, a cubic Bézier curve is automatically traced by putting the last control point P3 equal to P4. You can edit it in real time by dragging P3 as described above.

The GUI allows you to select the methods and the possible parameters described in the paper. Every time a parameter changes, the curve is updated.








                         






