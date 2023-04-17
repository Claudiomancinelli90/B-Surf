mkdir -p bin
mkdir -p output

COMPILE="g++ -std=c++1z -O2"

echo $COMPILE -o bin/test test.cpp
$COMPILE -o bin/test test.cpp

