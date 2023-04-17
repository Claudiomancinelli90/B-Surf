# Serialize
Header-only library for binary serialization and deserialization of simple data structures, in order to quickly load and/or write of the state of a program.  
`serialize.h` provides the basic tools to easily define the serialization code of custom data structures.  
In order to minimize disk access, a memory buffer is used during serialization.  
The library features the following built-in serializtion functions:
  * `serialize()`: write/read any POD (struct with no allocated resources).
  * `serialize_string()`: write/read a `std::string`.
  * `serialize_vector()`: write/read a `std::vector` of PODs.
  * `serialize_vector()`: write/read a `std::vector` of structs with custom serialization function.

# Example
In `example/test.cpp` the library is tested, showing its usage on a toy data structure.  
Build it with `build.sh`.

```C++
#include "serialize.h"

struct Object {
    std::string name;
    std::vector<int> vec;
    int i;
    float val;
};

// Define custom serialization function.
void serialize_object(Serializer& srl, Object& var) {
    serialize(srl, var.i);
    serialize(srl, var.val);
    serialize_string(srl, var.name);
    serialize_vector(srl, var.vec);
}

int main() {
    auto object   = Object{"Hello", {1,2,3,4}, 77, 12.0};
    auto filename = "object.bin"s;
    int  capacity = 64;
    
    // Save object into binary a binary file.
    auto writer = make_writer(filename, capacity);
    serialize_object(writer, object);
    close_serializer(writer);

    // Reload object from disk.
    auto object_reloaded = Object{};
    auto reader = make_reader(filename, capacity);
    serialize_object(reader, object_reloaded);
    close_serializer(reader);

    // object == object_reloaded
}
```


