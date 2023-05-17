#include "op/op.h"

#include <iostream>

void Op::not_implemented(const char* msg) {
    Error("Not implemented: " << msg);
}
