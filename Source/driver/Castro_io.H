#ifndef CASTRO_IO_H
#define CASTRO_IO_H

#include <iostream>
#include <sstream>

extern std::string inputs_name;

// Redirect std::cout to a temporary buffer.

class CoutRedirection {

public:
    CoutRedirection () : original_output(std::cout.rdbuf()) {
        std::cout.rdbuf(captured_output.rdbuf());
    }

    ~CoutRedirection () {
        std::cout.rdbuf(original_output);
    }

    std::string getCapturedOutput () {
        return captured_output.str();
    }

private:
    std::ostringstream captured_output;
    std::streambuf* original_output;

};

#endif
