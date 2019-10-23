#ifndef TESTHELPERS_H
#define TESTHELPERS_H

#include <gtest/gtest.h>

class Options;

namespace testing {

    AssertionResult IsValidPreprocessOutPut(const Options  &options);
}
#endif // TESTHELPERS_H
