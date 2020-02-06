#ifndef TESTHELPERS_H
#define TESTHELPERS_H

#include <gtest/gtest.h>

class Options;

namespace testing {

    // fixed effect file, fixed effect number
    struct FixedEffectsParams {
        unsigned int fixedEffectNumber = 0;
        std::string fixedFile = "";
    };

    AssertionResult IsValidPreprocessOutPut(const Options  &options);
}
#endif // TESTHELPERS_H
