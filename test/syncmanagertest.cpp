#include <gtest/gtest.h>

#include "asyncresult.h"
#include "hybridmpisyncmanager.h"

TEST(HybridMpiSyncManager, StreamTest)
{
    using namespace HybridMpi;

    auto result = std::make_shared<AsyncResult>();
    result->betaOld = 1;
    result->beta = 2;
    result->v = std::make_unique<MatrixXd>(MatrixXd::Zero(1, 3));
    result->component = 4;

    std::vector<ConstIndexResultPair> out;
    out.emplace_back(0, result);
    out.emplace_back(1, result);

    const auto size = sendSize(out);

    std::unique_ptr<char[]> buffer;
    {
        buffer.reset(new char[static_cast<unsigned long>(size)]);

        std::ostringstream stream;
        stream.rdbuf()->pubsetbuf(buffer.get(), size);

        writeResultPairs(out, &stream);
    }

    std::vector<IndexResultPair> in;

    {
        std::istringstream stream;
        stream.rdbuf()->pubsetbuf(buffer.get(), size);
        in = readResultPairs(&stream);
    }

    ASSERT_EQ(out.size(), in.size());
    for (size_t i = 0; i < out.size(); ++i) {
        const bool equal = *(out.at(i).second) == *(in.at(i).second);
        ASSERT_TRUE(equal);
    }
}
