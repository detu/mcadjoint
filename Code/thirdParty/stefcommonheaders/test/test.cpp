#undef NDEBUG
#include <argh.h> // Not needed, just to make sure that it can be found
#include <stefCommonHeaders/dbg.hpp>
#include <stefCommonHeaders/intPow.hpp>
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/functionTraits.hpp>
#include <stefCommonHeaders/eigenToCsv.hpp>
#include <stefCommonHeaders/simpleTiming.hpp>
#include <stefCommonHeaders/stefFenv.h>
#include <stefCommonHeaders/yetToInitialize.hpp>
#include <stefCommonHeaders/checkedNumber.hpp>
#include <acutest.h>
#include <string>
#include <cstdio>
#include <cmath>
#include <setjmp.h>



using namespace stefCommonHeaders;


template <typename F>
double evaluateAtZero(F realFunction) {
    static_assert(FunctorIsCompatibleToSignature<F, double(double)>::value, "Hm");
    return realFunction(0);
}

bool fileExists(const std::string& fileName) {
    std::ifstream infile(fileName.c_str());
    return infile.good();
}


void testCheckedNumberShouldThrow() {
    try {
        PositiveNumber<double> d = 0;
    } catch (const InvariantViolatedException& e) {
        TEST_CHECK(true);
        return;
    }

    TEST_CHECK_(false, "Should throw an InvariantViolatedException but didn't!");
}

void testCheckedNumberOperator() {
    try {
        PositiveNumber<double> d = 1;
        --d;
    } catch (const InvariantViolatedException& e) {
        TEST_CHECK(true);
        return;
    }

    TEST_CHECK_(false, "Should throw an InvariantViolatedException but didn't!");
}

void testCheckedNumberIntOperator() {
    try {
        PositiveNumber<int> d = 1;
        d &= 0;
    } catch (const InvariantViolatedException& e) {
        TEST_CHECK(true);
        return;
    }

    TEST_CHECK_(false, "Should throw an InvariantViolatedException but didn't!");
}

void testAssert() {
    try {
        ASSERT(false);
        TEST_CHECK_(false, "Assert should fail!");
    } catch (std::logic_error err) {
        TEST_CHECK(true);
        return;
    }
}



void testCheckedNumberShouldNotThrow() {
    try {
        PositiveNumber<double> d = 1;
    } catch (const InvariantViolatedException& e) {
        TEST_CHECK_(false, "Should not throw an InvariantViolatedException but did!");
        return;
    }

    TEST_CHECK(true);
}



void testIntPow() {
    TEST_CHECK(intPow<1>(2) == 2);
    TEST_CHECK(intPow<-1>(2.0) == 0.5);
    TEST_CHECK(intPow<0>(2) == 1);
    TEST_CHECK(intPow<3>(2) == 8);
}





void testWriteEigenToCsv() {
    const std::string csvFileName = "identityMatrix.csv";
    stefCommonHeaders::writeToCsv(csvFileName, Eigen::Matrix2d::Identity());
    TEST_CHECK_(fileExists(csvFileName), "File \"%s\" doesn't exist but it should!", csvFileName.c_str());
    std::remove(csvFileName.c_str());
}


void testSimpleTiming() {
    const auto justOne = [] () { return 1;};
    double timeInSeconds = -1;
    const int result = stefCommonHeaders::timeCall(timeInSeconds, justOne);

    TEST_CHECK_(timeInSeconds >= 0, "timeCall didn't set time");
    TEST_CHECK_(result == justOne(), "timeCall altered the result");
}

void testMedianTiming() {
    const auto justOne = [] () { return 1;};
    double timeInSeconds = -1;
    const int result = stefCommonHeaders::medianTimeCall(timeInSeconds, 10, justOne);

    TEST_CHECK_(timeInSeconds >= 0, "medianTimeCall didn't set time");
    TEST_CHECK_(result == justOne(), "medianTimeCall altered the result");
}


int signalExceptionCode = -99999;
sigjmp_buf env;

void writeOutSignalExceptionCode(
    int signal,
    siginfo_t* signalInfos,
    void* somethingIdontCareAbout
) {
    signalExceptionCode = signalInfos->si_code;
    siglongjmp(env, 1);
}



void checkHandlingOfFloatingPointException(const int exceptionToRaise, const int signalExceptionCodeToExpect) {
    if (sigsetjmp(env, 1) == 0) {
        feraiseexcept(exceptionToRaise);
    } else {
        TEST_CHECK_(signalExceptionCode == signalExceptionCodeToExpect,
            "Should get \"%s\" but instead got \"%s\"",
            StefFenv_GetFloatingPointExceptionDescription(signalExceptionCodeToExpect),
            StefFenv_GetFloatingPointExceptionDescription(signalExceptionCode));
    }
}



void testFloatingPointExceptions() {

    feenableexcept(FE_ALL_EXCEPT);
    StefFenv_RegisterFloatingPointExceptionHandler(writeOutSignalExceptionCode);

    const static int exceptionsToRaise[] = {
        FE_DIVBYZERO,
        FE_INEXACT,
        FE_INVALID,
        FE_OVERFLOW,
        FE_UNDERFLOW
    };

    const static int signalExceptionCodesToExpect[] = {
        FPE_FLTDIV,
        FPE_FLTRES,
        FPE_FLTINV,
        FPE_FLTOVF,
        FPE_FLTUND
    };

    for (int i = 0; i < sizeof(signalExceptionCodesToExpect) / sizeof(int); ++i) {
        checkHandlingOfFloatingPointException(exceptionsToRaise[i], signalExceptionCodesToExpect[i]);
        feclearexcept(FE_ALL_EXCEPT);

    }

    StefFenv_ClearFloatingPointExceptionHandler();
    fedisableexcept(FE_ALL_EXCEPT);

}


void testYetToInitialize() {
    YetToInitialize<int> i;

    try {
        std::cout << i;
    } catch (std::logic_error& e) {
        TEST_CHECK(true);
        return;
    }
    TEST_CHECK_(false, "Should throw a logic_error exception but didn't!");
}

TEST_LIST = {
    {"writeEigenToCsv", testWriteEigenToCsv},
    {"simpleTiming", testSimpleTiming},
    {"medianTiming", testMedianTiming},
    {"fpe", testFloatingPointExceptions},
    {"yetToInitialize", testYetToInitialize},
    {"intPow", testIntPow},
    {"checkedNumberShouldThrow", testCheckedNumberShouldThrow},
    {"checkedNumberShouldNotThrow", testCheckedNumberShouldNotThrow},
    {"checkedNumberOperator", testCheckedNumberOperator},
    {"CheckedNumberIntOperator", testCheckedNumberIntOperator},
    {"ASSERT", testAssert},
    {nullptr}
};
