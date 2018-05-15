#ifndef STEFCOMMONHEADERS_SIMPLETIMING_HPP
#define STEFCOMMONHEADERS_SIMPLETIMING_HPP
#pragma once
#include <type_traits> // for std::result_of
#include <utility> // for std::move
#include <cassert> // for assert
#include <vector> // for the median call
#include <algorithm> // for nth_element
#include <chrono> // for the timing functions

#if __cplusplus < 201103L
#error "This code needs at least C++11!"
#endif

namespace stefCommonHeaders {


/*!
 * \file
 * \author Stefano Weidmann
 * \version 1
 * \date 2016
 * \brief Functions timing a call. Threadsafe.
 * \details Takes away the hassle of using std::chrono if all you want is to measure the runtime (wall clock) in seconds of function calls.
 * \pre requires C++ 11
 */

/**************************************************************************************************
 *                                         INTERNAL STUFF                                         *
 **************************************************************************************************/

#ifdef __GNUC__
#define TIMECALL_INTERNAL_HIDDEN_NAMESPACE __attribute__((visibility("internal")))
#else
#define TIMECALL_INTERNAL_HIDDEN_NAMESPACE
#endif


namespace timeCallInternal TIMECALL_INTERNAL_HIDDEN_NAMESPACE {

    // make sure to use a steady clock (time monotonically increases even if thread sleeps)
    // to use it in multithreaded settings
    using SteadyClock = typename std::conditional<
        std::chrono::high_resolution_clock::is_steady,
        std::chrono::high_resolution_clock,
        std::chrono::steady_clock>::type;

    /*!
     * \internal
     * \details Structure abstracting away if the function in timeCall returns a value or not. Version returning a value.
     */
    template <typename T>
    struct ReturnValueCache {

        template <typename Callable, class... Args>
        inline
        ReturnValueCache(
            Callable toCall,
            Args... argsForTheCall):
                toReturn_M(toCall(argsForTheCall...)) {
        }
        

        // Nicer way to construct ReturnValueCache (could specify template arguments)
        template <typename Callable, class... Args>
        static inline
        ReturnValueCache<T>
        make(
            Callable toCall,
            Args... argsForTheCall
        ) {
            return ReturnValueCache<T>(toCall, argsForTheCall...);
        }

        inline
        T
        doReturn() const {
            return std::move(toReturn_M);
        }

        T toReturn_M;
    };

    /*!
     * \internal
     * \details Structure abstracting away if the function in timeCall returns a value or not. Version not returning anything.
     */
    template <>
    struct ReturnValueCache<void> {

        template <typename Callable, class... Args>
        inline
        ReturnValueCache(
            Callable toCall,
            Args... argsForTheCall
        ) {

            toCall(argsForTheCall...);
        }

        // Nicer way to construct ReturnValueCache (could specify template arguments)
        template <typename Callable, class... Args>
        static inline
        ReturnValueCache<void>
        make(
            Callable toCall,
            Args... argsForTheCall
        ) {
            return ReturnValueCache<void>(toCall, argsForTheCall...);
        }

        inline
        void
        doReturn() const {
            return;
        }
    };

}

/**************************************************************************************************
 *                            THE ACTUAL TIMECALL FAMILY OF FUNCTIONS                             *
 **************************************************************************************************/



    /*!
     * \brief Times a call in seconds.
     * \details Calls a function or function-like object with the specified object and writes the elapsed wall time into a variable.
     * \param measuredTimeInSeconds Where to write the elapsed wall time (in seconds) into.
     * \param toCall The function or function-like object to call.
     * \param argsForTheCall The arguments for the call.
     * \return Whatever a simple call to the function ("toCall(argsToCall...)") would return.
     * \exception Any exception the function itself raises is passed through
     * \pre toCall(argsForTheCall...) has to be a valid call
     * \post executes toCall(argsForTheCall...) including all side-effects. Writes the elapsed time in seconds into measuredTimeInSeconds.
     */
    template <typename Callable, class... Args>
    typename std::result_of<Callable(Args...)>::type
    timeCall(
        double& measuredTimeInSeconds,
        Callable toCall,
        Args... argsForTheCall
    ) {
        using namespace timeCallInternal;
        using namespace std::chrono;
        using TimePoint = SteadyClock::time_point;

        const TimePoint timeBeforeCall(SteadyClock::now());
        const auto returnValueOfCall(ReturnValueCache<typename std::result_of<Callable(Args...)>::type>::make(toCall, argsForTheCall...));
        const TimePoint timeAfterCall(SteadyClock::now());

        const duration<double> callDuration(timeAfterCall - timeBeforeCall);

        measuredTimeInSeconds = callDuration.count();

        return returnValueOfCall.doReturn();
    }

    /*!
     * \brief Times a call multiple times and write the median elapsed wall time. It returns the returnvalue of any one of those calls.
     * \details Calls a function or function-like object with the specified object and writes the median elapsed wall time into a variable.
     *          Returns the return value of one of the executed calls. It's undefined which call.
     * \param medianTimeInSeconds Where to write the median elapsed wall time (in seconds) into.
     * \param numberOfCalls How many times to call
     * \param toCall The function or function-like object to call.
     * \param argsForTheCall The arguments for the call.
     * \return Whatever a simple call to the function ("toCall(argsToCall...)") would return. It is not specified of what call you get the return value
     * \exception Any exception the function itself raises is passed through
     * \pre toCall(argsForTheCall...) has to be a valid call
     * \post executes toCall(argsForTheCall...) a number of times including all side-effects. Writes the median elapsed time in seconds into medianTimeInSeconds.
     */

    template <typename Callable, class... Args>
    typename std::result_of<Callable(Args...)>::type
    medianTimeCall(
        double& medianTimeInSeconds,
        const int numberOfCalls,
        Callable toCall,
        Args... argsForTheCall
    ) {

        using namespace timeCallInternal;
        assert(numberOfCalls >= 1);
        std::vector<double> runTimes;
        runTimes.reserve(numberOfCalls);

        double currentlyMeasuredTimeInSeconds;
        for (int callNumber = 0; callNumber < numberOfCalls - 1; ++callNumber) {
            (void) timeCall(currentlyMeasuredTimeInSeconds, toCall, argsForTheCall...);
            runTimes.push_back(currentlyMeasuredTimeInSeconds);
        }

        using ReturnType = typename std::result_of<Callable(Args...)>::type;


        const auto returnValueOfLastCall =
            ReturnValueCache<ReturnType>::make(timeCall<Callable, Args...>, currentlyMeasuredTimeInSeconds, toCall, argsForTheCall...);
        runTimes.push_back(currentlyMeasuredTimeInSeconds);

        assert(runTimes.size() == numberOfCalls);


        const auto usualMedianIterator = runTimes.begin() + (numberOfCalls + 1)/2;
        std::nth_element(runTimes.begin(), usualMedianIterator, runTimes.end());

        medianTimeInSeconds = *usualMedianIterator;
        if (numberOfCalls % 2 == 0) {
            const auto otherMedianIterator = usualMedianIterator + 1;
            std::nth_element(runTimes.begin(), otherMedianIterator, runTimes.end());
            medianTimeInSeconds += *otherMedianIterator;
            medianTimeInSeconds /= 2;
        }

        return returnValueOfLastCall.doReturn();
    }


}

#undef TIMECALL_INTERNAL_HIDDEN_NAMESPACE
#ifdef __GNUC__
#pragma GCC poison timeCallInternal
#endif
#endif