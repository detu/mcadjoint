#pragma once
#include <stdexcept>
namespace stefCommonHeaders {

    #ifdef NDEBUG
    template <typename T>
    using YetToInitialize = T;
    #else


    template <typename T>
    class YetToInitialize {
    public:
        YetToInitialize(
        ): haveBeenInitialized(false) {
        }

        T&
        operator=(const T& other) {
            haveBeenInitialized = true;
            value = other;
            return value;
        }

        // From https://stackoverflow.com/questions/18100297/how-can-i-use-stdenable-if-in-a-conversion-operator
        // Enable the conversion operator only for fundamental types or pointers or references or arrays to them

        using BareT = typename std::remove_cv<
            typename std::remove_reference<
                typename std::remove_pointer<
                    typename std::remove_all_extents<T>::type
                >::type
            >::type
        >::type;

        using TIfItsFundamentalOtherwiseVoid = typename std::conditional<std::is_fundamental<BareT>::value, T, void>::type;

        operator TIfItsFundamentalOtherwiseVoid() const {
            failIfUninitialized();
            return (TIfItsFundamentalOtherwiseVoid) value;
        }

    private:

        void
        failIfUninitialized() const {
            if (!haveBeenInitialized) {
                throw std::logic_error("This variable has been accessed without being initialized!");
            }
        }

        T value;
        bool haveBeenInitialized;
    };
    #endif
}
