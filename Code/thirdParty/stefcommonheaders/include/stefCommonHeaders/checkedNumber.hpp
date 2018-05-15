#pragma once
#include <type_traits>
#include <exception>
#include <cassert>
#include <string>
#include <functional>
#include <iostream>

namespace stefCommonHeaders {


    class InvariantViolatedException: public std::exception {
    public:
        template <typename Number>
        InvariantViolatedException(
            const std::string& invariantDescription,
            const Number& value
        ): exceptionMessage_M("Invariant \"" + invariantDescription + "\" violated for value " + std::to_string(value)) {
        }

        const char*
        what(
        ) const noexcept {
            return exceptionMessage_M.c_str();
        }

    private:
        const std::string exceptionMessage_M;
    };


    #ifndef NO_CHECKED_NUMBERS
    template <typename Number>
    class CheckedNumber {
    public:
        template <typename T>
        using Invariant = std::function<bool(T)>;



        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value, OtherNumber>::type>
        CheckedNumber(
             const CheckedNumber<OtherNumber> & other
        ): value_M(other.value_M), invariant_M(other.invariant_M), invariantDescription_M(other.invariantDescription_M) {
        }



        operator Number(
        ) const {
            return value_M;
        }

        CheckedNumber<Number>& operator =(const Number& other) {
            value_M = other;
            checkInvariant();
            return *this;
        }


        CheckedNumber<Number>& operator +=(const Number& other) {
            value_M += other;
            checkInvariant();
            return *this;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value, void>::type>
        Number operator +(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M + otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value, void>::type>
        Number operator -(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M - otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value, void>::type>
        Number operator /(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M / otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value, void>::type>
        Number operator *(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M * otherNumber.value_M;
        }

        CheckedNumber<Number>& operator -=(const Number& other) {
            value_M -= other;
            checkInvariant();
            return *this;

        }

        CheckedNumber<Number>& operator *=(const Number& other) {
            value_M *= other;
            checkInvariant();
            return *this;

        }

        CheckedNumber<Number>& operator /=(const Number& other) {
            value_M /= other;
            checkInvariant();
            return *this;
        }

        CheckedNumber<Number>& operator++() {
            return operator+=(1);
        }

        CheckedNumber<Number> operator++(int) {
            const Number oldValue = value_M;
            operator++();
            return *this;
        }

        CheckedNumber<Number>& operator--() {
            return operator-=(1);
        }

        CheckedNumber<Number> operator--(int) {
            const Number oldValue = value_M;
            operator--();
            return *this;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value && std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        Number operator &(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M & otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value && std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        Number operator |(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M | otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value && std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        Number operator %(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M % otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value && std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        Number operator <<(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M << otherNumber.value_M;
        }

        template <typename OtherNumber, typename = typename std::enable_if<std::is_convertible<OtherNumber, Number>::value && std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        Number operator >>(const CheckedNumber<OtherNumber>& otherNumber) {
            return value_M >> otherNumber.value_M;
        }

        template <typename OtherNumber = Number, typename = typename std::enable_if<std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        CheckedNumber<OtherNumber> operator %=(const OtherNumber& other) {
            value_M %= other;
            checkInvariant();
            return *this;
        }

        template <typename OtherNumber = Number, typename = typename std::enable_if<std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        CheckedNumber<OtherNumber> operator |=(const OtherNumber& other) {
            value_M |= other;
            checkInvariant();
            return *this;
        }

        template <typename OtherNumber = Number, typename = typename std::enable_if<std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        CheckedNumber<OtherNumber> operator &=(const OtherNumber& other) {
            value_M &= other;
            checkInvariant();
            return *this;
        }

        template <typename OtherNumber = Number, typename = typename std::enable_if<std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        CheckedNumber<OtherNumber> operator <<=(const OtherNumber& other) {
            value_M <<= other;
            checkInvariant();
            return *this;
        }

        template <typename OtherNumber = Number, typename = typename std::enable_if<std::is_integral<OtherNumber>::value && std::is_integral<Number>::value, void>::type>
        CheckedNumber<OtherNumber> operator >>=(const OtherNumber& other) {
            value_M >>= other;
            checkInvariant();
            return *this;
        }
    
    protected:
        CheckedNumber(
            Number value,
            const Invariant<Number>& invariant,
            const std::string& invariantDescription
        ): value_M(value), invariant_M(invariant), invariantDescription_M(invariantDescription) {
            checkInvariant();
        }
    private:
        const Invariant<Number> invariant_M;
        Number value_M;
        const std::string invariantDescription_M;


        void
        checkInvariant() {
            if (!invariant_M(value_M)) {
                throw InvariantViolatedException(invariantDescription_M, value_M);
            }
        }


    };

    template <typename Number>
    std::ostream& operator<<(std::ostream& os, const CheckedNumber<Number>& number) {
        return os << Number(number);
    }

    template <typename Number>
    std::istream& operator>>(std::istream& is, CheckedNumber<Number>& number) {
        Number n;
        is >> n;
        if (is.good()) {
            number = n;
        }
        return is;
    }

    template <typename Number>
    class PositiveNumber : public CheckedNumber<Number> {
    public:
        PositiveNumber(const Number& number):
            CheckedNumber<Number>(number, positive, positiveDescription_C) {
            }
    private:
        static bool positive(const Number& number) {
            return number > 0;
        }

        const static std::string positiveDescription_C;
    };

    template <typename Number>
    const std::string PositiveNumber<Number>::positiveDescription_C = "> 0";

    template <typename Number>
    class NonNaNNumber: public CheckedNumber<Number> {
    public:
        NonNaNNumber(const Number& number):
              CheckedNumber<Number>(number, nonnan, nonnanDescription_C) {
        }

        NonNaNNumber():
              NonNaNNumber(0) {}

        NonNaNNumber(const NonNaNNumber<Number>& other):
              NonNaNNumber<Number>(Number(other)) {
        }


    private:
        static bool nonnan(const Number& number) {
            return ! std::isnan(number);
        }

        const static std::string nonnanDescription_C;
    };

    template <typename Number>
    const std::string NonNaNNumber<Number>::nonnanDescription_C = "!= NaN";


    #else
    template <typename T>
    using CheckedNumber = T;

    template <typename T>
    using PositiveNumber = T;
    #endif


}
