#pragma once
#include <type_traits>

// using enable_if because partial specialization of function templates isn't allowed


namespace stefCommonHeaders {
   // prototypes
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent < 0, Scalar>::type
   intPow(const Scalar x);

   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 0, Scalar>::type
   intPow(const Scalar x);

   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<(exponent > 2) && (exponent % 2 == 0), Scalar>::type
   intPow(const Scalar x);

   template <int exponent, typename Scalar>
   constexpr typename std::enable_if< (exponent > 1) && (exponent % 2 == 1), Scalar>::type
   intPow(const Scalar x);


   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 1, Scalar>::type
   intPow(const Scalar x);

   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 2, Scalar>::type
   intPow(const Scalar x);


   // negative exponents
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent < 0, Scalar>::type
   intPow(const Scalar x) {
      static_assert(std::is_floating_point<Scalar>::value, "Negative exponents are only allowed with floating point numbers!");
      return Scalar(1.0) / intPow<-exponent, Scalar>(x);
   }

   // base cases
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 0, Scalar>::type
   intPow(const Scalar x) {
      return Scalar(1.0);
   }

   // not strictly necessary but cuts the recursion a little shorter
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 1, Scalar>::type
   intPow(const Scalar x) {
      return x;
   }


   // not strictly necessary but cuts the recursion a little shorter
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<exponent == 2, Scalar>::type
   intPow(const Scalar x) {
      return x*x;
   }


   // recursive cases
   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<(exponent > 2) && (exponent % 2 == 0), Scalar>::type
   intPow(const Scalar x) {
      return intPow<2, Scalar>(intPow<exponent / 2, Scalar>(x));
   }

   template <int exponent, typename Scalar>
   constexpr typename std::enable_if<(exponent > 1) && (exponent % 2 == 1), Scalar>::type
   intPow(const Scalar x) {
      return x * intPow<2, Scalar>(intPow<exponent / 2, Scalar>(x));
   }
}
