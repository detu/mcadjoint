//
// Created by Stefano Weidmann on 17.12.17.
//

#ifndef STEFCOMMONHEADERS_FUNCTIONTRAITS_HPP
#define STEFCOMMONHEADERS_FUNCTIONTRAITS_HPP
#pragma once

#include <tuple>
#include <type_traits>
#include <functional>

namespace stefCommonHeaders {

    template <class F>
    struct FunctionTraits;

// function pointer
    template <class R, class... Args>
    struct FunctionTraits<R(*)(Args...)> : public FunctionTraits<R(Args...)> {
    };

    template <class R, class... Args>
    struct FunctionTraits<R(Args...)> {
        using ReturnType = R;

        static constexpr unsigned arity = sizeof...(Args);

        template <unsigned N>
        struct Argument {
            static_assert(N < arity, "Error: invalid parameter index.");
            using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
        };

        template <unsigned N>
        using ArgumentType = typename Argument<N>::type;

    };


    template <typename CandidateSignature, typename SignatureToBeCompatibleTo>
    struct FunctorIsCompatibleToSignature {
        constexpr static bool value = std::is_convertible<CandidateSignature, std::function<SignatureToBeCompatibleTo>>::value;
    };
}



#endif //STEFCOMMONHEADERS_FUNCTIONTRAITS_HPP
