#pragma once
#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <glog/logging.h>
#include <maki/maki.h>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;

// constexpr string literal to pass string literals as template arguments,
// from Kevin Hartman - https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
template <size_t N>
struct StringLiteral
{
    constexpr StringLiteral(const char (&msg)[N])
    {
        std::copy_n(msg, N, value);
    }

    char value[N];
};

template <typename... Args>
constexpr std::string print_variadic_types(Args... args)
{
    std::string typestr;
    ([&]
     {
        if (typestr.size() > 0) {
            typestr.append(", ");
        }
        typestr.append(typeid(std::forward<Args>(args)).name()); }(), ...);
    return typestr;
}

// void impl from https://medium.com/@barryrevzin/without-form-and-void-cfc62091d4d6
struct Void
{
    Void() = default;
    Void(Void const &) = default;
    Void(Void &&) = default;
    Void &operator=(Void const &) = default;
    Void &operator=(Void &&) = default;

    template <typename Arg, typename... Args,
              std::enable_if_t<!std::is_base_of_v<Void, std::decay_t<Arg>>, int> = 0>
    explicit Void(Arg &&, Args &&...) {}
};

template <typename T>
using wrap_void_t = std::conditional_t<std::is_void_v<T>, Void, T>;

template <typename T>
using unwrap_void_t = std::conditional_t<std::is_same_v<std::decay_t<T>, Void>, void, T>;

// callable returning non-void

template <
    class Function, class... Args,
    typename Result = std::invoke_result_t<Function, Args...>,
    std::enable_if_t<!std::is_void_v<Result>, int> = 0>
constexpr Result invoke_void(Function &&f, Args &&...args)
{
    return std::invoke(std::forward<Function>(f), std::forward<Args>(args)...);
}

// callable returning void

template <
    class Function, class... Args,
    typename Result = std::invoke_result_t<Function, Args...>,
    std::enable_if_t<std::is_void_v<Result>, int> = 0>
constexpr Void invoke_void(Function &&f, Args &&...args)
{
    std::invoke(std::forward<Function>(f), std::forward<Args>(args)...);
    return Void();
}

// timing contexts --------------------------------------------------------------

#if DO_TIMING
// timing region of code that only acts on non-local references
#define timing_context(msg, ...)                           \
    {                                                      \
        auto t1 = high_resolution_clock::now();            \
        __VA_ARGS__                                        \
        auto t2 = high_resolution_clock::now();            \
        auto ms_int = duration_cast<nanoseconds>(t2 - t1); \
        LOG(INFO) << "@timing "                            \
                  << msg                                   \
                  << " " << ms_int.count() << " ns";       \
    }
#else
// timing is not being performed, set DO_TIMING in CMakeLists.txt to enable.
#define timing_context(msg, ...) __VA_ARGS__
#endif

// timing functions -------------------------------------------------------------

#if DO_TIMING // auto result = invoke_void(std::forward<decltype(&__VA_ARGS__)>(__VA_ARGS__), std::forward<decltype(args)>(args)...);
#define time_this_function(...) \
    ([](auto &&...args) -> decltype(auto) { \
            auto t1 = high_resolution_clock::now(); \
            auto result = __VA_ARGS__(std::forward<decltype(args)>(args)...); \
            auto t2 = high_resolution_clock::now(); \
            auto ms_int = duration_cast<nanoseconds>(t2 - t1); \
            LOG(INFO) << "@timing " \
                    << #__VA_ARGS__ \
                    << " " << ms_int.count() << " ns"; \
            return result; })
#else
#define time_this_function(...) __VA_ARGS__
#endif

#if DO_TIMING
#define time_this_member(obj, mem) time_this_function(std::bind(&std::remove_cvref_t<decltype(obj)>::mem, obj))
#else
#define time_this_member(obj, mem) obj.mem
#endif

// timed operator() member ------------------------------------------------------------

template <class Op, StringLiteral msg>
struct operator_timing_wrapper : public Op
{
    // overload constructors
    template <class... Args>
    operator_timing_wrapper(Args &&...args) : Op(std::forward<Args>(args)...) {}

    // return types for given Arg set
    template <typename... Args>
    using ReturnType = decltype(std::declval<Op>()(std::declval<Args>()...));

    template <class... Args>
    using ConstOperatorType = ReturnType<Args...> (Op::*)(Args...) const;

    // overload operator()
    template <class... Args>
    auto operator()(Args &&...args) const
    {
        auto t1 = high_resolution_clock::now();
        auto result = invoke_void(static_cast<ConstOperatorType<Args...>>(&Op::operator()), *this, std::forward<Args>(args)...);
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<nanoseconds>(t2 - t1);
        LOG(INFO) << "@timing "
                  << msg.value
                  //<< "(" << print_variadic_types(std::forward<Args>(args)...) << ") "
                  << " " << ms_int.count() << " ns";
        return result;
    }
};

#if DO_TIMING
#define time_this_operator(...) operator_timing_wrapper<__VA_ARGS__, #__VA_ARGS__>
#else
#define time_this_operator(...) __VA_ARGS__
#endif