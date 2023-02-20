#ifndef I2L_OPTIONAL_H
#define I2L_OPTIONAL_H

#if __has_include(<optional>)
#include <optional>
using std::optional;
using std::nullopt;
#elif __has_include(<experimental/optional>)
#include <experimental/optional>
using std::experimental::optional;
using std::experimental::nullopt;
#endif

#endif
