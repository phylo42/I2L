#ifndef I2L_HASH_MAP_H
#define I2L_HASH_MAP_H

#include <tsl/hopscotch_map.h>

namespace i2l
{
    template<typename... Args>
    using hash_map = tsl::hopscotch_map<Args...>;
}
#endif
