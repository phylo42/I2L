#ifndef I2L_SERIALIZATION_H
#define I2L_SERIALIZATION_H

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <utility>
#include <fstream>
#include <memory>
#include "phylo_kmer_db.h"
#include "phylo_node.h"
#include "version.h"

namespace fs = boost::filesystem;

namespace i2l
{
    /// A struct with definitions for different version of the serialization protocol
    struct protocol
    {
        /// IPK v0.1.x
        static const unsigned int v0_1_x = 2;

        /// IPK v0.2.x + v0.3.1
        static const unsigned int v0_2_WITHOUT_POSITIONS = 3;
        static const unsigned int v0_2_WITH_POSITIONS = v0_2_WITHOUT_POSITIONS + 1;

        /// IPK v0.4.0, added tree indexing
        static const unsigned int v0_4_0_WITHOUT_POSITIONS = 5;
        static const unsigned int v0_4_0_WITH_POSITIONS = v0_4_0_WITHOUT_POSITIONS + 1;

        /// IPK v0.4.1, added kmer ordering and filter values
        static const unsigned int v0_4_1_WITHOUT_POSITIONS = 7;
        static const unsigned int v0_4_1_WITH_POSITIONS = v0_4_1_WITHOUT_POSITIONS + 1;

#ifdef KEEP_POSITIONS
        static const unsigned int EARLIEST_INDEX = v0_4_0_WITH_POSITIONS;
        static const unsigned int CURRENT = v0_4_1_WITH_POSITIONS;
#else
        static const unsigned int EARLIEST_INDEX = v0_4_0_WITHOUT_POSITIONS;
        static const unsigned int CURRENT = v0_4_1_WITHOUT_POSITIONS;
        static const unsigned int CURRENT_WITH_POSITIONS = CURRENT + 1;
#endif
        static const unsigned int ERROR = 42;
    };

    i2l::phylo_kmer_db load(const std::string& filename, float mu=1.0, float user_epsilon=0.0f,
                            size_t max_entries=std::numeric_limits<size_t>::max());
    i2l::phylo_kmer_db load_compressed(const std::string& filename, float user_mu,
                                       float user_omega, size_t max_entries);
    i2l::phylo_kmer_db load_uncompressed(const std::string& filename, float user_mu,
                                         float user_omega, size_t max_entries);

    void save(const i2l::phylo_kmer_db& db, const std::string& filename, bool uncompressed = false);
    void save_compressed(const i2l::phylo_kmer_db& db, const std::string& filename);
    void save_uncompressed(const i2l::phylo_kmer_db& db, const std::string& filename);

    using oarchive = boost::archive::binary_oarchive;

}


#endif
