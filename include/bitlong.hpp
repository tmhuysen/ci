

#ifndef RDM_BITLONG_HPP
#define RDM_BITLONG_HPP

#include <boost/dynamic_bitset.hpp>


namespace bmqc {

/** Give the next long permutation in a lexicographical sense.
*
*      examples:
*          011 -> 101
*          101 -> 110
*
* Taken from (http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation)
* Credits: Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
*
* (using longs is up to 10x to 40x faster than the equivalent bitset algorithm without the bitset conversion.)
*/

size_t next_long_permutation(const size_t &v);


/** Generate the smallest (in lexicographical sense) long, given a total number of ones (n_ones).
*
*      example:
*           nones=2 -> ...011
*/
size_t smallest_long(size_t n_ones);


/** Generate all bitset permutations in a lexicographical sense, based on a total length and a given number of ones (nones).
*/
std::vector<size_t> all_long_permutations(size_t length, size_t nones);


/** Concatenate two sets of longs, producing every possible combination between both. set1 is major.
*
*      example:
*          concatenate_bitset_sets {000}, {011, 101, 110} -> {000011, 000101, 000110}
*/
std::vector<size_t>
concatenate_long_sets(const std::vector<size_t> &set1, const std::vector<size_t> &set2,
                      size_t length_set2);


/** Count the number of differences between two longs.
*
*      example:
*          count_differences(001, 100) -> 2
*          count_differences(001, 011) -> 1
*/
//__builtin_popcount((u ^ v))


/** Slice a given bitset from least significant to most significant bit (i.e. from right to left) This is NOT as they are printed!
*
*     example:
*          u: 0100011101
*          bitslice(u, 2, 5) is NOT equal to 0001
*          bitslice(u, 2, 5) => 0100[0111]01 -> 0111
*/
size_t longslice(const size_t &u, size_t i, size_t j);

//add phaseCheck



/**
*  Applies the annihilation operator (without sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool annihilation(size_t &u, size_t i);

/**
*  Applies the creation operator (without sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool creation(size_t &u, size_t i);


/**
* Evaluates the sign change between two (spin separated) base vectors
* (only for vectors changed by one creation and one annihilation operators)
*/

int phase_factor(const size_t &u, const size_t &v);

/**
* Evaluates the sign change for a creation or annihilation operator.
*/

int operator_phase(const size_t &u, size_t i);

/**
*  Applies the annihilation operator (with sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool annihilation_s(size_t &u, size_t i, int &sign);

/**
*  Applies the creation operator (with sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool creation_s(size_t &u, size_t i, int &sign);

}

#endif //RDM_BITLONG_HPP
