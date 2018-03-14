#include "bitlong.hpp"
#include <boost/math/special_functions.hpp>



//#FIXME : ADD error throws

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

size_t bmqc::next_long_permutation(const size_t &v) {
    size_t v_ul = v;

    // t gets v's least significant 0 bits set to 1
    size_t t = v_ul | (v_ul - 1);
    // Next set to 1 the most significant bit to change and set to 0 the least significant ones,
    // and add the necessary 1 bits.
    size_t w = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(v_ul) + 1));

    return  w;
}

/** Generate the smallest (in lexicographical sense) long, given a total number of ones (n_ones).
 *
 *      example:
 *           nones=2 -> ...011
 */

size_t bmqc::smallest_long(size_t n_ones) {
    size_t v = 1;
    for (size_t i = 1; i < n_ones; i++) {
        v <<= 1;
        v += 1;
    }
    return v;
}

/** Generate all bitset permutations in a lexicographical sense, based on a total length and a given number of ones (nones).
*/

std::vector<size_t> bmqc::all_long_permutations(size_t size, size_t n_ones) {
    // First, create an empty set of all possible bitset permutations
    //      The number of different permutations of the 1's in the bitset is equal to
    //          L choose N
    //      with L being the size of the bitset, and N being the number of ones in the bitset
    auto nperm = static_cast<size_t>(boost::math::binomial_coefficient<double>(static_cast<unsigned>(size), static_cast<unsigned>(n_ones))); // total number of permutations
    std::vector<size_t> permutations (nperm);  // std::vector only accepts size_t

    // The first permutation is the smallest lexicographical one
    permutations[0] = smallest_long(n_ones);

    for (size_t i = 1; i < nperm; i++) {
        permutations[i] = next_long_permutation(permutations[i - 1]);
    }

    return permutations;
}

/** Concatenate two sets of longs, producing every possible combination between both. set1 is major.
*
*      example:
*          concatenate_bitset_sets {000}, {011, 101, 110} -> {000011, 000101, 000110}
*/


std::vector<size_t> bmqc::concatenate_long_sets(const std::vector<size_t> &set1, const std::vector<size_t> &set2,
                                                       size_t length_set2) {
    std::vector<size_t> longset (set1.size()*set2.size());
    size_t i = 0;
    for (const size_t& bitset1 : set1) {
        for (const size_t& bitset2: set2) {
            longset[i] = (bitset1 << length_set2) + bitset2;
            i++;
        }
    }
    return longset;
}

/** Slice a given bitset from least significant to most significant bit (i.e. from right to left) This is NOT as they are printed!
*
*     example:
*          u: 0100011101
*          bitslice(u, 2, 5) is NOT equal to 0001
*          bitslice(u, 2, 5) => 0100[0111]01 -> 0111
*/

size_t bmqc::longslice(const size_t &u, size_t i, size_t j) {
    size_t v = u >> i;
    v &= smallest_long(j-i+1);;
    return v;
}


/**
*  Applies the annihilation operator (without sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/

bool bmqc::annihilation(size_t &u, size_t i) {
    size_t op = 1<<i;
    if (u & op){
        u ^= op;
        return true;
    }else{
        return false;
    }
}


/**
*  Applies the creation operator (without sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool bmqc::creation(size_t &u, size_t i) {
    size_t op = 1<<i;
    if (u & op){
        return false;
    }else{
        u ^= op;
        return true;
    }
}

/**
 * Evaluates the sign change between two (spin separated) base vectors
 * (only for vectors changed by one creation and one annihilation operators)
 */

int bmqc::phase_factor(const size_t& u, const size_t& v){
    int phase_factor = 1;

    auto changed_indices = u ^ v;
    auto first_site_index = __builtin_ctz(changed_indices);
    auto second_site_index = __builtin_ctz(changed_indices>>(first_site_index+1));

    auto no_indices_in_between = second_site_index;

    if (no_indices_in_between == 0 || u == v )  {
        // If there are no indices in between, the phase factor is (+1), since no equal-spin electrons are passed over
    } else {
        // In the region [first_site_index+1, second_site_index+first_site_index], both basis vectors bf1 and bf2 are equal, so it doesn't matter which one we take to take a bitslice

        size_t subbitset = longslice(u, first_site_index + 1, second_site_index + first_site_index);
        if (__builtin_popcount(subbitset) % 2 == 0) {
            // For even number of ones, the phase factor is +1
            phase_factor = 1;
        } else {
            // For odd number of ones, the phase factor is -1
            phase_factor = -1;
        }
    }

    return phase_factor;

}

/**
* Evaluates the sign change for a creation or annihilation operator.
*/
int bmqc::operator_phase(const size_t &u, size_t i) {
    size_t t = u & ((1 << i)-1);
    int q = __builtin_popcount(t);
    if(q%2){
        return -1;
    }else{
        return 1;
    }
}

/**
*  Applies the annihilation operator (with sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/

bool bmqc::annihilation_s(size_t &u, size_t i, int &sign) {
    if(annihilation(u,i)){
        sign *= operator_phase(u,i);
        return true;
    } else{
        return false;
    }
}


/**
*  Applies the creation operator (with sign evaluation)
*  returns bool according to succes.
*
*  WARNING: this does not create a copy of the long, it is altered if succesfull.
*/
bool bmqc::creation_s(size_t &u, size_t i, int &sign) {
    if (creation(u,i)){
        sign *= operator_phase(u,i);
        return true;
    } else{
        return false;
    }

}


