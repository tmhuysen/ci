#ifndef CI_BASEHAMILTONIANSOLVER_HPP
#define CI_BASEHAMILTONIANSOLVER_HPP



#include <iostream>

#include <boost/math/special_functions.hpp>
#include <Eigen/Dense>
#include <bmqc.hpp>

#include "State.hpp"



namespace ci {


/**
 *  A base (i.e. abstract) class to solve the eigenvalue problem for the considered Hamiltonian.
 */
class BaseHamiltonianSolver {
protected:
    const size_t dimension;  // the dimension of the Hamiltonian matrix (i.e. the number of occupation number basis functions)

    bool is_solved = false;
    std::vector<ci::State> ground_states;


    // PROTECTED METHODS
    /**
     *  Adds a @param state to the this->ground_states in the following manner:
     *      if the state's eigenvalue is equal to the eigenvalue of current ground states, it is added because then it is a degenerate ground state
     *      if this state's eigenvalue is lower to the eigenvalue of the current ground states, then it replaces the current ground states.
     *
     *  @returns if either of the previously described cases happened: @returns false if no updating was done.
     */
    bool addToGroundStates(const ci::State& state);


public:
    // GETTERS
    std::vector<ci::State> get_ground_states() const { return this->ground_states; }


    // PURE VIRTUAL METHODS
    /**
     *  Solves the eigenvalue problem of the Hamiltonian.
     */
    virtual void solve() = 0;

    /**
     *  Adds @param value to the current value of the Hamiltonian matrix at position (@param index1, @param index2).
     */
    virtual void addElement(double value, size_t index1, size_t index2) = 0;  // TODO: for the Davidson Hamiltonian, this shouldn't be a pure virtual function


//    /**
//     * @param nbf the dimensions of the hamiltonian.
//     * @return Hamiltonian based on the memory it would require.
//     */
//    static Hamiltonian* make_hamiltonian(size_t nbf);
};


}  // namespace ci



#endif // CI_BASEHAMILTONIANSOLVER_HPP
