#include "BaseHamiltonianSolver.hpp"



namespace ci {

/*
 *  PROTECTED METHODS
 */

/**
 *  Adds a @param state to the this->ground_states in the following manner:
 *      if the state's eigenvalue is equal to the eigenvalue of current ground states, it is added because then it is a degenerate ground state
 *      if this state's eigenvalue is lower to the eigenvalue of the current ground states, then it replaces the current ground states.
 *
 *  @returns if either of the previously described cases happened: @returns false if no updating was done.
 */
bool BaseHamiltonianSolver::addToGroundStates(const ci::State& state) {

    if (state.isDegenerate(this->ground_states[0])) {
        this->ground_states.push_back(state);
        return true;
    } else if (state.hasLowerEnergy(this->ground_states[0])) {
        this->ground_states = std::vector<ci::State> {state};
        return true;
    }

    return false;
}


}  // namespace ci




/*
 * STATIC METHODS
 */

//doci::Hamiltonian *doci::Hamiltonian::make_hamiltonian(size_t nbf) {
//    return new doci::DenseHamiltonian(nbf);
//}
