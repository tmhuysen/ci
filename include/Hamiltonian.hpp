#ifndef CI_HAMILTONIAN_HPP
#define CI_HAMILTONIAN_HPP

#include "State.hpp"
#include "StorageType.hpp"

#include <bmqc.hpp>
#include <boost/math/special_functions.hpp>
#include <iostream>
#include <Eigen/Dense>

namespace doci {

class Hamiltonian {
private:
    static Hamiltonian* (*pFactoryMethods[(int)StorageType::STORAGE]) (size_t nbf); //Static array of pointer to factory methods.

    static Hamiltonian* make_dense(size_t nbf); //return a pointer to a new dense instance of the hamiltonian
    static Hamiltonian* make_davidson(size_t nbf); //return a pointer to a new sparse instance of the hamiltonian

protected:
	size_t nbf; //number of basis fucntions in the Hamiltionian aka the dimensions.
	Eigen::VectorXd eigenvalues; //might not keep these in the abstract
	Eigen::MatrixXd eigenvectors; //might not keep these in the abstract
	std::vector<doci::State> groundstates= { doci::State (std::numeric_limits<double>::max(), Eigen::VectorXd()) }; //vector of the all states with the ground energy.

	/**
	 * Adds a @param state to the groundstates vector of our Hamiltonian,
	 * if this state's eigenvalue is equal to the eigenvalue of current groundstates, it is added.
	 * if this state's eigenvalue is lower then it replaces the current groundstates.
	 */
	void groundStates(doci::State state);

public:

	const std::vector<doci::State>& getGroundstates() const;

	//Virtuals
	/**
	 * Adds a value as an element of the hamiltonian matrix.
	 * Virtual: when implementing this class one can opt for many ways to represent the hamiltonian.
	 */
	virtual void add(double value, size_t index1, size_t index2)=0;
	/**
	 *  Solves the eigenvalue problem of the hamiltonian
	 */
	virtual bool solve()=0;

	/**
	 * @param nbf the dimensions of the hamiltonian.
	 * @return Hamiltonian based on the memory it would require.
	 */
	static Hamiltonian* make_hamiltonian(size_t nbf);
	static Hamiltonian* make_hamiltonian(size_t nbf, StorageType type);
};



}  // namespace doci


#endif // CI_HAMILTONIAN_HPP
