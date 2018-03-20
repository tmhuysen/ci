#ifndef CI_SOLVERTYPE_HPP
#define CI_SOLVERTYPE_HPP



namespace ci {
namespace solver {

/**
 *  An enum class for the implemented solver types.
 */
enum class SolverType {
    DENSE,
    SPARSE,
    DAVIDSON
};


}  // namespace solver
}  // namespace ci



#endif  // CI_SOLVERTYPE_HPP
