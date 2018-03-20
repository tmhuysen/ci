#ifndef CI_BASESOLVER_HPP
#define CI_BASESOLVER_HPP



namespace ci {
namespace solver {


/**
 *  The abstract and base solver class.
 */
class BaseSolver {
public:

    /*
     *  DESTRUCTOR
     */
    virtual ~BaseSolver() = default;


    /*
     *  PURE VIRTUAL FUNCTIONS
     */
    /**
     *  Solve the eigenvalue problem associated to the type of solver.
     */
    virtual void solve() = 0;
};


}  // namespace solver
}  // namespace ci



#endif  // CI_BASESOLVER_HPP
