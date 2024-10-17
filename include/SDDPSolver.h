/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SDDPSolver class, implementing the Solver interface,
 * for multistage linear stochastic programming problems defined by the
 * SDDPBlock. This is a wrapper for the SDDP method implemented by the
 * STochastic OPTimization library (StOpt):
 *
 * https://gitlab.com/stochastic-control/StOpt
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright Copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SDDPSolver
#define __SDDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <Eigen/Dense>
#include "BlockSolverConfig.h"
#include "ScenarioSimulator.h"
#include "SDDPBlock.h"
#include "Solver.h"
#include "StOpt/sddp/OptimizerSDDPBase.h"
#include "StOpt/sddp/SDDPFinalCut.h"

#ifdef USE_MPI
#include <boost/mpi/communicator.hpp>
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

class BendersBFunction;
class StochasticBlock;

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPSolver_CLASSES Classes in SDDPSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SDDPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an SDDP solver for multistage linear stochastic programming problems
/**
 * The SDDPSolver class derives from Solver and implements the stochastic dual
 * dynamic programming (SDDP) method for multistage linear stochastic
 * problems. In fact, this works as a wrapper for the SDDP solver implemented
 * by the STochastic OPTimization (StOpt) library, which is publicly available
 * at https://gitlab.com/stochastic-control/StOpt. This SDDPSolver can be
 * attached to an SDDPBlock whose subproblems are linear, i.e., one that has
 * the following form:
 *
 * \f[
 *    \min_{\substack{x_0 \in \mathbb{R}^{n_0} \\ A_0 x_0 + B_0 x_{-1} = b_0\\
 *                    x_0 \ge 0}} c_0^{\top}x_0 +
 *    \mathbb{E} \left \lbrack
 *    \min_{\substack{x_1 \in \mathbb{R}^{n_1} \\ A_1 x_1 + B_1 x_0 = b_1\\
 *                    x_1 \ge 0}} c_1^{\top}x_1 +
 *    \mathbb{E} \left \lbrack \dots +
 *    \mathbb{E} \left \lbrack
 *    \min_{\substack{x_{T-1} \in \mathbb{R}^{n_{T-1}} \\
 *          A_{T-1} x_{T-1} + B_{T-1} x_{T-2} = b_{T-1}\\
 *                    x_{T-1} \ge 0}} c_{T-1}^{\top}x_{T-1}
 *    \right\rbrack \right\rbrack\right\rbrack,
 * \f]
 *
 * where \f$ T \f$ is called the time horizon and \f$ \xi = \{ (b_t,
 * c_t, A_t, B_t) \}_{t \in \{1, \dots, T-1\}} \f$ is a stochastic
 * process. This means that some (or all) the components of the
 * matrices \f$ A_t \f$ and \f$ B_t \f$ and the vectors \f$ b_t \f$
 * and \f$ c_t \f$ may be random variables. Notice that \f$ x_{-1} \f$
 * and \f$ (b_0, c_0, A_0, B_0) \f$, which we denote by \f$ \xi_0 \f$,
 * are deterministic. The term \f$ B_0 x_{-1} \f$ in the first stage
 * problem could be disregarded (i.e., we could have \f$ B_0 = 0 \f$
 * or \f$ x_{-1} = 0 \f$ without loss of generality), but we keep them
 * in order to have all subproblems with the same structure, which
 * will facilitate our approach.
 *
 * We shall distinguish two types of random variables: the convex and
 * the non-convex random variables.
 *
 * - The *convex random variables* are those that can only appear in
 *   the right-hand side of the constraints, i.e., they can only be
 *   part of the vectors \f$ b_t \f$ for \f$ t \in \{1, \dots, T-1\}
 *   \f$.
 *
 * - The *non-convex random variables* are those that can only appear
 *   in the left-hand side of the constraints or in the objective
 *   function, i.e., they can only be part of the matrices \f$ A_t \f$
 *   or the vectors \f$ c_t \f$ for \f$ t \in \{1, \dots, T-1\} \f$.
 *
 * Among the convex random variables, we further consider two
 * types. The *time-related random variables* are those that depend on
 * some random data of the previous stage. The *time-independent
 * random variables* are those that do not depend on random variables
 * of previous stages. Notice that we do not allow non-convex random
 * variables to be time-related. We also denote the \f$t\f$-th random
 * variable of the stochastic process \f$ \xi \f$ as
 *
 * \f[  \xi_t = ( \omega_t , \nu_t ) \f]
 *
 * where \f$ \omega_t \f$ is a convex random vector and \f$ \nu_t \f$
 * is a non-convex random vector, for \f$ t \in \{1, \dots, T-1\}
 * \f$. The random vector \f$ \omega_t \f$, for \f$ t \in \{1, \dots,
 * T-1\} \f$, is divided into time-related and time-independent random
 * vectors as follows:
 *
 * \f[
 *    \omega_t = ( \omega_t^{\text{dep}} , \omega_t^{\text{ind}} ).
 * \f]
 *
 * This means that the first components of the random vector \f$
 * \omega_t \f$ are time-related random variables and the last
 * components are time-independent random variables. We also define
 * \f$ \omega_{-1}^{\text{dep}} \f$ in order for the first stage
 * problem to (artificially) present the same dependence as the other
 * subproblems do. For each \f$ t \in \{0, \dots, T-1\}\f$, we call
 *
 * \f[
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t}\\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}})
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}}) =
 *      \mathbb{E}
 *        \left\lbrack
 *          V_{t+1}(x_t, \xi_{t+1}) \mid \omega_t^{\text{dep}}
 *        \right\rbrack
 * \f]
 *
 * is the (expected value) cost-to-go function (also called value
 * function, future value function, future cost function), with \f$
 * \mathcal{V}_{T} \equiv 0 \f$ and
 *
 * \f[
 *    V_{t}(x_{t-1}, \xi_{t}) =
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t} \\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{V}_{t+1}(x_t, \omega_t^{\text{dep}})
 * \f]
 *
 * with given \f$ x_{-1} \f$ and (deterministic)
 * \f$ \xi_0\f$. We consider an approximation to the problem
 * associated with stage \f$ t \in \{0, \dots, T-1\} \f$ as the problem
 *
 * \f[
 *    \min_{\substack{x_t \in \mathbb{R}^{n_t} \\
 *                    A_t x_t + B_t x_{t-1} = b_t\\
 *                    x_t \ge 0}} c_t^{\top}x_t +
 *    \mathcal{P}_{t+1}(x_t)
 * \f]
 *
 * where \f$ \mathcal{P}_{t+1}(x_t) \f$ is a polyhedral
 * function, i.e., it is a function of the form
 *
 * \f[
 *    \mathcal{P}_{t+1}(x_t) = \max_{i \in \{1,\dots,k_t\}}
 *                                     \{ d_{t,i}^{\top}x_t + e_{t,i} \}
 * \f]
 *
 * with \f$ d_{t,i} \in \mathbb{R}^{n_t} \f$ and \f$ e_{t,i} \in
 * \mathbb{R} \f$ for each \f$ i \in \{1,\dots,k_t\} \f$.
 */

class SDDPSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// public enum for the possible return values of compute()
 /** Public enum "extending" Solver::sol_type with more detailed values
  * specific to SDDPSolver. */

 enum sol_type_SDDP_S {
  kCurveCross = sol_type::kLastSolverError
  ///< backward and forward curves are crossing
  /**< The convergence of the method is checked every #intNStepConv iterations
   * by computing a "special" forward value. See #intNStepConv for
   * details. The kCurveCross status is returned when the difference between
   * the most recent backward value and the "special" forward value. When the
   * sign of this difference changes with respect to that that was computed
   * for the first time, the method returns this status.
   */
 };  // end( sol_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "int" type that the SDDP solver has, besides those
  * defined in Solver. The value intLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum int_par_type_SDDP_S {

  intNStepConv = int_par_type_S::intLastAlgPar ,
  ///< Frequency at which the convergence is checked
  /**< The method stops when either the maximum number of iterations
   * is reached (an iteration performs one backward pass and one (or
   * two) forward pass(es); see parameter #intMaxIter for the maximum
   * number of iterations) or if the method converges according to
   * the criterion explained in the following. The convergence
   * criterion is checked every intNStepConv iterations. Thus, the
   * method stops if, at some iteration \f$ k \f$ that is a multiple
   * of intNStepConv, the following condition holds:
   *
   *   \f[
   *     \left | \frac{ \bar{z}_{k} - \underline{z}_{k} }
   *                  { \bar{z}_{k} } \right |         \le \epsilon.
   *   \f]
   *
   * where \f$ \epsilon \f$ is the value defined by the #dblAccuracy
   * parameter. \f$ \underline{z}_{k} \f$ is the lower bound obtained
   * at iteration \f$ k \f$, from the backward pass (such a lower
   * bound is computed at every iteration and requires no extra
   * effort). The upper bound \f$ \bar{z}_{k} \f$, on the other hand,
   * is only computed at iterations that are multiple of the value
   * given by the intNStepConv parameter. At such an iteration, an
   * extra forward pass is performed, in which the number of
   * simulations considered is given by the value of the parameter
   * #intNbSimulCheckForConv. \f$ \bar{z}_{k} \f$ is the upper bound
   * computed by this forward pass. The default value for
   * intNStepConv is 1. */

  intPrintTime ,
  ///< Indicates whether computational time should be displayed
  /**< This parameter determines whether the computational time spent
   * at each step of the backward and forward passes should be
   * displayed. If the value of this parameter is zero, the
   * computation time is not displayed. Otherwise, the computational
   * time is displayed at every step. The default value for
   * intPrintTime is 1. */

  intNbSimulCheckForConv ,
  ///< Number of simulations considered when checking convergence
  /**< This parameter determines the number of simulations that must
   * be considered during the forward pass that computes the upper
   * bound used for checking convergence. See the comments for the
   * parameter #intNStepConv for more details. The default value for
   * intNbSimulCheckForConv is 1. */

  intNbSimulBackward ,
  ///< Number of simulations considered in the backward pass
  /**< This parameter determines the number of simulations that must be
   * considered during the backward pass. By default, this number is equal to
   * the number of scenarios. */

  intNbSimulForward ,
  ///< Number of simulations considered in the forward pass
  /**< This parameter determines the number of simulations that must be
   * considered during the forward pass. By default, this number is 1. */

  intOutputFrequency ,
  ///< The frequency at which file outputs are performed
  /**< This parameter determines the frequency at which file outputs are
   * performed (saving the approximations to the future cost functions or
   * serializing an SDDPSolverState). If it is positive, these file outputs
   * are performed every #intOutputFrequency iterations and once at the end of
   * compute() to the files specified by the #strOutputFile, #strStateFile,
   * and #strRandomCutsFile parameters. The default value for this parameter
   * is 0 (i.e., no file output is performed). */

  intFirstStageScenarioId ,
  ///< The id of the scenario to be considered at the first stage
  /**< This parameter specifies the id of the scenario that must be considered
   * while solving the subproblem at the first stage. If it is negative, it
   * means that no scenario must be set while solving the sub-problem at the
   * first stage (i.e., the data for that subproblem has already been set,
   * except possibly the initial state). If it is nonnegative, it must be a
   * number between 0 and the total number of scenarios minus 1. By default,
   * its value is 0. */

  intLastAlgPar
  ///< First allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes
   * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "double" type that the SDDPSolver has, besides those
  * defined in Solver. The value dblLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum dbl_par_type_SDDP_S {

  dblAccuracy = dbl_par_type_S::dblLastAlgPar ,
  ///< relative accuracy for declaring a solution optimal
  /**< The algorithmic parameter for setting the *relative* accuracy required
   * to the solution of the SDDPBlock. Please see the comments of the
   * #intNStepConv parameter for a detailed explanation of its meaning. The
   * default value for dblAccuracy is 1.0e-4. */

  dblLastAlgPar
  ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend
   * the set of double algorithmic parameters. */

 };  // end( dbl_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters
 /** Public enum describing the different types of algorithmic
  * parameters of "string" type that the SDDPSolver has, besides those
  * defined in Solver. The value strLastAlgPar is provided so that the
  * list can be easily further extended by derived classes. */

 enum str_par_type_SDDP_S {

  strRegressorsFilename = str_par_type_S::strLastAlgPar ,
  ///< name of the file that will store regressors
  /**< Name of the file in which the regressors will be stored. The
   * default value is "regressors.sddp". */

  strCutsFilename ,
  ///< name of the file that will store the cuts
  /**< Name of the file in which the cuts will be stored. The
   * default value is "cuts.sddp". */

  strVisitedStatesFilename ,
  ///< name of the file that will store the visited states
  /**< Name of the file in which the visited states will be stored.
   * The default value is "visited_states.sddp". */

  strInnerBC ,
  ///< name of the file containing the default BlockConfig for the inner Block
  /**< Name of the file containing the default BlockConfig that will be
   * applied to the inner Block of each BendersBFunction. By default, this is
   * empty. */

  strInnerBSC ,
  ///< name of the file containing the default BlockSolverConfig for inner Block
  /**< Name of the file containing the default BlockSolverConfig that will be
   * applied to the inner Block of each BendersBFunction. By default, this is
   * empty. */

  strOutputFile ,
  ///< name of the file to which the future cost functions will be output
  /**< Name of the file to which the approximations to the future cost
   * functions are output. See #intOutputFrequency for controlling if and when
   * these approximations are output. By default, the name of this file is
   * empty, which means that the future cost functions will not be output. */

  strStateFile ,
  ///< name of the file in which the SDDPSolverState will be serialized
  /**< Name of the file in which the SDDPSolverState will be serialized. See
   * #intOutputFrequency for controlling if and when these approximations are
   * output. By default, the name of this file is empty, which means that no
   * state will be serialized. */

  strRandomCutsFile ,
  ///< name of the file to which the random cuts will be output
  /**< A random cut is a cut associated with a particular scenario. This
   * parameter indicates the name of the file to which the random cuts will be
   * output. See #intOutputFrequency for controlling if and when the random
   * cuts are output. See serialize_random_cuts() for a description of the
   * format of the output file. By default, the name of this file is empty,
   * which means that the random cuts are not output (and not stored). */

  strFilenameSuffix ,
  ///< suffix to be added to a filename every other iteration
  /**< This is the suffix that will be added to a non-empty filename (given by
   * #strOutputFile, #strStateFile, and #strRandomCutsFile) every other
   * iteration in which a file output is performed (see
   * #intOutputFrequency). This parameter can be useful, for instance, if this
   * Solver is running on an unreliable system, which may crash while the
   * output is being performed. By using a suffix, at least some not so old
   * data will be available. For instance, suppose that #intOutputFrequency >
   * 0, #strStateFile = "state.nc4", and #strSuffix = ".0". Then, the first
   * time the State is serialized, it will be serialized in the file called
   * "state.nc4". The second time, it will be serialized into
   * "state.nc4.0". The third time it will be serialized again into
   * "state.nc4" and so on. By default, this is empty. */

  strSubSolverLogFilePrefix ,
  ///< prefix of the names of the files for the logs of the sub-Solvers
  /**< Prefix of the names of the files in which the logs of the sub-Solvers
   * will be output. For instance, suppose this prefix is "logfile-". Then, if
   * the SDDPBlock has a single sub-Block at each stage, the log of the
   * sub-Solver associated with time instant t will be output into a file
   * called "logfile-t". If the SDDPBlock has n > 1 sub-Blocks at each stage,
   * then the log of the sub-Solver associated with time instant t and i-th
   * sub-Block at stage t will be output into a file called "logfile-t-i". By
   * default, this is empty. */

  strDirOUT,
  ///< path where the results are written
  /**< this can be empty
   * in that case results are written where solver is launched */
  
  strLastAlgPar
  ///< first allowed new string parameter for derived classes
  /**< Convenience value for easily allow derived classes
   * to extend the set of string algorithmic parameters. */

 };  // end( str_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-int parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-int type that SDDPSolver has in addition to these of
  * Solver. The value vintLastAlgPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum vint_par_type_SDDP_S {

  vintMeshDiscretization = vint_par_type_S::vintLastAlgPar ,
  ///< number of meshes in each direction
  /**< The algorithmic parameter for setting the discretization of the meshes,
   * which is a vector containing the number of meshes (number of steps) in
   * each direction. The size of this vector must be equal to the dimension of
   * a simulation particle and its i-th component contains the number of
   * meshes (number of steps) at the i-th direction. */

  vintLastAlgPar
  ///< first allowed new vector-of-int parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * vector-of-int parameters. */

 };  // end( vint_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-double parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-double type that SDDPSolver has in addition to these of
  * Solver. The value vdblLastAlgPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum vdbl_par_type_SDDP_S {

  vdblLastStageCuts = vdbl_par_type_S::vdblLastAlgPar ,
  ///< the cuts to be used at the last stage
  /**< The parameter for setting the cuts to be used at the last time
  * instant. Each cut is represented by a vector whose size is equal to m + 1,
  * where m is the number of state variables. The first m elements of a cut
  * are the coefficients for the state variables and the last element is the
  * constant term of that cut. The vector #vdblLastStageCuts can store
  * multiple cuts and its size must be a multiple of m + 1. If it is non-empty
  * and has size k * (m + 1), then it contains k cuts and the i-th cut is
  * given by the elements between the indices i * ( m + 1 ) and ( i + 1 ) * (
  * m + 1 ) - 1. By default, this vector is empty. */

  vdblInitialState ,
  ///< the initial state for the first stage problem
  /**< The parameter for setting the initial state, i.e., the initial state to
   * be considered in the subproblem of the first stage. The size of this
   * vector must be equal to the size of the initial state and the i-th
   * element in this vector will be the value of the i-th initial state
   * variable of the first stage subproblem. If this vector is empty, no
   * initial state is set for the first stage problem. By default, this vector
   * is empty. */

  vdblLastAlgPar
  ///< first allowed new vector-of-double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * vector-of-double parameters. */

 };  // end( vdbl_par_type_SDDP_S )

/*--------------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend class SDDPOptimizer;
 friend class SDDPSolverState;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPSolver -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPSolver
 *  @{ */

 /// constructor
 SDDPSolver( void ) {

  v_events.resize( max_event_number() );

  // Set the parameters to their default values

  set_default_parameters();

  // SDDPOptimizer

  sddp_optimizer = std::make_shared< SDDPOptimizer >( this );
 }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPSolver() {
  delete f_inner_block_config;
  delete f_inner_block_solver_config;
  delete f_get_var_solution_config;
 }

/*--------------------------------------------------------------------------*/

 void set_Block( Block * block ) override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set a given integer (int) numerical parameter
 /** Set a given integer (int) numerical parameter. Besides
  * considering the integer parameters defined in #int_par_type_S,
  * this function also accepts the following parameters:
  *
  * - #intMaxIter
  *
  * - #intNStepConv
  *
  * - #intPrintTime
  *
  * - #intNbSimulCheckForConv
  *
  * - #intNbSimulBackward
  *
  * - #intNbSimulForward
  *
  * - #intFirstStageScenarioId
  *
  * Please refer to the #int_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const int value ) override {
  switch( par ) {
   case( intMaxIter ): maximum_number_iterations = value; return;
   case( intNStepConv ): convergence_frequency = value; return;
   case( intPrintTime ): print_cpu_time = value; return;
   case( intNbSimulCheckForConv ):
    number_simulations_for_convergence = value; return;
   case( intNbSimulBackward ):
    sddp_optimizer->set_number_simulations_backward( value ); return;
   case( intNbSimulForward ):
    sddp_optimizer->set_number_simulations_forward( value ); return;
   case( intLogVerb ): log_verbosity = value; return;
   case( intOutputFrequency ): output_frequency = value; return;
   case( intFirstStageScenarioId ):
    first_stage_scenario_index = value; return;
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set a given float (double) numerical parameter
 /** Set a given float (double) numerical parameter. Besides
  * considering the integer parameters defined in #dbl_par_type_S,
  * this function also accepts the following parameters:
  *
  * - #dblAccuracy
  *
  * Please refer to the #dbl_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const double value ) override {
  if( par == dblAccuracy ) {
   accuracy = value;
   return;
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set a given string parameter
 /** Set a given string parameter. Besides considering the integer
  * parameters defined in #str_par_type_S, this function also accepts
  * the following parameters:
  *
  * - #strRegressorsFilename
  *
  * - #strCutsFilename
  *
  * - #strVisitedStatesFilename
  *
  * - #strInnerBC
  *
  * - #strInnerBSC
  *
  * - #strOutputFile
  *
  * - #strStateFile
  *
  * - #strRandomCutsFile
  *
  * - #strFilenameSuffix
  *
  * - #strSubSolverLogFilePrefix
  *
  * - #strDirOUT [""]: path where results are written
  *
  * Please refer to the #str_par_type_SDDP_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const std::string & value ) override {
  switch( par ) {
   case( strRegressorsFilename ): regressors_filename = value; return;
   case( strCutsFilename ): cuts_filename = value; return;
   case( strVisitedStatesFilename ): visited_states_filename = value; return;
   case( strInnerBC ): f_inner_block_config_filename = value; return;
   case( strInnerBSC ): f_inner_block_solver_config_filename = value; return;
   case( strOutputFile ): f_output_filename = value; return;
   case( strStateFile ): f_state_filename = value; return;
   case( strRandomCutsFile ): f_random_cuts_filename = value; return;
   case( strFilenameSuffix ): f_filename_suffix = value; return;
   case( strDirOUT): f_dir_out_pathname = value; return;
   case( strSubSolverLogFilePrefix ): {
    f_sub_solver_filename_prefix = value;
    return;
   }
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-int paramaters of SDDPSolver
 /** Set a given vector-of-int paramater. Besides considering the vector-of-int
  * parameters defined in #vint_par_type_S, this function also accepts the
  * following parameters:
  *
  * - #vintMeshDiscretization
  *
  * Please refer to the #vint_par_type_SDDP_S enumeration for a detailed
  * description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( idx_type par , std::vector< int > && value ) override {
  if( par == vintMeshDiscretization ) {
   mesh_discretization = value;
   return;
  }
  Solver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the vector-of-double paramaters of SDDPSolver
 /** Set a given vector-of-double paramater. Besides considering the
  * vector-of-double parameters defined in #vdbl_par_type_S, this function
  * also accepts the following parameters:
  *
  * - #vdblLastStageCuts
  *
  * - #vdblInitialState
  *
  * Please refer to the #vdbl_par_type_SDDP_S enumeration for a detailed
  * description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( idx_type par , std::vector< double > && value ) override {
  switch( par ) {
   case( vdblLastStageCuts ): last_stage_cuts = value; return;
   case( vdblInitialState ): initial_state = value; return;
  }
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

 /// set the whole set of parameters of this SDDPSolver in one blow
 /** This method sets the whole set of parameters of this SDDPSolver in one
  * blow using a ComputeConfig object.
  *
  * Besides considering all the parameters of an SDDPSolver, it can also be
  * used to configure the inner Block of every BendersBFunction by means of
  * the extra Configuration (ComputeConfig::f_extra_Configuration). If the
  * pointer to the extra Configuration is not nullptr, it can be any of the
  * following:
  *
  * - a pointer to a BlockConfig, which will be used to configure the inner
  *   Block of the BendersBFunction at every stage;
  *
  * - a pointer to a BlockSolverConfig, which will be used to configure the
  *   Solver of the inner Block of the BendersBFunction at every stage;
  *
  * - a pointer to a SimpleConfiguration< std::vector< Configuration * > >.
  *
  * In the last case, the first element of the vector, if present, must be
  * either nullptr or a pointer to a BlockConfig. The second element, if
  * present, must be either nullptr or a pointer to a BlockSolverConfig. These
  * will be used to configure the inner Block of the BendersBFunction at every
  * stage and their Solver. Finally, the third element, if present, must be
  * either nullptr or a pointer to a Configuration. This Configuration will be
  * used to retrieve the Solution from the inner Block of the
  * BendersBFunction, at every stage, after it is solved. This Configuration
  * will be passed to get_var_solution() of the inner Solver. The relevant
  * part of the Solution of the inner Block is the values of the active
  * Variables of the PolyhedralFunction. Thus, this Configuration can be used
  * to specify that only that portion of the Solution should be retrieved.
  *
  * If the extra Configuration is not any of the specified above, an exception
  * is thrown.
  *
  * Here, we are assuming that the same Configuration can be applied to the
  * inner Block of the BendersBFunction at all stages. However, in principle,
  * the inner Block of the BendersBFunction at different stages could require
  * different Configuration. If this case ever happens, the implementation of
  * this method should be adapted to take it into consideration.
  *
  * If the given pointer to the ComputeConfig \p scfg is nullptr, then the
  * Configuration of this SDDPSolver is reset to its default one.
  *
  * It is important to notice that every Configuration provided by \p scfg is
  * cloned (see Configuration::clone()) and, therefore, the caller is
  * responsible for destroying all these Configuration and the Configuration
  * pointed by \p scfg.
  *
  * @param scfg a pointer to a ComputeConfig.
  */

 void set_ComputeConfig( ComputeConfig *scfg = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the SDDPSolver
 *  @{ */

 /// get the number of int parameters
 /** Get the number of int parameters.
  *
  * @return The number of int parameters.
  */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of double parameters
 /** Get the number of double parameters.
  *
  * @return The number of double parameters.
  */

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of string parameters
 /** Get the number of string parameters.
  *
  * @return The number of string parameters.
  */

 idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of vector-of-int parameters
 /** Get the number of vector-of-int  parameters.
  *
  * @return The number of vector-of-int parameters.
  */

 idx_type get_num_vint_par( void ) const override {
  return( idx_type( vintLastAlgPar ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of vector-of-double parameters
 /** Get the number of vector-of-double  parameters.
  *
  * @return The number of vector-of-double parameters.
  */

 idx_type get_num_vdbl_par( void ) const override {
  return( idx_type( vdblLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of an int parameter
 /** Get the default value of the int parameter with given index.  Please see
  * the #int_par_type_SDDP_S and #int_par_type_S enumerations for a detailed
  * explanation of the possible parameters. This function returns the
  * following values depending on the desired parameter:
  *
  * - #intNStepConv: 1
  *
  * - #intPrintTime: 1
  *
  * - #intNbSimulCheckForConv: 1
  *
  * - #intNbSimulBackward: given by
  *   SDDPOptimizer::get_dflt_number_simulations_backward()
  *
  * - #intNbSimulForward: 1
  *
  * - #intLogVerb: 0
  *
  * - #intOutputFrequency: 0
  *
  * - #intFirstStageScenarioId: 0
  *
  * For any other parameter, see Solver::get_dflt_int_par().
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 int get_dflt_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intNStepConv ): return 1;
   case( intPrintTime ): return 1;
   case( intNbSimulCheckForConv ): return 1;
   case( intNbSimulBackward ):
    return sddp_optimizer->get_dflt_number_simulations_backward();
   case( intNbSimulForward ): return 1;
   case( intLogVerb ): return 0;
   case( intOutputFrequency ): return 0;
   case( intFirstStageScenarioId ): return 0;
  }
  return Solver::get_dflt_int_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the default value of a double parameter
 /** Get the default value of the double parameter with given index.
  * Please see the #dbl_par_type_SDDP_S and #dbl_par_type_S
  * enumerations for a detailed explanation of the possible
  * parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 double get_dflt_dbl_par( const idx_type par ) const override {
  if( par == dblAccuracy ) return 1.0e-4;
  return Solver::get_dflt_dbl_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the default value of a string parameter
 /** Get the default value of the string parameter with given index.
  * Please see the #str_par_type_SDDP_S and #str_par_type_S
  * enumerations for a detailed explanation of the possible
  * parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::string & get_dflt_str_par( const idx_type par ) const override {

  static const std::vector< std::string > default_values =
   { "regressors.sddp" , "cuts.sddp" , "visited_states.sddp" ,
     "", "" , "" , "" , "" , "" , "" , "" };

  if( par >= str_par_type_S::strLastAlgPar && par < strLastAlgPar )
   return default_values[ par - str_par_type_S::strLastAlgPar ];

  return Solver::get_dflt_str_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of a vector-of-int parameter
 /** Get the default value of the vector-of-int parameter with given index.
  * Please see the #vint_par_type_SDDP_S and #vint_par_type_S enumerations for
  * a detailed explanation of the possible parameters. This function returns
  * the following values depending on the desired parameter:
  *
  * - #vintMeshDiscretization: an empty vector
  *
  * For any other parameter, see Solver::get_dflt_vint_par().
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::vector< int > & get_dflt_vint_par( const idx_type par )
  const override {
  const static std::vector< int > empty;

  if( par == vintMeshDiscretization ) {
   return empty;
  }

  return Solver::get_dflt_vint_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of a vector-of-double parameter
 /** Get the default value of the vector-of-double parameter with given index.
  * Please see the #vdbl_par_type_SDDP_S and #vdbl_par_type_S enumerations for
  * a detailed explanation of the possible parameters. This function returns
  * the following values depending on the desired parameter:
  *
  * - #vdblMeshDiscretization: an empty vector
  *
  * - #vdblInitialState: an empty vector
  *
  * For any other parameter, see Solver::get_dflt_vdbl_par().
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::vector< double > & get_dflt_vdbl_par( const idx_type par )
  const override {
  const static std::vector< double > empty;

  if( par == vdblLastStageCuts || par == vdblInitialState ) {
   return empty;
  }

  return Solver::get_dflt_vdbl_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get a specific integer (int) numerical parameter
 /** Get a specific integer (int) numerical parameter. Please see the
  * #int_par_type_SDDP_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 int get_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intMaxIter ): return maximum_number_iterations;
   case( intNStepConv ): return convergence_frequency;
   case( intPrintTime ): return print_cpu_time;
   case( intNbSimulCheckForConv ): return number_simulations_for_convergence;
   case( intNbSimulBackward ):
    return sddp_optimizer->get_number_simulations_backward();
   case( intNbSimulForward ):
    return sddp_optimizer->get_number_simulations_forward();
   case( intLogVerb ): return log_verbosity;
   case( intOutputFrequency ): return output_frequency;
   case( intFirstStageScenarioId ): return first_stage_scenario_index;
  }
  return( Solver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific float (double) numerical parameter
 /** Get a specific float (double) numerical parameter. Please see the
  * #dbl_par_type_SDDP_S and #dbl_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 double get_dbl_par( const idx_type par ) const override {
  if( par == dblAccuracy )
   return accuracy;
  return( get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific string numerical parameter
 /** Get a specific string numerical parameter. Please see the
  * #str_par_type_SDDP_S and #str_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::string & get_str_par( const idx_type par ) const override {
  switch( par ) {
   case( strRegressorsFilename ): return regressors_filename;
   case( strCutsFilename ): return cuts_filename;
   case( strVisitedStatesFilename ): return visited_states_filename;
   case( strInnerBC ): return f_inner_block_config_filename;
   case( strInnerBSC ): return f_inner_block_solver_config_filename;
   case( strOutputFile ): return f_output_filename;
   case( strStateFile ): return f_state_filename;
   case( strRandomCutsFile ): return f_random_cuts_filename;
   case( strFilenameSuffix ): return f_filename_suffix;
   case( strSubSolverLogFilePrefix ): return f_sub_solver_filename_prefix;   
   case( strDirOUT ): return f_dir_out_pathname;
  }
  return Solver::get_str_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific vector-of-int parameter
 /** Get a specific vector-of-int parameter. Please see the
  * #vint_par_type_SDDP_S and #vint_par_type_S enumerations for a detailed
  * explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::vector< int > & get_vint_par( const idx_type par ) const override {
  if( par == vintMeshDiscretization )
   return mesh_discretization;
  return Solver::get_vint_par( par );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get a specific vector-of-double parameter
 /** Get a specific vector-of-double parameter. Please see the
  * #vdbl_par_type_SDDP_S and #vdbl_par_type_S enumerations for a detailed
  * explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::vector< double > & get_vdbl_par( const idx_type par )
  const override {
  switch( par ) {
   case( vdblLastStageCuts ): return last_stage_cuts;
   case( vdblInitialState ): return initial_state;
  }
  return Solver::get_vdbl_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with given string \p name
 /** This method takes a string, which is assumed to be the name of an int
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it. The method is given a void
  * implementation (throwing exception), rather than being pure virtual, so
  * that derived classes not having any int parameter do not have to bother
  * with implementing it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "intNStepConv" ) return intNStepConv;
  if( name == "intPrintTime" ) return intPrintTime;
  if( name == "intNbSimulCheckForConv" ) return intNbSimulCheckForConv;
  if( name == "intNbSimulBackward" ) return intNbSimulBackward;
  if( name == "intNbSimulForward" ) return intNbSimulForward;
  if( name == "intOutputFrequency" ) return intOutputFrequency;
  if( name == "intFirstStageScenarioId" ) return intFirstStageScenarioId;
  return Solver::int_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the double parameter with given string name
 /** This method takes a string, which is assumed to be the name of a double
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  if( name == "dblAccuracy" ) return dblAccuracy;
  return Solver::dbl_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the string parameter with given string name
 /** This method takes a string, which is assumed to be the name of a string
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type str_par_str2idx( const std::string & name ) const override {
  if( name == "strRegressorsFilename" ) return strRegressorsFilename;
  if( name == "strCutsFilename" ) return strCutsFilename;
  if( name == "strVisitedStatesFilename" ) return strVisitedStatesFilename;
  if( name == "strInnerBC" ) return strInnerBC;
  if( name == "strInnerBSC" ) return strInnerBSC;
  if( name == "strOutputFile" ) return strOutputFile;
  if( name == "strStateFile" ) return strStateFile;
  if( name == "strRandomCutsFile" ) return strRandomCutsFile;
  if( name == "strFilenameSuffix" ) return strFilenameSuffix;
  if( name == "strSubSolverLogFilePrefix" ) return strSubSolverLogFilePrefix;
  if( name == "strDirOUT" ) return strDirOUT;
  return Solver::str_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the vector-of-int parameter with given string name
 /** This method takes a string, which is assumed to be the name of a
  * vector-of-int parameter, and returns its index, i.e., the integer value
  * that can be used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type vint_par_str2idx( const std::string & name ) const override {
  if( name == "vintMeshDiscretization" ) return vintMeshDiscretization;
  return Solver::vint_par_str2idx( name );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the index of the vector-of-double parameter with given string name
 /** This method takes a string, which is assumed to be the name of a
  * vector-of-double parameter, and returns its index, i.e., the double value
  * that can be used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type vdbl_par_str2idx( const std::string & name ) const override {
  if( name == "vdblLastStageCuts" ) return vdblLastStageCuts;
  if( name == "vdblInitialState" ) return vdblInitialState;
  return Solver::vdbl_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/
 /// returns the string name of the int parameter with given index
 /** This method takes an int parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & int_par_idx2str( const idx_type idx ) const override {

  static const std::vector< std::string > parameter_names =
   { "intNStepConv", "intPrintTime", "intNbSimulCheckForConv" ,
     "intNbSimulBackward" , "intNbSimulForward" , "intOutputFrequency" ,
     "intFirstStageScenarioId"  };

  if( idx >= int_par_type_S::intLastAlgPar && idx < intLastAlgPar )
   return parameter_names[ idx - int_par_type_S::intLastAlgPar ];

  return Solver::int_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the double parameter with given index
 /** This method takes a double parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & dbl_par_idx2str( const idx_type idx ) const override {
  static const std::string dblAccuracy_name = "dblAccuracy";
  if( idx == dblAccuracy ) return dblAccuracy_name;
  return Solver::dbl_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the string parameter with given index
 /** This method takes a string parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & str_par_idx2str( const idx_type idx ) const override {

  static const std::vector< std::string > parameter_names =
   { "strRegressorsFilename", "strCutsFilename", "strVisitedStatesFilename" ,
     "strInnerBC" , "strInnerBSC" , "strOutputFile" , "strStateFile" ,
     "strRandomCutsFile", "strFilenameSuffix" , "strSubSolverLogFilePrefix",
     "strDirOUT"	 };

  if( idx >= str_par_type_S::strLastAlgPar && idx < strLastAlgPar )
   return parameter_names[ idx - str_par_type_S::strLastAlgPar ];

  return Solver::str_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the vector-of-int parameter with given index
 /** This method takes a vector-of-int parameter index, i.e., the integer
  * value that can be used in [set/get]_par() [see above] to set/get it, and
  * returns its "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & vint_par_idx2str( const idx_type idx ) const override {
  static const std::vector< std::string > parameter_names =
   { "vintMeshDiscretization" };
  if( idx >= vint_par_type_S::vintLastAlgPar && idx < vintLastAlgPar )
   return parameter_names[ idx - vint_par_type_S::vintLastAlgPar ];
  return Solver::vint_par_idx2str( idx );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the string name of the vector-of-double parameter with given index
 /** This method takes a vector-of-double parameter index, i.e., the double
  * value that can be used in [set/get]_par() [see above] to set/get it, and
  * returns its "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & vdbl_par_idx2str( const idx_type idx ) const override {
  static const std::vector< std::string > parameter_names =
   { "vdblLastStageCuts" , "vdblInitialState" };
  if( idx >= vdbl_par_type_S::vdblLastAlgPar && idx < vdblLastAlgPar )
   return parameter_names[ idx - vdbl_par_type_S::vdblLastAlgPar ];
  return Solver::vdbl_par_idx2str( idx );
 }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Set event handlers
 *
 *  SDDPSolver manages the following events:
 *
 * - eEverykIteration, called just after a new forward pass begins.
 *
 * Events have to be set with set_event_handler() for them to be called.
 * @{ */

 /// register a new event handler, returning its id
 /** The new event handler is added at the back of v_events[ type ]. As the &&
  * tells, the event handler becomes property of the SDDPSolver, which is
  * completely OK if, as one expects, it is defined via a lambda function. The
  * method returns a unique id for the handler, which can (and must) be later
  * used to remove the handler before it becomes invalid. Note that the
  * handler is type-specific, i.e., two event handlers of different types can
  * have the same id; in other words, the "real" id is the pair ( type , id
  * ). An exception is thrown if the SDDPSolver is not capable of handling
  * this type or event for whatever reason, among which that it has exhausted
  * the available maximum number of event handlers slots for the given
  * type. */

 EventID set_event_handler( int type , EventHandler && event ) override {
  if( type != eEverykIteration )
   throw( std::invalid_argument( "SDDPSolver:set_event_handler: unsupported "
                                 "event type " + std::to_string( type ) ) );

  if( v_events[ type ].size() > std::numeric_limits< EventID >::max() )
   throw( std::invalid_argument( "SDDPSolver:set_event_handler: too many event "
                                 "handlers for type" + std::to_string( type ) ) );

  EventID id = v_events[ type ].size();
  v_events[ type ].push_back( std::move( event ) );

  return id ;
 }

/*--------------------------------------------------------------------------*/

 /// unregister an existing event handler
 /** Removes the event handler with the given id from the list of those
  * registered for the given type. If there is no event handler with the given
  * id for the given type, exception will be thrown. */
 void reset_event_handler( int type , EventID id ) override;

/*--------------------------------------------------------------------------*/

 /// returns the maximum number of event types supported by the SDDPSolver
 [[nodiscard]] virtual EventID max_event_number() const override {
  return e_last_event_type;
 }

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR HANDLING THE State OF THE SDDPSolver ------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the State of the SDDPSolver
 *  @{ */

 State * get_State( void ) const override;

/*--------------------------------------------------------------------------*/

 void put_State( const State & state ) override;

/*--------------------------------------------------------------------------*/

 void put_State( State && state ) override;

/*--------------------------------------------------------------------------*/

 void serialize_State( netCDF::NcGroup & group ,
		       const std::string & sub_group_name = "" )
  const override;

/*--------------------------------------------------------------------------*/

 /// serialize the State of this SDDPSolver
 /** This function serializes the State of this SDDPSolver in the file with
  * the given name. If \p filename is empty, then no operation is performed.
  *
  * @param filename The name of the file in which the State of this SDDPSolver
  *        will be serialized. */

 void serialize_State( const std::string & filename ) const;

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// (try to) solve the model encoded in the SDDPBlock
 /**
  */

 int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 /// this function has an empty implementation (and thus does nothing)
 /** This function has an empty implementation (and thus does nothing). The
  * solution produced by the SDDPSolver is in the form of cuts to the
  * PolyhedralFunction of the SDDPBlock. These cuts are added during
  * compute() and therefore no task needs to be performed here. */

 void get_var_solution( Configuration *solc = nullptr ) override { }

/*--------------------------------------------------------------------------*/

 double get_lb( void ) override;

/*--------------------------------------------------------------------------*/

 double get_ub( void ) override;

/*--------------------------------------------------------------------------*/

 double get_backward_value( void ) const {
  return backward_value;
 }

/*--------------------------------------------------------------------------*/

 double get_forward_value( void ) const {
  return forward_value;
 }

/*--------------------------------------------------------------------------*/

 /// output the future cost functions
 /** This function outputs the approximations to the future cost functions to
  * the file with the given name. If \p filename is empty, then no operation
  * is performed.
  *
  * @param filename The name of the file to which the functions will be
  *        output. */

 void output_future_cost_functions( const std::string & filename ) const;

/*--------------------------------------------------------------------------*/

 /// output the cuts and the SDDPSolverState to files
 /** This function outputs the approximations to the future cost functions (if
  * #strOutputFile is non-empty) and the random cuts (if #strRandomCutsFile is
  * non-empty), and serializes the SDDPSolverState (if #strStateFile is
  * non-empty). If #strFilenameSuffix is non-empty, then it is added to the
  * filename every other iteration that this function is called. */

 void file_output() const;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE SDDPSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the state of the SDDPSolver
 *  @{ */

/*--------------------------------------------------------------------------*/

 /// returns the number of iterations performed in the last call to compute()
 /** This function returns the number of iterations performed by the solver
  * during the last call to compute().
  *
  * @return The number of iterations performed by the solver during the last
  *         call to compute(). */

 int get_number_iterations_performed() const {
  return number_iterations_performed;
 }

/*--------------------------------------------------------------------------*/

 /// returns the time horizon of the problem associated with the SDDPBlock
 /** This function returns the time horizon of the problem associated with the
  * SDDPBlock with which this SDDPSolver is attached.
  *
  * @return The time horizon of the problem associated with the SDDPBlock. */

 SDDPBlock::Index get_time_horizon() const;

/*--------------------------------------------------------------------------*/

/// returns the solution associated with the subproblem at the given stage
/** This function returns the solution of the subproblem associated with the
 * given \p stage, which is part of the state variables of the next stage.
 *
 * @param stage The stage whose solution is required.
 *
 * @return The array containing the solution of the subproblem at the given
 *         stage. */

 template< class T = Eigen::ArrayXd >
 T get_solution( SDDPBlock::Index stage ,
                 SDDPBlock::Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 /// returns the index of the scenario for the first stage (if any)
 /** The first stage problem may not have any scenario associated with it. In
  * this case, this function returns Inf< Index >(). Otherwise, it returns the
  * index of the scenario that must be considered in the first stage.
  *
  * @return The index of the scenario to be considered in the first stage (if
  *         any). */

 Index get_first_stage_scenario_index() const {
  if( first_stage_scenario_index < 0 )
   return Inf< Index >();
  return first_stage_scenario_index;
 }

/*--------------------------------------------------------------------------*/

 /// indicates whether random cuts must be stored
 /** A random cut is a cut associated with a particular scenario. This
  * function indicates whether the random cuts that are produced must be
  * stored in the SDDPBlock.
  *
  * @return true if and only if random cuts must be stored in the
  *         SDDPBlock. */

 bool store_random_cuts() const {
  return ! f_random_cuts_filename.empty();
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// returns a pointer to the BendersBFunction associated with the given stage
 /** This function returns a pointer to the BendersBFunction associated with
  * the given \p stage, which must be an integer between 0 and
  * get_time_horizon() - 1.
  *
  * @param stage An index between 0 and get_time_horizon() - 1. */

 BendersBFunction * get_benders_function
 ( SDDPBlock::Index stage , SDDPBlock::Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 /// returns true if a valid mesh discretization has been provided
 bool mesh_provided() const {
  return( ( ! mesh_discretization.empty() ) &&
          std::all_of( mesh_discretization.begin() ,
                       mesh_discretization.end() ,
                       []( auto i ) { return( i > 0 ); } ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED CLASSES ----------------------------*/
/*--------------------------------------------------------------------------*/

 class SDDPOptimizer : public StOpt::OptimizerSDDPBase {

 public:

  /// constructor taking a pointer to an SDDPSolver
  /** Constructs an SDDPOptimizer associated with the given SDDPSolver.
   *
   * @param solver A pointer to an SDDPSolver. This parameter is optional and
   *        its default value is nullptr. */

  SDDPOptimizer( SDDPSolver * solver = nullptr ) {
   sddp_solver = solver;
   simulator_backward = std::make_shared< ScenarioSimulator >( true );
   simulator_forward = std::make_shared< ScenarioSimulator >( false );
  }

/*--------------------------------------------------------------------------*/

  /// prepares this SDDPOptimizer for another call to the StOpt SDDP solver
  /** This function prepares this SDDPOptimizer for another call to the StOpt
   * SDDP solver. This function must be called within SDDPSolver::compute()
   * before StOpt is invoked. */

  virtual void reset() {
   current_iteration = 0;
   previous_pass_was_backward = false;
  }

/*--------------------------------------------------------------------------*/

  Eigen::ArrayXd oneStepBackward
  ( const StOpt::SDDPCutOptBase & sddp_cut ,
    const std::tuple< std::shared_ptr< Eigen::ArrayXd > , int , int > & state ,
    const Eigen::ArrayXd & particle , const int & simulation_id ,
    const Index scenario_index , const bool scenario_must_be_set ,
    const Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

  Eigen::ArrayXd oneStepBackward
  ( const StOpt::SDDPCutOptBase & p_linCut,
    const std::tuple< std::shared_ptr< Eigen::ArrayXd >, int, int > & p_aState,
    const Eigen::ArrayXd & p_particle, const int & p_isample) const override;

/*--------------------------------------------------------------------------*/

  double oneStepForward
  ( const Eigen::ArrayXd & particle , Eigen::ArrayXd & state ,
    Eigen::ArrayXd & state_to_store , const StOpt::SDDPCutOptBase & sddp_cut ,
    const int & simulation_id , const Index scenario_index ,
    const bool scenario_must_be_set , const Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

  double oneStepForward
  ( const Eigen::ArrayXd &p_aParticle , Eigen::ArrayXd &p_state ,
    Eigen::ArrayXd &p_stateToStore , const StOpt::SDDPCutOptBase &p_linCut ,
    const int &p_isimu ) const override;

/*--------------------------------------------------------------------------*/

  /// updates this SDDPOptimizer for a new stage
  /** This function updates this SDDPOptimizer for a new stage. The \p date
   * and \p date_next parameters have different meanings in the backward and
   * forward steps:
   *
   * - In a backward step, the optimization subproblem that must be solved is
   *   that associated with the stage \p date_next. The parameter \p date
   *   indicates, therefore, the previous stage.
   *
   * - In a forward step, the optimization subproblem that must be solved is
   *   that associated with the stage \p date. The parameter \p date_next
   *   indicates, therefore, the next stage.
   *
   *   We assume that the given arguments correspond to consecutive stages,
   *   i.e., \p date_next == \p date + 1. Moreover, we assume that -1 <= \p
   *   date < T, where T is the time horizon. Also, the following conditions
   *   should hold:
   *
   * -# -1 <= \p date <= T - 2 for any backward step;
   *
   * -#  0 <= \p date <= T - 1 for any forward step.
   *
   * @param date A stage.
   *
   * @param date_next Another stage. */

  void updateDates( const double & date, const double & date_next ) override;

/*--------------------------------------------------------------------------*/

  /// returns an initial state for the subproblem at the given stage
  /** This function must return an initial state for the optimization
   * subproblem associated with the given date. If the given date is t, then
   * the optimization subproblem associated with time t has variables x_t and
   * depends on the state (x_{t-1}, w_{t-1}^{dep}). Thus, this function must
   * return an array containing values for x_{t-1} and w_{t-1}^{dep}. Notice
   * that the subvector w_{t-1}^{dep} is only present in the state if some
   * random variables of the subproblem associated with time t depend on the
   * random variables w_{t-1}^{dep} of the subproblem at stage t-1.
   *
   * @param stage The stage for which an initial state must be provided.
   *
   * @return An initial state for the optimization subproblem associated with
   *         the given stage. */

  Eigen::ArrayXd oneAdmissibleState( const double & stage ) override;

/*--------------------------------------------------------------------------*/

  /// return the size of the state vector
  /** This function returns the size of the state vector. It assumes
   * that the states of all stages have the same size. If the state
   * at a time t is given by (x_t, w_t^{dep}), then this function
   * should return the size of x_t plus the size of w_t^{dep}.
   *
   * @return The size of the state vector. */

  int getStateSize() const override;

/*--------------------------------------------------------------------------*/

  /// returns the simulator for the backward pass
  /** This function returns the simulator that is used during the
   * backward pass.
   *
   * @return The simulator associated with the backward pass. */

  std::shared_ptr< StOpt::SimulatorSDDPBase >
  getSimulatorBackward() const override {
   return simulator_backward;
  }

/*--------------------------------------------------------------------------*/

  /// returns the simulator for the forward pass
  /** This function returns the simulator that is used during the forward
   * pass.
   *
   * @return The simulator associated with the forward pass. */

  std::shared_ptr< StOpt::SimulatorSDDPBase >
  getSimulatorForward() const override {
   return simulator_forward;
  }

/*--------------------------------------------------------------------------*/

  /// set the SDDPSolver with which this SDDPOptimizer will be associated
  /** This method is used to set the (pointer to the) SDDPSolver with which
   * this SDDPOptimizer will be associated.
   *
   * @param solver A pointer to an SDDPSolver. */

  void set_solver( SDDPSolver * solver ) {
   sddp_solver = solver;
  }

/*--------------------------------------------------------------------------*/

  void set_scenarios( const ScenarioSet & scenario_set ) {
   simulator_backward->set_scenarios( scenario_set );
   simulator_forward->set_scenarios( scenario_set );

   if( simulator_backward->getNbSimul() == 0 )
    simulator_backward->set_number_simulations
     ( get_dflt_number_simulations_backward() );

   if( simulator_forward->getNbSimul() == 0 )
    simulator_forward->set_number_simulations( 1 );
  }

/*--------------------------------------------------------------------------*/

  int get_number_simulations_backward() const {
   return simulator_backward->getNbSimul();
  }

/*--------------------------------------------------------------------------*/

  int get_number_simulations_forward() const {
   return simulator_forward->getNbSimul();
  }

/*--------------------------------------------------------------------------*/

  int get_number_simulations( bool backward ) const {
   if( backward )
    return get_number_simulations_backward();
   return get_number_simulations_forward();
  }

/*--------------------------------------------------------------------------*/

  void set_number_simulations_backward( int number_simulations ) {
   simulator_backward->set_number_simulations( number_simulations );
  }

/*--------------------------------------------------------------------------*/

  void set_number_simulations_forward( int number_simulations ) {
   simulator_forward->set_number_simulations( number_simulations );
  }

/*--------------------------------------------------------------------------*/

  int get_dflt_number_simulations_backward() const {
   return simulator_backward->get_number_scenarios();
  }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/

  /// returns the current stage for a backward pass
  Index get_current_backward_stage() const {
   return date_next;
  }

/*--------------------------------------------------------------------------*/

  /// returns the current stage for a forward pass
  Index get_current_forward_stage() const {
   return date;
  }

/*--------------------------------------------------------------------------*/

  /// returns the index of the backward scenario associated with simulation_id
  /** This function returns the index of the backward scenario associated with
   * the given \p simulation_id.
   *
   * @param simulation_id The index of a simulation, which must be an integer
   *        between 0 and get_number_simulations_backward() - 1.
   *
   * @return The index of the backward scenario associated with the given \p
   *         simulation_id. */

  Index get_backward_scenario_index( Index simulation_id ) const {
   if( get_current_backward_stage() == 0 )
    return sddp_solver->get_first_stage_scenario_index();
   return simulator_backward->get_scenario_index( simulation_id );
  }

/*--------------------------------------------------------------------------*/

  /// returns the index of the forward scenario associated with simulation_id
  /** This function returns the index of the forward scenario associated with
   * the given \p simulation_id.
   *
   * @param simulation_id The index of a simulation, which must be an integer
   *        between 0 and get_number_simulations_forward() - 1.
   *
   * @return The index of the forward scenario associated with the given \p
   *         simulation_id. */

  Index get_forward_scenario_index( Index simulation_id ) const {
   if( get_current_forward_stage() == 0 )
    return sddp_solver->get_first_stage_scenario_index();
   return simulator_forward->get_scenario_index( simulation_id );
  }

/*--------------------------------------------------------------------------*/

  /// returns the index of the scenario associated with simulation_id
  /** This function returns the index of the scenario associated with the
   * given \p simulation_id. If backward is true, then it returns the backward
   * scenario index associated with \p simulation_id. Otherwise, it returns
   * the forward scenario index associated with \p simulation_id.
   *
   * @param simulation_id The index of a simulation, which must be an integer
   *        between 0 and get_number_simulations( backward ) - 1.
   *
   * @param backward Indicates whether the desired scenario index is
   *        associated with the backward simulator.
   *
   * @return The index of the forward scenario associated with the given \p
   *         simulation_id. */

  Index get_scenario_index( Index simulation_id , bool backward ) const {
   if( backward )
    return get_backward_scenario_index( simulation_id );
   return get_forward_scenario_index( simulation_id );
  }

/*--------------------------------------------------------------------------*/

  Index get_num_sub_blocks_per_stage() const {
   return static_cast< SDDPBlock * >( sddp_solver->f_Block )->
    get_num_sub_blocks_per_stage();
  }

/*--------------------------------------------------------------------------*/

  /// returns the pointer to a PolyhedralFunction
  /** This function returns a pointer to the i-th PolyhedralFunction of a
   * sub-Block of the given \p stage.
   *
   * @param stage A number between 0 and get_time_horizon() - 1.
   *
   * @param i If there are more than one PolyhedralFunction per sub-Block, this
   *        parameter informs the index of the desired PolyhedralFunction at
   *        the given \p stage.
   *
   * @param sub_block_index The index of the sub-Block, which must be an
   *        integer between 0 and get_num_sub_blocks_per_stage() - 1.
   *
   * @return A pointer to the i-th PolyhedralFunction of the given \p stage. */

  PolyhedralFunction * get_polyhedral_function
  ( Index stage , Index i = 0 , Index sub_block_index = 0 ) const {
   return static_cast< SDDPBlock * >( sddp_solver->f_Block )->
    get_polyhedral_function( stage , i , sub_block_index );
  }

/*--------------------------------------------------------------------------*/

  bool mesh_provided() const {
   return sddp_solver->mesh_provided();
  }

/*--------------------------------------------------------------------------*/

  int get_output_frequency() const {
   return sddp_solver->output_frequency;
  }

/*--------------------------------------------------------------------------*/

  double date;
  double date_next;
  std::shared_ptr< ScenarioSimulator > simulator_forward;
  std::shared_ptr< ScenarioSimulator > simulator_backward;

  /// the number of the current iteration
  mutable int current_iteration = 0;

  /// indicates whether the previous pass was backward
  mutable bool previous_pass_was_backward = false;

  SDDPSolver * sddp_solver;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

  /// checks if the linearization is correct
  void check_linearization( double objective_value , double alpha , double gy ,
                            const Eigen::ArrayXd & linearization ,
                            Index sub_block_index ) const;

 };   // end( class SDDPOptimizer )

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// Initial state
 /** The initial state at the beginning of the simulation. */
 std::vector< double > initial_state;

 /// Number of meshes in each direction
 /** This vector stores the mesh discretization in each direction. The i-th
  * component of this vector contains the number of meshes (number of steps)
  * at direction i. */
 std::vector< int > mesh_discretization;

 /// The cuts to be used at the last time instant
 /** This vector stores the cuts to be used at the last time instant. Each cut
  * is represented by a vector whose size is equal to m + 1, where m is the
  * number of state variables. The m first elements of a cut are the
  * coefficients for the state variables and the last element is the constant
  * of that cut. The vector #last_stage_cuts can store multiple cuts and its
  * size must be a multiple of m + 1. If it is non-empty and has size k * (m +
  * 1), then it contains k cuts and the i-th cut is given by the elements
  * between the indices i * ( m + 1 ) and ( i + 1 ) * ( m + 1 ) - 1. */
 std::vector< double > last_stage_cuts;

 /// Number of iterations performed by the method
 /** Number of iterations performed by the method at the last call of
  * compute(). */
 int number_iterations_performed;

 /// Accuracy achieved by the method
 /** Accuracy achieved by the method at the last call to compute(),
  * which is given by
  *
  * | backwardValue - forwardValue | / forwardValue
  *
  * where backwardValue is the value of the last backward pass and
  * forwardValueForConv is the value obtained in the forward pass
  * when checking for convergence.
  */
 double accuracy_achieved;

 /// It indicates the level of verbosity of the log
 int log_verbosity = 0;

 /// Index of the scenario to be considered at the first stage
 int first_stage_scenario_index;

 /// Name of the file in which regressors will be stored
 std::string regressors_filename;

 /// Name of the file in which the cuts will be stored
 std::string cuts_filename;

 /// Name of the file in which the visited states will be stored
 std::string visited_states_filename;

 /// Name of the file to which the future cost functions are output
 std::string f_output_filename;

 /// Prefix for the name of the file in which the State is saved
 std::string f_state_filename;

 /// Prefix for the name of the file to which the random cuts are output
 std::string f_random_cuts_filename;
 
 /// The path where to write the results
 std::string f_dir_out_pathname = "";

 /// The suffix to be added to an output filename every other iteration
 std::string f_filename_suffix = "";

 /// It indicates whether a suffix must be added to an output filename
 /** This variable indicates whether a suffix must be added to the name of an
  * output file (e.g., the file containing the future cost functions or the
  * file containing an SDDPSolverState). If this variable is false, then the
  * name of the file will be the given one (#f_output_filename for the future
  * cost functions and #f_state_filename for the SDDPSolverState). If this
  * variable is true, then #f_filename_suffix will be added to the name of the
  * file every other iteration in which the output is performed. */
 mutable bool f_add_suffix = false;

 /// Prefix of the names of the files for the logs of the sub-Solvers
 std::string f_sub_solver_filename_prefix = "";

 /// Name of the default BlockConfig file for the inner Blocks
 std::string f_inner_block_config_filename;

 /// Name of the default BlockSolverConfig file for the inner Blocks
 std::string f_inner_block_solver_config_filename;

 /// Default BlockConfig for the inner Blocks
 BlockConfig * f_inner_block_config = nullptr;

 /// Default BlockSolverConfig for the inner Blocks
 BlockSolverConfig * f_inner_block_solver_config = nullptr;

 /// Configuration to be passed to get_var_solution() of the inner Solver
 Configuration * f_get_var_solution_config = nullptr;

 /// Maximum number of iterations that the method should perform
 int maximum_number_iterations;

 /** Frequency in which the convergence check is performed. The
  * convergence is checked every "convergence_frequency" steps of the
  * method. See the comments about the intNStepConv parameter for
  * more details. */
 int convergence_frequency;

 /** Indicates whether the CPU time spent at each backward and
  * forward steps should be printed. */
 bool print_cpu_time;

 /** The number of simulations (number of forward passes called) when
  * we have to check the convergence by comparing the outcome given
  * by the forward pass and the one given by the backward pass. */
 int number_simulations_for_convergence;

 /// The frequency in which the future cost functions are output
 int output_frequency;

 /// Periodicity of eEverykIteration events
 int f_handle_events_every_k_iter;

 /** Relative accuracy for declaring a solution optimal. See the
  * comments about the dblAccuracy parameter for more details. */
 double accuracy;

 /// Status returned by compute()
 int status = kUnEval;

 /// Pointer to the SDDPOptimizer
 std::shared_ptr< SDDPOptimizer > sddp_optimizer;

#ifdef USE_MPI
 /// MPI communicator
 /** A function to set it could be implemented, but, for now, it is just the
  * default communicator. */
 boost::mpi::communicator mpi_communicator;
#endif

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// add cuts to the sub-Block at the given stage
 /** This function adds cuts to the sub-Block with index \p sub_block_index at
  * the given \p stage.
  *
  * @param cuts An Eigen::ArrayXXd containing the cuts to be added. It must be
  *        a matrix with as many columns as there are cuts to be added and the
  *        number of rows must be equal to one plus the number of Variable
  *        defined in the BendersBlock associated with stage \p stage.
  *
  * @param stage The stage at which cuts should be updated, which must be an
  *        integer between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the sub-Block, which must be an
  *        integer between 0 and get_num_sub_blocks_per_stage() - 1.
  *
  * @param range The indices of the cuts in \p cuts that should be added. */

 void add_cuts( const Eigen::ArrayXXd & cuts , SDDPBlock::Index stage ,
                SDDPBlock::Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 void set_state( const Eigen::ArrayXd & state , SDDPBlock::Index stage ,
                 SDDPBlock::Index  sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 void set_scenario( SDDPBlock::Index scenario_id , SDDPBlock::Index stage ,
                    SDDPBlock::Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification();

/*--------------------------------------------------------------------------*/

 /// solves the sub-Block associated with the given stage
 /** This function solves the sub-Block with index \p sub_block_index at the
  * given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the sub-Block, which must be an
  *        integer between 0 and get_num_sub_blocks_per_stage() - 1. */

 double solve( SDDPBlock::Index stage , SDDPBlock::Index sub_block_index );

/*--------------------------------------------------------------------------*/

 /// sets the parameters of the SDDPSolver to their default values
 /** This function sets the parameters of the SDDPSolver to their default
  * values. */

 void set_default_parameters() {

  // int

  maximum_number_iterations = get_dflt_int_par( intMaxIter );
  convergence_frequency = get_dflt_int_par( intNStepConv );
  print_cpu_time = get_dflt_int_par( intPrintTime );
  number_simulations_for_convergence =
   get_dflt_int_par( intNbSimulCheckForConv );
  output_frequency = get_dflt_int_par( intOutputFrequency );
  first_stage_scenario_index = get_dflt_int_par( intFirstStageScenarioId );

  // double

  accuracy = get_dflt_dbl_par( dblAccuracy );

  // string

  regressors_filename = get_dflt_str_par( strRegressorsFilename );
  cuts_filename = get_dflt_str_par( strCutsFilename );
  visited_states_filename = get_dflt_str_par( strVisitedStatesFilename );
  f_inner_block_config_filename = get_dflt_str_par( strInnerBC );
  f_inner_block_solver_config_filename = get_dflt_str_par( strInnerBSC );
  f_output_filename = get_dflt_str_par( strOutputFile );
  f_state_filename = get_dflt_str_par( strStateFile );
  f_random_cuts_filename = get_dflt_str_par( strRandomCutsFile );
  f_dir_out_pathname = get_dflt_str_par( strDirOUT );
  f_filename_suffix = get_dflt_str_par( strFilenameSuffix );
  f_sub_solver_filename_prefix = get_dflt_str_par( strSubSolverLogFilePrefix );
  f_handle_events_every_k_iter = Solver::get_dflt_int_par( intEverykIt );

  // vector of int

  mesh_discretization = get_dflt_vint_par( vintMeshDiscretization );

  // vector of double

  last_stage_cuts = get_dflt_vdbl_par( vdblLastStageCuts );
  initial_state = get_dflt_vdbl_par( vdblInitialState );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the value of the last backward pass
 double backward_value;

 /// the value obtained during the forward pass when checking for convergence
 double forward_value;

 /// the number of cuts already present at each stage when compute() is called
 std::vector< Index > number_initial_cuts;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

};   // end( class SDDPSolver )


/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SDDPSolverState --------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe the "internal state" of an SDDPSolver
/** Derived class from State to describe the "internal state" of an
 * SDDPSolver. An SDDPSolverState is formed by the data that characterizes the
 * cuts associated with each stage. For each stage t, the data of the
 * PolyhedralFunction associated with that stage that is relevant are the
 * following:
 *
 * - the cuts (a.k.a. rows) present in the PolyhedralFunction;
 *
 * - the global lower bound (if it is convex; upper bound if it is concave) on
 *   the value of the PolyhedralFunction;
 *
 * - the "verse" of the PolyhedralFunction, i.e., if it is convex or
 *   concave. */

class SDDPSolverState : public State {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

public:

/*------------- CONSTRUCTING AND DESTRUCTING SDDPSolverState ---------------*/

 /// constructor, doing everything
 /** Constructor of SDDPSolverState. If a pointer to an SDDPSolver is
  * provided, then it immediately copies its "internal state". If nullptr is
  * passed (as by default), then an "empty" SDDPSolverState is constructed
  * that can only be filled by calling deserialize(). */

 SDDPSolverState( const SDDPSolver * solver = nullptr );

/*--------------------------------------------------------------------------*/

 /// de-serialize an SDDPSolverState out of a netCDF::NcGroup
 /** De-serialize an SDDPSolverState out of the given netCDF::NcGroup; see
  * SDDPSolverState::serialize() for a description of the format.
  *
  * @param group The netCDF::NcGroup out of which this SDDPSolverState will be
  *        de-serialized. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 ///< destructor

 virtual ~SDDPSolverState() { }

/*---------- METHODS DESCRIBING THE BEHAVIOR OF A SDDPSolverState ----------*/

 /// serialize an SDDPSolverState into a netCDF::NcGroup
 /** This method serializes this SDDPSolverState into the provided
  * netCDF::NcGroup, so that it can later be read back by deserialize().
  *
  * After this SDDPSolverState is serialized, \p group will have the dimension
  * "TimeHorizon", containing the time horizon, and, for each t in {0, ...,
  * TimeHorizon - 1}, the following data of the PolyhedralFunction associated
  * with stage t:
  *
  * - The dimension "PolyFunction_sign_t" (actually a bool), which contains
  *   the "verse" of the PolyhedralFunction, i.e., true for a convex
  *   max-function and false for a concave min-function (encoded in the
  *   obvious way, i.e., zero for false, nonzero for true). This dimension is
  *   optional: if it is not provided, true is assumed.
  *
  * - The dimension "PolyFunction_NumRow_t", containing the number of rows of
  *   the A matrix. This dimension is optional: if it is not provided, then 0
  *   (no rows) is assumed.
  *
  * - The dimension "PolyFunction_NumVar_t", containing the number of columns
  *   of the A matrix, i.e., the number of active variables.
  *
  * - The variable "PolyFunction_A_t", of type netCDF::NcDouble() and indexed
  *   over both the dimensions "PolyFunction_NumRow_t" and
  *   "PolyFunction_NumVar_t" (in this order); it contains the (row-major)
  *   representation of the matrix A. This variable is only optional if
  *   "PolyFunction_NumRow_t" == 0.
  *
  * - The variable "PolyFunction_b_t", of type netCDF::NcDouble() and indexed
  *   over the dimension "PolyFunction_NumRow_t", which contains the vector
  *   b. This variable is only optional if "PolyFunction_NumRow_t" == 0.
  *
  * - The scalar variable "PolyFunction_lb_t", of type netCDF::NcDouble() and
  *   not indexed over any dimension, which contains the global lower (if
  *   PolyFunction_sign_t == true, upper otherwise) bound on the value of the
  *   PolyhedralFunction over all the space. This variable is optional: if it
  *   is not provided, it means that no finite lower (upper) bound exist,
  *   i.e., the lower (upper) bound is -(+)
  *   Inf< PolyhedralFunction::FunctionValue >(). */

 void serialize( netCDF::NcGroup & group ) const override;

/*-------------------------------- FRIENDS ---------------------------------*/

 friend class SDDPSolver;  // make SDDPSolver friend

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "SDDPSolverState [" << this << "] with";
  if( v_b.empty() )
   output << " no cuts.";
  else {
   output << std::endl;
   for( Index t = 0 ; t < v_b.size() ; ++t )
    output << v_b[ t ].size() << " cuts for stage " << t << std::endl;
  }
 }

/*--------------------------- PROTECTED FIELDS -----------------------------*/

 std::vector< bool > v_is_convex;
 ///< true if the PolyhedralFunction is a convex one

 std::vector< Index > v_num_var;
 ///< the number of variables of each PolyhedralFunction

 std::vector< PolyhedralFunction::MultiVector > v_A;
 ///< the A matrix of each PolyhedralFunction

 std::vector< PolyhedralFunction::RealVector > v_b;
 ///< the b vector of each PolyhedralFunction

 std::vector< PolyhedralFunction::FunctionValue > v_bound;
 ///< the global (lower or upper) bound of each PolyhedralFunction

/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/

private:

/*--------------------------- PRIVATE METHODS ------------------------------*/

 static void serialize
 ( netCDF::NcGroup & group , Index t , Index num_var , bool is_convex ,
   PolyhedralFunction::FunctionValue bound ,
   const PolyhedralFunction::MultiVector & A ,
   const PolyhedralFunction::RealVector & b );

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SMSpp_insert_in_factory_h;

};  // end( class( SDDPSolverState ) )

/** @} end( group( SDDPSolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
