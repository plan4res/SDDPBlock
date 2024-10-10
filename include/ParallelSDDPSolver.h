/*--------------------------------------------------------------------------*/
/*--------------------- File ParallelSDDPSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ParallelSDDPSolver class, which derives from SDDPSolver
 * and implements a parallel version of SDDPSolver.
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

#ifndef __ParallelSDDPSolver
#define __ParallelSDDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SDDPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ParallelSDDPSolver_CLASSES Classes in ParallelSDDPSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ParallelSDDPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an SDDPSolver where the subproblems can be solved in parallel
/** The ParallelSDDPSolver derives from SDDPSolver to implement the
 * computation of the subproblems in parallel.
 *
 * Each inner iteration of the SDDP method, in both forward and backward
 * passes, is associated with a particular stage and considers a set of
 * scenarios. The subproblem associated with that stage must be solved for
 * each scenario in that set of scenarios. (Each subproblem also depends on an
 * initial state, which will be omitted here for simplicity.) In the simplest
 * case, the SDDPBlock contains a single sub-Block for each stage. In this
 * case, this single sub-Block must be solved, in sequence, for each
 * scenario. Solving a sub-Block for a particular scenario requires changing
 * the data of the sub-Block to reflect that scenario. Thus, in this case,
 * this *sequential* procedure can be summarized as follows:
 *
 *   for each scenario S in the given set of scenarios, do
 *     update the sub-Block according to S
 *     solve the sub-Block
 *
 * The SDDPBlock can also contain more than one sub-Block for each stage. This
 * opens the possibility of solving the subproblems in parallel, each in a
 * different thread and associated with a different sub-Block. Suppose that
 * the SDDPBlock contains B sub-Blocks at stage t. All of these sub-Blocks are
 * identical, except perhaps for the realization of their stochastic
 * data. Considering that at the beginning of the inner iteration all
 * sub-Blocks are considered "available" to be "used", we could have the
 * following procedure to solve the subproblems for each scenario:
 *
 *   for each scenario S in the given set of scenarios, do
 *     select and lock an available sub-Block, let us say sub-Block B_i
 *     update B_i according to S
 *     solve B_i
 *     unlock B_i
 *
 * Supposing that a process running ParallelSDDPSolver::compute() has M
 * threads available, it would be possible to allocate min( B , N , M )
 * subproblems to the available threads, where N is the number of scenarios
 * to be considered at an inner iteration of the SDDP method. */

class ParallelSDDPSolver : public SDDPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend class ParallelSDDPOptimizer;

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ CONSTRUCTING AND DESTRUCTING ParallelSDDPSolver -------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ParallelSDDPSolver
 *  @{ */

 /// constructor
 ParallelSDDPSolver() {
  sddp_optimizer = std::make_shared< ParallelSDDPOptimizer >( this );
 }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~ParallelSDDPSolver() {}

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED CLASSES -----------------------------*/
/*--------------------------------------------------------------------------*/

 class ParallelSDDPOptimizer : public SDDPSolver::SDDPOptimizer {

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 public:

  /// constructor taking a pointer to a ParallelSDDPSolver
  /** Constructs a ParallelSDDPOptimizer associated with the given
   * ParallelSDDPSolver.
   *
   * @param solver A pointer to a ParallelSDDPSolver. This parameter is
   *        optional and its default value is nullptr. */

  ParallelSDDPOptimizer( ParallelSDDPSolver * solver = nullptr ) :
   SDDPSolver::SDDPOptimizer( solver ) { }

/*--------------------------------------------------------------------------*/

  void reset() override {
   SDDPSolver::SDDPOptimizer::reset();
   new_stage = true;
   cuts_synchronized = false;
  }

/*--------------------------------------------------------------------------*/

  Eigen::ArrayXd oneStepBackward
  ( const StOpt::SDDPCutOptBase & p_linCut,
    const std::tuple< std::shared_ptr< Eigen::ArrayXd >, int, int > & p_aState,
    const Eigen::ArrayXd & p_particle, const int & p_isample) const override;

/*--------------------------------------------------------------------------*/

  double oneStepForward
  ( const Eigen::ArrayXd &p_aParticle, Eigen::ArrayXd &p_state,
    Eigen::ArrayXd &p_stateToStore,
    const StOpt::SDDPCutOptBase &p_linCut,
    const int &p_isimu ) const override;

/*--------------------------------------------------------------------------*/

  void updateDates( const double & date, const double & date_next ) override {
   SDDPSolver::SDDPOptimizer::updateDates( date , date_next );
   new_stage = true;
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/

  /// returns true if, and only if, the sub-Block with the given index is locked
  bool is_locked( Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

  /// lock the sub-Block with the given index
  void lock_sub_block( Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

  /// returns the index of a sub-Block at the given stage with a scenario
  /** This function returns the index of a sub-Block of SDDPBlock at the given
   * \p stage whose data is set to the scenario whose index is \p
   * scenario_index. If no sub-Block at the given \p stage is set to the
   * scenario specified, Inf< Index >() is returned.
   *
   * @param stage A stage between 0 and time_horizon - 1.
   *
   * @param scenario_index The index of a scenario.
   *
   * @return The index of a sub-Block of SDDPBlock at the given \p stage whose
   *         data is set to the scenario whose index is \p scenario_index. If
   *         no such a sub-Block exists, Inf< Index >() is returned. */

  Index get_sub_block_with_scenario( Index stage ,
                                     Index scenario_index ) const;

/*--------------------------------------------------------------------------*/

  /// locks a sub-Block and returns its index
  /** This function locks a sub-Block at the given \p stage in order for the
   * subproblem associated with the given simulation id to be solved. The
   * parameter \p backward indicates whether the current pass is a backward
   * one.
   *
   * @param stage The stage at which a sub-Block must be locked.
   *
   * @param simulation_id The id of the simulation for which the subproblem
   *        will be solved.
   *
   * @param backward It indicates whether the current pass is a backward one.
   *
   * @return The index of the sub-Block that was locked. */

  Index lock( Index stage , Index simulation_id , bool backward ) const;

/*--------------------------------------------------------------------------*/

  /// unlocks the sub-Block with the given index
  /** This function unlocks the sub-Block at the given \p stage whose index is
   * \p sub_block_index. A sub-Block is always locked in order for it to be
   * solved for a particular scenario. Thus, \p scenario_index is the index of
   * the scenario for which the sub-Block was solved. Moreover, \p
   * scenario_was_set indicates whether the scenario was set to the sub-Block
   * (which is usually the case). */

  void unlock( Index stage , Index sub_block_index , Index scenario_index ,
               bool scenario_was_set ) const;

/*--------------------------------------------------------------------------*/

  /// indicates whether the given scenario is present in the given sub-Block
  /** This function returns true if, and only if, the data of the sub-Block
   * with index \p sub_block_index at the given \p stage is currently set to
   * the scenario whose index is \p scenario_index. */

  bool is_scenario_set( Index stage , Index sub_block_index ,
                        Index scenario_index ) const;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

  Index find_available_sub_block( Index simulation_id , Index stage ) const;

/*--------------------------------------------------------------------------*/

  void synchronize_cuts() const;

/*--------------------------------------------------------------------------*/

  void prepare_new_stage( Index stage , bool backward ) const;

/*--------------------------------------------------------------------------*/

  /// This is the waiting time before trying to acquire the lock again
  double waiting_time = 0.0;

  /// It indicates whether the current subproblem is the first one of his stage
  mutable bool new_stage = true;

  /// This vector indicates whether each sub-Block is locked
  mutable std::vector< bool > locked;

  /// It indicates whether the cuts of all sub-Blocks are synchronized
  mutable bool cuts_synchronized = false;

  /// It indicates the index of the sub-Block that is preferred by a simulation
  /** For each simulation i, preferred_block[ i ] indicates the index of the
   * sub-Block that has higher priority to be associated with simulation i. If
   * preferred_block[ i ] is Inf, then simulation i has no preferred
   * sub-Block. */
  mutable std::vector< Index > preferred_block;

  /// It indicates whether each sub-Block is not "preferred" by any simulation
  mutable std::vector< bool > non_reserved_block;

  /// Indices of sub-Blocks which are "preferred" by some simulation
  mutable std::vector< Index > reserved_blocks;

  /// For each stage, it gives the index of the last sub-Block solved
  mutable std::vector< Index > last_sub_block_solved;

  /// Stores the indices of the scenarios currently set in each sub-Block
  /* For each stage t, scenario_currently_set is a vector with size N_t (the
   * number of sub-Block at stage t). scenario_currently_set[ t ][ i ] is the
   * index of the scenario that is currently set in the i-th Block (if no
   * scenario is set, it stores Inf).
   *
   * Note: this control only works if the set of scenarios does not change in
   * any way. If it changes, we must deal with a Modification, which we
   * currently do not do. We can make it empty in add_Modification() if any
   * Modification arrives. */
  mutable std::vector< std::vector< Index > > scenario_currently_set;


/*--------------------------------------------------------------------------*/

 };   // end( class ParallelSDDPOptimizer )

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

};   // end( class ParallelSDDPSolver )

/** @} end( group( ParallelSDDPSolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ParallelSDDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ParallelSDDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
