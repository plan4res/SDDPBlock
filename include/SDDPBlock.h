/*--------------------------------------------------------------------------*/
/*------------------------- File SDDPBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of SDDPBlock, a class for representing a multistage stochastic
 * programming problem specifically designed to be solved by an SDDP solver.
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

#ifndef __SDDPBlock
 #define __SDDPBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "Objective.h"
#include "PolyhedralFunction.h"
#include "ScenarioSimulator.h"
#include "ScenarioSet.h"
#include "StochasticBlock.h"
#include "StOpt/sddp/SimulatorSDDPBase.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SDDPBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// SDDPBlock, representing a multistage stochastic programming problem
/** The SDDPBlock is a class that derives from Block and represents a
 * multistage stochastic programming problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{X}_0} f_0(x_0) +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_1 \in \mathcal{X}_1} f_1(x_1) +
 *   \mathbb{E} \left \lbrack \dots +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{X}_{T-1}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack,
 * \f]
 *
 * where T is called the time horizon, \f$\mathcal{X}_t \equiv
 * \mathcal{X}_t(x_{t-1}, \xi_t) \subseteq \mathbb{R}^{n_t}\f$ for each
 * \f$t \in \{0, \dots, T-1\}\f$, and \f$ \xi = \{ \xi_t \}_{t \in \{1, \dots,
 * T-1\}} \f$ is a stochastic process. Notice that \f$ x_{-1} \f$ and \f$
 * \xi_0 \f$ are deterministic. For each \f$ t \in \{0, \dots, T-1\}\f$, we
 * call
 *
 * \f[
 *   \min_{x_t \in \mathcal{X}_t} f_t(x_t) +
 *   \mathcal{V}_{t+1}(x_t)
 * \f]
 *
 * the problem associated with stage \f$ t \f$, where
 *
 * \f[
 *   \mathcal{V}_{t+1}(x_t) =
 *    \mathbb{E}
 *      \left\lbrack
 *        V_{t+1}(x_t, \xi_{t+1})
 *      \right\rbrack
 * \f]
 *
 * is the (expected value) cost-to-go function (also called value function,
 * future value function, future cost function), with \f$ \mathcal{V}_{T}
 * \equiv 0 \f$ and
 *
 * \f[
 *
 *    V_{t}(x_{t-1}, \xi_{t}) =
 *    \min_{x_t \in \mathcal{X}_t} f_t(x_t) +
 *    \mathcal{V}_{t+1}(x_t)
 * \f]
 *
 * with given \f$ x_{-1} \f$ and (deterministic) \f$ \xi_0\f$. We consider an
 * approximation to the problem associated with stage \f$ t \in \{0, \dots,
 * T-1\} \f$ as the problem
 *
 * \f[
 *    \min_{x_t \in \mathcal{X}_t} f_t(x_t) +
 *    \mathcal{P}_{t+1}(x_t)
 *    \qquad (1)
 * \f]
 *
 * where \f$ \mathcal{P}_{t+1}(x_t) \f$ is a polyhedral function, i.e., it is
 * a function of the form
 *
 * \f[
 *    \mathcal{P}_{t+1}(x_t) = \max_{i \in \{1,\dots,k_t\}}
 *                                     \{ d_{t,i}^{\top}x_t + e_{t,i} \}
 * \f]
 *
 * with \f$ d_{t,i} \in \mathbb{R}^{n_t} \f$ and \f$ e_{t,i} \in \mathbb{R}
 * \f$ for each \f$ i \in \{1,\dots,k_t\} \f$.
 *
 * An SDDPBlock is then characterized by the following:
 *
 * - It has a time horizon T.
 *
 * - It has T sub-Blocks, each one being a StochasticBlock. The t-th sub-Block
 *   represents an approximation to the problem associated with stage t as
 *   defined in (1). See the note below for the case in which it may have more
 *   than T sub-Blocks.
 *
 * - It has pointers to "T - 1" PolyhedralFunction. The t-th
 *   PolyhedralFunction represents the function \f$ \mathcal{P}_{t+1} \f$ in
 *   (1) and, therefore, must be defined in the t-th sub-Block of this
 *   SDDPBlock or in any of the sub-Blocks of that sub-Block, recursively.
 *
 * - It has a set of scenarios \f$\mathcal{S}\f$. Each scenario in
 *   \f$\mathcal{S}\f$ is represented by a vector of double and spans all the
 *   time horizon T. Each vector is divided into T parts, each one being
 *   associated with a stage of the multistage problem. Let \f$S\f$ denote a
 *   vector representing a scenario in \f$\mathcal{S}\f$. Then, \f$S\f$ is
 *   defined as
 *
 *   \f[
 *     S = ( S_0 , \dots, S_{T-1} )
 *   \f]
 *
 *   where \f$S_t\f$ is a sub-vector of \f$S\f$ with size \f$s_t\f$, for each
 *   \f$t \in \{ 0, \dots, T-1 \}\f$, and is associated with the sub-problem
 *   at stage \f$t\f$, i.e., it provides data for the \f$t\f$-th sub-Block of
 *   this SDDPBlock. We say that \f$S_t\f$ represents the \f$t\f$-th
 *   sub-scenario of the scenario represented by \f$S\f$.
 *
 *   We assume that the sub-scenarios are organized in such a way that related
 *   random data appear in contiguous areas of the sub-scenario. For instance,
 *   suppose that the random data is associated with demand, inflow, and wind
 *   power. In this case, the data related to demand should be a contiguous
 *   sub-vector \f$D_t\f$ of the sub-scenario associated with stage \f$t\f$,
 *   as well as that related to inflow (\f$F_t\f$) and wind power
 *   (\f$W_t\f$). In this example, the sub-scenario \f$S_t\f$ could be
 *   organized as
 *
 *   \f[
 *   S_t = ( D_t , F_t , W_t ).
 *   \f]
 *
 *   We say that this sub-scenario has three groups of related random
 *   data. The order in which the groups of related random data appear in
 *   \f$S_t\f$ is not relevant. We could have, for instance,
 *
 *   \f[
 *   S_t = ( W_t , D_t , F_t ).
 *   \f]
 *
 *   But the sub-scenario associated with stage \f$t\f$ must respect the same
 *   order for each scenario in \f$\mathcal{S}\f$.
 *
 * In the simplest case, an SDDPBlock has \f$ T \f$ sub-Blocks, the \f$t\f$-th
 * one being an StochasticBlock associated with stage \f$t\f$. However, it may
 * be interesting to have more than one sub-Block associated with each stage
 * in some circumstances. The solution process implemented by SDDPSolver, for
 * instance, involves the solution of a number of subproblems, each of them
 * associated with a stage \f$ t \f$, a particular sub-scenario \f$ S_t^i \f$,
 * and some initial state. There are at least two clear situations under which
 * the presence of multiple sub-Blocks for each stage can be beneficial to
 * SDDPSolver.
 *
 * -# To change the sub-Block as little as possible.
 *
 *    Whenever a subproblem associated with a stage \f$ t \f$ must be solved,
 *    the data of the Block associated with that subproblem must be updated
 *    according to some sub-scenario \f$ S_t^i \f$ and some initial state. In
 *    the case in which the SDDPBlock has \f$ T \f$ sub-Blocks, this means
 *    that its \f$t\f$-th sub-Block must be updated every time a particular
 *    scenario and state is considered. In order to allow reoptimization,
 *    SMS++ is designed to deal with changes in the data of a Block by means
 *    of its Modification mechanism. However, one would expect, in general,
 *    that the less a Block is modified, the faster it can be reoptimized. In
 *    the ideal case, an SDDPBlock would have as much sub-Blocks for each
 *    stage as there are scenarios. In this case, each sub-Block would be
 *    associated with a particular scenario, and the data of each of these
 *    sub-Blocks that depend on the scenarios would be updated only once, in
 *    the beginning. Of course, the sub-Block must still be modified every
 *    time before it is solved, because it also depends on the initial
 *    state. But in general, most of the data in a Block that needs to be
 *    updated is dependent on the scenarios. Therefore, having one sub-Block
 *    associated with each scenario would imply that only the initial state of
 *    the sub-Block must be updated, and the reoptimization could be expected
 *    to be faster (not to mention the process of modifying the scenario of a
 *    sub-Block that would also be avoided by itself).
 *
 * -# To allow parallelization in a shared-memory multiprocessing system.
 *
 *    At each iteration, and for each stage, SDDPSolver must solve a set of
 *    subproblems, each one associated with some scenario and some initial
 *    state. If SDDPBlock has only one sub-Block for each stage, then the
 *    process of solving all those subproblems is inevitably sequential, as
 *    solving a particular subproblem requires changing the data of that
 *    sub-Block. The presence of multiple sub-Blocks per stage, however, makes
 *    it possible to parallelize this procedure. Suppose, for instance, that
 *    at every iteration, for each stage, SDDPSolver must solve N subproblems,
 *    each one associated with some scenario and some initial state. Suppose
 *    also that SDDPBlock has B sub-Blocks per stage and a process running
 *    ParallelSDDPSolver::compute() has M threads available. In this case, it
 *    would be possible to allocate min( B , N , M ) subproblems to the
 *    available threads.
 *
 * The decision about the number of Blocks per stage must be well thought out,
 * as it depends, in particular, on the memory resources available. Besides
 * the memory required to store multiple sub-Blocks per stage, one has also to
 * take into account the memory required by the Solver attached to the inner
 * Block of each sub-Block of this SDDPBlock.
 */

class SDDPBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPBlock ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPBlock
 *  @{ */

 /// constructor
 /** Constructs an SDDPBlock with the given \p father Block. The input
  * parameter has a default value (\c nullptr), so that this can be used as
  * the void constructor.
  *
  * @param father A pointer to the father Block of this SDDPBlock.
  */
 SDDPBlock( Block * father = nullptr ) : Block( father ) { }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPBlock() {
  for( auto & block : v_Block )
   delete block;
  v_Block.clear();
 }

/*--------------------------------------------------------------------------*/
 /// loads SDDPBlock out of an istream - not implemented yet

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "SDDPBlock::load: method not implemented yet." ) );
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize an SDDPBlock out of netCDF::NcGroup
 /** The method takes a netCDF::NcGroup supposedly containing all the
  * information required to de-serialize the SDDPBlock. Besides the mandatory
  * "type" attribute of any :Block, the group must contain the following:
  *
  * - The "TimeHorizon" dimension, containing the time horizon.
  *
  * - The description of the sub-Blocks of the SDDPBlock. This is given by the
  *   sub-groups "StochasticBlock" and "StochasticBlock_t", for each t in {0,
  *   ..., "TimeHorizon - 1"}. These sub-groups are optional, but they cannot
  *   be all absent. If "StochasticBlock_t" is not provided for some t in {0,
  *   ..., "TimeHorizon - 1"}, then the "StochasticBlock" group must be
  *   provided and contain a complete description of the t-th sub-Block of
  *   this SDDPBlock. If the "StochasticBlock" group is not provided, then
  *   "StochasticBlock_t" must be provided for each t in {0, ..., "TimeHorizon
  *   - 1"} and contain a complete description of the t-th sub-Block of this
  *   SDDPBlock.
  *
  *   If "StochasticBlock_t" is provided but the description of its inner
  *   Block is not provided, then the "StochasticBlock" group must be provided
  *   and contain the description of an inner Block of a StochasticBlock. In
  *   this case, the description of the inner Block provided in the
  *   "StochasticBlock" group will be used to construct the inner Block of the
  *   StochasticBlock described by the "StochasticBlock_t" group.
  *
  *   If "StochasticBlock_t" is provided but the description of its vector of
  *   DataMapping is not provided, then if the "StochasticBlock" group is
  *   provided and contains a description of a vector of DataMapping, then it
  *   is used to construct the vector of DataMapping of the StochasticBlock
  *   described by the "StochasticBlock_t" group.
  *
  * - The "NumPolyhedralFunctionsPerSubBlock" dimension, containing the number
  *   of PolyhedralFunction that are present in each sub-Block. This dimension
  *   is optional. If it is not provided, then we assume that there is a
  *   single PolyhedralFunction in each sub-Block.
  *
  * - The AbstractPath group containing the description of a vector of
  *   AbstractPath as described in the AbstractPath class. The number of
  *   AbstractPath must be equal to either "NumPolyhedralFunctionsPerSubBlock"
  *   or "NumPolyhedralFunctionsPerSubBlock * TimeHorizon". If the number of
  *   AbstractPath is "NumPolyhedralFunctionsPerSubBlock" then the i-th
  *   PolyhedralFunction of each sub-Block is given by the i-th
  *   AbstractPath. If the number of AbstractPath is
  *   "NumPolyhedralFunctionsPerSubBlock * TimeHorizon" then the i-th
  *   PolyhedralFunction of a sub-Block at stage t is given by the
  *   AbstractPath at position "i + t * NumPolyhedralFunctionsPerSubBlock" for
  *   each t in {0, ..., TimeHorizon-1} and i in {0, ...,
  *   NumPolyhedralFunctionsPerSubBlock-1}. An AbstractPath associated with a
  *   stage t is taken with respect to the inner Block of a sub-Block for
  *   stage t of this SDDPBlock.
  *
  * - The "NumSubBlocksPerStage" dimension, containing the number of
  *   sub-Blocks that must be constructed for each stage. This dimension is
  *   optional. If it is not provided, then we assume that there is a single
  *   sub-Block for each stage.
  *
  * - The description of a ScenarioSet, as specified in the comments to
  *   ScenarioSet::deserialize().
  *
  * - The "InitialState" variable, a one-dimensional array of type
  *   netCDF::NcDouble, containing an initial state for the first stage
  *   problem.
  *
  * - The "StateSize" variable, of type netCDF::Uint and being either a scalar
  *   or a one-dimensional array indexed over "TimeHorizon" dimension,
  *   specifying the sizes of the states at each stage. If this variable is a
  *   scalar, then all states are assumed to have the same size given by
  *   "StateSize". If it is an array then, for each t in {0, ...,
  *   TimeHorizon-1}, StateSize[t] contains the size of the final state at
  *   stage t. The state being a vector, its size is the dimension of the
  *   space in which it lies.
  *
  * - The "AdmissibleState" variable, of type netCDF::NcDouble, containing an
  *   admissible state for each stage. An admissible state for a stage is a
  *   feasible final state for the problem at that stage. If "StateSize" is
  *   scalar and "AdmissibleState" has dimension "StateSize", then all stages
  *   are assumed to have the same admissible state given by
  *   "AdmissibleState". Otherwise, "AdmissibleState" contains the
  *   concatenation of the states for all stages as follows.
  *
  *   - If "StateSize" is scalar then, for each t in {0, ..., TimeHorizon -
  *     1}, an admissible state for stage t is given by
  *
  *       (AdmissibleState[s_t], ..., AdmissibleState[s_t + StateSize - 1]),
  *
  *     where s_t = t * StateSize. In this case, "AdmissibleState" must have
  *     size "TimeHorizon * StateSize".
  *
  *   - If "StateSize" is a one-dimensional array then, for each t in {0, ...,
  *     "TimeHorizon - 1"}, an admissible state for stage t is given by
  *
  *       (AdmissibleState[s_t], ..., AdmissibleState[s_t + StateSize[t] - 1]),
  *
  *     where s_t = \f$ \sum_{i=0}^{t-1} \f$ StateSize[i]. In this case,
  *     "AdmissibleState" must have size equal to
  *
  *     \f[
  *        \sum_{i=0}^{\text{TimeHorizon} - 1} \text{StateSize}[i].
  *     \f]
  *
  * @param group A netCDF::NcGroup holding the data describing this SDDPBlock.
  */

 void deserialize( const netCDF::NcGroup & group ) override {
  // TimeHorizon

  Index time_horizon;
  ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                        time_horizon , false );

  // NumSubBlocksPerStage

  Index NumSubBlocksPerStage;
  if( ::SMSpp_di_unipi_it::deserialize_dim
      ( group , "NumSubBlocksPerStage" , NumSubBlocksPerStage ) ) {
   num_sub_blocks_per_stage = NumSubBlocksPerStage;
  }

  // StochasticBlock

  v_Block.reserve( time_horizon * num_sub_blocks_per_stage );

  for( Index i = 0 ; i < time_horizon ; ++i )
   for( Index j = 0 ; j < num_sub_blocks_per_stage ; ++j )
    v_Block.push_back( deserialize_sub_Block( group , i ) );

  // PolyhedralFunctions

  auto path_group = group.getGroup( "AbstractPath" );

  auto paths = AbstractPath::vector_deserialize( path_group );

  if( ! ::SMSpp_di_unipi_it::deserialize_dim
      ( group , "NumPolyhedralFunctionsPerSubBlock" ,
        num_polyhedral_per_sub_block ) ) {
   num_polyhedral_per_sub_block = 1;
  }

  if( paths.size() != num_polyhedral_per_sub_block * time_horizon &&
      ! ( paths.size() == num_polyhedral_per_sub_block && time_horizon > 1 ) ) {
   if( num_polyhedral_per_sub_block == 1 )
    throw( std::invalid_argument(
     "SDDPBlock::deserialize: The number of AbstractPath to "
     "PolyhedralFunction must be either equal to 1 or equal to "
     "the time horizon." ) );
   else
    throw( std::invalid_argument(
     "SDDPBlock::deserialize: The number of AbstractPath to "
     "PolyhedralFunction must be either equal to K or equal to K "
     "times the time horizon, where K is the number of "
     "PolyhedralFunction per sub-Block." ) );
  }

  v_polyhedral_functions.clear();
  v_polyhedral_functions.reserve
   ( num_sub_blocks_per_stage * num_polyhedral_per_sub_block * time_horizon );

  for( Index t = 0 ; t < time_horizon ; ++t ) {
   for( Index j = 0 ; j < num_sub_blocks_per_stage ; ++j ) {
    auto reference_block = get_sub_Block( t , j )->get_nested_Block( 0 );
    assert( reference_block );
    for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
     Index path_index = num_polyhedral_per_sub_block * t + i;
     if( paths.size() == num_polyhedral_per_sub_block )
      path_index = i;
     auto polyhedral_function = dynamic_cast< PolyhedralFunction * >
      ( paths[ path_index ].get_element< Function >( reference_block ) );
     if( ! polyhedral_function )
      throw( std::invalid_argument(
       "SDDPBlock::deserialize: PolyhedralFunction for stage "
       + std::to_string( t ) + " was not found." ) );
     v_polyhedral_functions.push_back( polyhedral_function );
    }
   }
  }

  // Scenarios

  scenario_set.deserialize( group );

  // Initial state

  ::SMSpp_di_unipi_it::deserialize( group , "InitialState" ,
                                    initial_state , false );

  // StateSize

  std::vector< Index > state_size;

  ::SMSpp_di_unipi_it::deserialize( group , "StateSize" , { time_horizon } ,
                                    state_size , false , true );

  bool state_size_is_scalar = ( state_size.size() == 1 );
  if( state_size.size() == 1 )
   state_size.resize( time_horizon , state_size[ 0 ] );
  else if( state_size.size() != time_horizon )
   throw( std::logic_error( "SDDPBlock::deserialize: 'StateSize' must be "
                            "either a scalar or an array with size "
                            "'TimeHorizon'." ) );

  // AdmissibleState

  ::SMSpp_di_unipi_it::deserialize( group , "AdmissibleState" ,
                                    admissible_states , false );

  if( state_size_is_scalar ) {
   if( admissible_states.size() != state_size[ 0 ] &&
       admissible_states.size() != time_horizon * state_size[ 0 ] )
    throw( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                             "array has an invalid size." ) );

   if( admissible_states.size() != time_horizon * state_size[ 0 ] ) {
    std::vector< double > state = admissible_states;
    admissible_states.reserve( time_horizon * state_size[ 0 ] );
    for( Index t = 1 ; t < time_horizon ; ++t )
     admissible_states.insert( admissible_states.cend() ,
                               state.cbegin() , state.cend() );
   }
  }
  else if( admissible_states.size() !=
           std::accumulate( state_size.begin() , state_size.end() ,
                            decltype( state_size )::value_type( 0 ) ) ) {
   throw( std::logic_error( "SDDPBlock::deserialize: 'AdmissibleState' "
                            "array has an invalid size." ) );
  }

  // Construct the vector admissible_state_begin

  admissible_state_begin.resize( time_horizon );
  if( time_horizon > 0 )
   admissible_state_begin.front() = 0;
  for( Index t = 1 ; t < time_horizon ; ++t )
   admissible_state_begin[ t ] =
    admissible_state_begin[ t - 1 ] + state_size[ t - 1 ];

  Block::deserialize( group );

 }

/*--------------------------------------------------------------------------*/

 /// deserialize the random cuts
 /** This function deserializes the random cuts out of the file whose path is
  * given by \p filename. If the path to the file is empty, then no operation
  * is performed. The file must have the following netCDF format:
  *
  * - The netCDF dimension "TimeHorizon" containing the number of stages.
  *
  * - The netCDF dimension "NumberScenarios" containing the number of
  *   scenarios.
  *
  * - The netCDF group "PolyhedralFunction_t_s", for each t in {0, ...,
  *   TimeHorizon - 1} and s in {0, ..., NumberScenarios - 1}, containing the
  *   serialization of the PolyhedralFunction representing the random cuts
  *   associated with stage t and scenario s.
  *
  * @param filename The path to the file containing the netCDF description of
  *        the random cuts. */

 void deserialize_random_cuts( const std::string & filename );

/*--------------------------------------------------------------------------*/

 /// sets the number of sub-Blocks for each stage
 /** This function sets the number of sub-Blocks that must be constructed at
  * each stage. If this function is invoked after the sub-Blocks of this
  * SDDPBlocks have been constructed, it has no effect. In particular, it has
  * no effect if it is invoked after deserialize() is invoked.
  *
  * @param n The number of sub-Blocks that must be constructed at each stage.
  */
 void set_num_sub_blocks_per_stage( Index n ) {
  if( v_Block.empty() )
   num_sub_blocks_per_stage = n;
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR Saving THE DATA OF THE SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Saving the data of the SDDPBlock
 *  @{ */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// serialize an SDDPBlock into a netCDF::NcGroup
 /** Serialize an SDDPBlock into a netCDF::NcGroup with the format
  * explained in the comments of the deserialize() function.
  *
  * @param group The NcGroup in which this SDDPBlock will be serialized. */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------------------------------------------------------------*/

 /// serialize the random cuts
 /** This function serializes the random cuts in the file with the given name
  * (path). If \p filename is empty, then no operation is performed. The file
  * will have the following netCDF format:
  *
  * - The dimension "TimeHorizon" containing the number of stages.
  *
  * - The dimension "NumberScenarios" containing the number of scenarios.
  *
  * - The group "PolyhedralFunction_t_s", for each t in {0, ..., TimeHorizon -
  *   1} and s in {0, ..., NumberScenarios - 1}, containing the serialization
  *   of the PolyhedralFunction representing the random cuts associated with
  *   stage t and scenario s. Each individual group is optional. If the group
  *   "PolyhedralFunction_t_s" is not provided, then the PolyhedralFunction
  *   associated with stage t and scenario s will not be loaded (which means
  *   it will have no cuts).
  *
  * @param filename The name of the file in which the random cuts will be
  *        serialized. */

 void serialize_random_cuts( const std::string & filename ) const;

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SDDPBlock
    @{ */

 /// returns the time horizon
 /** This function returns the time horizon associated with this SDDPBlock.
  *
  * @return The time horizon.
  */
 virtual Index get_time_horizon() const {
  return v_Block.size() / num_sub_blocks_per_stage;
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of PolyhedralFunction
 /** This function returns the vector of PolyhedralFunction associated with
  * this SDDPBlock.
  *
  * @return The vector of PolyhedralFunction.
  */
 const std::vector< PolyhedralFunction * > &
 get_polyhedral_functions() const {
  return v_polyhedral_functions;
 }

/*--------------------------------------------------------------------------*/

 /// returns a PolyhedralFunction
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

  if( ! num_polyhedral_per_sub_block )
   // Well, this is a funny SDDPBlock that has no PolyhedralFunction.
   return nullptr;

  assert( stage < get_time_horizon() );
  assert( sub_block_index < num_sub_blocks_per_stage );
  assert( i < num_polyhedral_per_sub_block );

  const auto index = ( stage * num_sub_blocks_per_stage + sub_block_index )
   * num_polyhedral_per_sub_block + i;

  return v_polyhedral_functions[ index ];
 }

/*--------------------------------------------------------------------------*/

 /// returns the number of PolyhedralFunction in each sub-Block
 /** This function returns the number of PolyhedralFunction present in each
  * sub-Block.
  *
  * @return The number of PolyhedralFunction in each sub-Block.
  */
 Index get_num_polyhedral_function_per_sub_block() const {
  return num_polyhedral_per_sub_block;
 }

/*--------------------------------------------------------------------------*/

 /// returns the number of sub-Blocks for each stage
 /** This function returns the number of sub-Blocks for each stage.
  *
  * @return The number of sub-Blocks for each stage.
  */
 Index get_num_sub_blocks_per_stage() const {
  return num_sub_blocks_per_stage;
 }

/*--------------------------------------------------------------------------*/

 /// returns a sub-Block of this SDDPBlock
 /** This function returns the sub-Block of index \p sub_block_index at the
  * given \p stage of this SDDPBlock. The given \p stage must be an integer
  * between 0 and get_time_horizon() - 1 and the index of the sub-Block must
  * be an integer between 0 and get_num_sub_blocks_per_stage() 0 - 1. If any
  * of them is an invalid index, an exception is thrown.
  *
  * @param stage The stage of the desired sub-Block. It must be a number
  *        between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the desired sub-Block at the given \p
  *        stage. It must be a number between 0 and
  *        get_num_sub_blocks_per_stage() - 1.
  *
  * @return The sub-Block of this SDDPBlock associated with the given \p stage
  *         and having index \p sub_block_index.
  */
 virtual StochasticBlock * get_sub_Block( Index stage ,
                                          Index sub_block_index = 0 ) const;

/*--------------------------------------------------------------------------*/

 /// returns the admissible state associated with the given \p stage
 /** This function returns an iterator to the vector containing the admissible
  * state for the given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @return An iterator to the vector containing the admissible state for the
  *         given \p stage. */

 std::vector< double >::const_iterator get_admissible_state( Index stage ) const {
  assert( stage < get_time_horizon() );
  return std::next( admissible_states.cbegin() ,
                    admissible_state_begin[ stage ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of the admissible state associated with the given \p stage
 /** This function returns the size of the admissible state associated with
  * the given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @return The size of the admissible state for the given \p stage. */

 Index get_admissible_state_size( Index stage ) const {
  assert( stage < get_time_horizon() );
  if( stage == get_time_horizon() - 1 )
   return admissible_states.size() - admissible_state_begin[ stage ];
  else
   return admissible_state_begin[ stage + 1 ] - admissible_state_begin[ stage ];
 }

/*--------------------------------------------------------------------------*/

 /// returns the set of scenarios
 /** This function returns the set of scenarios. */
 const ScenarioSet & get_scenario_set() const {
  return scenario_set;
 }

/*--------------------------------------------------------------------------*/

 /// returns the initial state for the first stage problem
 /** This function returns the initial state for the first stage problem. */
 const std::vector< double > & get_initial_state() const {
  return initial_state;
 }

/*--------------------------------------------------------------------------*/

 /// returns a random cut
 /** This function returns the PolyhedralFunction representing the random cut
  * associated with the given \p stage and the scenario whose index is \p
  * scenario_index. The \p stage argument must be between 0 and
  * get_time_horizon() - 1 while \p scenario_index must be between 0 and
  * get_scenario_set().size() - 1. */

 PolyhedralFunction & get_random_cut( Index stage , Index scenario_index ) {
  if( stage >= random_cuts.size() )
   throw( std::invalid_argument( "SDDPBlock::get_random_cut: no random cut "
                                 "for stage " + std::to_string( stage ) ) );

  if( scenario_index >= random_cuts[ stage ].size() )
   throw( std::invalid_argument
          ( "SDDPBlock::get_random_cut: no random cut for scenario index " +
            std::to_string( scenario_index ) + "." ) );

  return random_cuts[ stage ][ scenario_index ];
 }

/*--------------------------------------------------------------------------*/

 /// returns the sense of the Objective of the SDDPBlock
 /** This function returns the sense of the Objective of the SDDPBlock, which
  * is defined to be the sense of the Objective of its first inner Block. If
  * this SDDPBlock has no inner Block, this function returns
  * Objective::eUndef.
  *
  * @return the sense of the Objective of the first inner Block of this
  *         SDDPBlock if there is one. Otherwise, it returns
  *         Objective::eUndef. */

 int get_objective_sense() const override;

/**@} ----------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of an SDDPBlock
 * @{ */

 /// add cuts to a sub-Block at the given stage
 /** This function adds cuts to the sub-Block with index \p sub_block_index at
  * the given \p stage. The parameters must satisfy the following
  * requirements:
  *
  * 1. \p A must be a matrix with as many columns as there are cuts to be
  *    added and the number of rows must be equal to the number of Variable
  *    defined in the BendersBlock associated with stage \p stage.
  *
  * 2. \p b must be a vector whose size is equal to the number of
  *    rows of \p A. The cuts are given by Ax + b.
  *
  * 3. \p stage must be an integer between 0 and get_time_horizon() - 1.
  *
  * 4. \p sub_block_index must be an integer between 0 and
  *    get_num_sub_blocks_per_stage() - 1.
  *
  * @param A A matrix containing the coefficients of the cuts to be added.
  *
  * @param b A vector containing the constants of the cuts to be added.
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param sub_block_index The index of the sub-Block whose cuts will be
  *        updated.
  *
  * @param remove_current_cuts If true, all cuts currently part of sub-Block
  *        with index \p sub_block_index at the given \p stage are removed
  *        before the given cuts are added. If false, the current cuts are
  *        kept. */

 void add_cuts( PolyhedralFunction::MultiVector && A ,
                PolyhedralFunction::RealVector && b , Index stage ,
                Index sub_block_index ,
                Index number_cuts_to_keep = Inf< Index >() );

/*--------------------------------------------------------------------------*/
 /// add cuts to all sub-Blocks at the given stage
 /** This function adds cuts to all subs-Blocks at the given \p stage. The
  * parameters must satisfy the following requirements:
  *
  * 1. \p A must be a matrix with as many columns as there are cuts to be
  *    added and the number of rows must be equal to the number of Variable
  *    defined in the BendersBlock associated with stage \p stage.
  *
  * 2. \p b must be a vector whose size is equal to the number of
  *    rows of \p A. The cuts are given by Ax + b.
  *
  * 3. \p stage must be an integer between 0 and get_time_horizon() - 1.
  *
  * @param A A matrix containing the coefficients of the cuts to be added.
  *
  * @param b A vector containing the constants of the cuts to be added.
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param remove_current_cuts If true, all cuts currently part of sub-Block
  *        with index \p sub_block_index at the given \p stage are removed
  *        before the given cuts are added. If false, the current cuts are
  *        kept. */

 void add_cuts( PolyhedralFunction::MultiVector && A ,
                PolyhedralFunction::RealVector && b , Index stage ,
                Index number_cuts_to_keep = Inf< Index >() ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i ) {
   auto A_ = A;
   auto b_ = b;
   add_cuts( std::move( A_ ) , std::move( b_ ) , stage , i ,
             number_cuts_to_keep );
  }
 }

/*--------------------------------------------------------------------------*/
  /// removes all cuts from each PolyhedralFunction
  /** This function removes all cuts from each PolyhedralFunction. */

 void remove_cuts() {
  for( Index stage = 0 ; stage < get_time_horizon() ; ++stage )
   for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i )
    for( Index sub_block_index = 0 ;
         sub_block_index < num_sub_blocks_per_stage ; ++sub_block_index )
     get_polyhedral_function( stage , i , sub_block_index )->delete_rows();
 }

/*--------------------------------------------------------------------------*/
 /// store the given random cut
 /** This function store the random cut given by \p coefficients and \p alpha,
  * which must be associated with the given \p stage and with the scenario
  * whose index is \p scenario_index.
  *
  * @param coefficients The coefficients of the cut.
  *
  * @param alpha The constant of the cut.
  *
  * @param stage The stage (a number between 0 and get_time_horizon() - 1)
  *        associated with the given cut.
  *
  * @param scenario_index The index (a number between 0 and
  *        get_scenario_set().size() - 1) of the scenario associated with the
  *        given cut. */

 void store_random_cut( std::vector< double > && coefficients , double alpha ,
                        Index stage , Index scenario_index );

/*--------------------------------------------------------------------------*/
 /// returns the number of cuts currently present at the given \p stage
 /** This function returns the number of cuts currently present in the i-th
  * PolyhedralFunction of the sub-Block with index \p sub_block_index at the
  * given \p stage. The parameter \p i is the index of the PolyhedralFunction
  * to which the cuts should be added (its default value is 0).
  *
  * @param stage The stage whose cuts should be updated.
  *
  * @param sub_block_index The index of a sub-Block at the given \p stage.
  *
  * @param i The index of the PolyhedralFunction in the indicated sub-Block.
  *
  * @return The number of cuts currently present in the i-th
  *         PolyhedralFunction of the sub-Block with index \p sub_block_index
  *         at the given \p stage. */

 Index get_number_cuts( Index stage , Index sub_block_index , Index i = 0 )
  const {

  if( stage >= get_time_horizon() )
   throw( std::invalid_argument( "SDDPBlock::get_num__cuts: invalid stage "
                                 "index: " + std::to_string( stage ) ) );

  return get_polyhedral_function( stage , i , sub_block_index )->get_nrows();
 }

/*--------------------------------------------------------------------------*/

 /// returns the current future cost of the given sub-Block at the given stage
 double get_future_cost( Index stage , Index sub_block_index ) const;

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage. The size of the \p values array parameter must be equal to the
  * number N of state Variable of the problem at the given \p stage, so that
  * the value of the i-th state Variable will be values( i ), for each i in
  * {0, ..., N-1}.
  *
  * @param values The Eigen::ArrayXd containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. This must be an
  *              integer between 0 and get_time_horizon() - 1.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_state( const Eigen::ArrayXd & values , Index stage ,
                 Index sub_block_index );

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of all sub-Blocks at the given stage
 /** This function sets the values of the state Variable of all sub-Blocks at
  * the given \p stage. The size of the \p values array parameter must be
  * equal to the number N of state Variable of the problem at the given \p
  * stage, so that the value of the i-th state Variable will be values( i ),
  * for each i in {0, ..., N-1}.
  *
  * @param values The Eigen::ArrayXd containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. This must be an
  *              integer between 0 and get_time_horizon() - 1. */

 void set_state( const Eigen::ArrayXd & values , Index stage ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i )
   set_state( values , stage , i );
 }

/*--------------------------------------------------------------------------*/

 /// returns the values of the state Variable of the problem at the given stage
 /** This function returns the current values of the state Variable of the
  * problem associated with the sub-Block with index \p sub_block_index at the
  * given \p stage.
  *
  * @param stage The stage whose state Variable values are desired.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1.
  *
  * @return The current values of the state Variable of the problem at the
  *         given \p stage. */

 std::vector< double > get_state( Index stage ,
                                  Index sub_block_index = 0 ) const;

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage. The size of the \p values array parameter must be equal to the
  * number N of state Variable of the problem at the given \p stage, so that
  * the value of the i-th state Variable will be values[ i ], for each i in
  * {0, ..., N-1}.
  *
  * @param values The vector containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_state( const std::vector< double > & values , Index stage ,
                 Index sub_block_index );

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of all sub-Blocks at the given stage
 /** This function sets the values of the state Variable of all sub-Blocks at
  * the given \p stage. The size of the \p values array parameter must be
  * equal to the number N of state Variable of the problem at the given \p
  * stage, so that the value of the i-th state Variable will be values[ i ],
  * for each i in {0, ..., N-1}.
  *
  * @param values The vector containing the values of the Variable.
  *
  * @param stage The stage whose state Variable must be set. */

 void set_state( const std::vector< double > & values , Index stage ) {
  for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i )
   set_state( values , stage , i );
 }

/*--------------------------------------------------------------------------*/

 /// sets the values of the state Variable of the problem at the given stage
 /** This function sets the values of the state Variable of the problem
  * associated with the sub-Block with index \p sub_block_index at the given
  * \p stage, according to the admissible state of this SDDPBlock.
  *
  * @param stage The stage whose state must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_admissible_state( Index stage , Index sub_block_index = 0 );

/*--------------------------------------------------------------------------*/

 /// updates the sub-Block at the given stage for the given scenario
 /** This function updates the sub-Block whose index is \p sub_block_index at
  * the given \p stage for the given \p scenario.
  *
  * @param scenario_id The id of the scenario that must be set.
  *
  * @param sub_block_index The index of the sub-Block at the given \p
  *        stage. This must be an integer between 0 and
  *        get_num_sub_blocks_per_stage() - 1. */

 void set_scenario( Index scenario_id , Index stage ,
                    Index sub_block_index = 0 ) {
  auto sub_scenario_begin = scenario_set.
   sub_scenario_begin( scenario_id , stage );

  try {
   get_sub_Block( stage , sub_block_index )->set_data( sub_scenario_begin );
  }
  catch( const std::exception & e ) {
   std::cout << "SDDPBlock::set_scenario: exception while setting scenario "
             << scenario_id << " of stage " << stage << ".\n"
             << e.what() << std::endl;
   std::exit( EXIT_FAILURE );
  }
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/
 /// Pointers to the PolyhedralFunction of each sub-Block
 /** This vector stores the pointers to the PolyhedralFunction of each
  * sub-Block of this SDDPBlock. The pointer to the i-th PolyhedralFunction of
  * the j-th sub-Block of stage t is located at position
  *
  *   ( t * num_sub_blocks_per_stage + j ) * num_polyhedral_per_sub_block + i
  */
 std::vector< PolyhedralFunction * > v_polyhedral_functions;

 /// Number of PolyhedralFunctions for each sub-Block
 Index num_polyhedral_per_sub_block = 1;

 /// Number of sub-Blocks for each stage
 Index num_sub_blocks_per_stage = 1;

 /// Simulator for the forward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_forward;

 /// Simulator for the backward step of the SDDP method
 std::shared_ptr< ScenarioSimulator > simulator_backward;

 /// The set of scenarios
 ScenarioSet scenario_set;

 /// The start index of each admissible state
 /** For each t in {0, ..., TimeHorizon - 1}, admissible_state_begin[ t ] is
  * the index in vector #admissible_states at which the admissible state for
  * stage t begins. */

 std::vector< Index > admissible_state_begin;

 /// A vector containing the concatenation of states for each time instant
 std::vector< double > admissible_states;

 /// An initial state for the first stage problem
 std::vector< double > initial_state;

 /// Random cuts for each stage and each scenario
 /** This boost::multi_array stores the random cuts for all stages and all
  * scenarios. A random cut is a cut associated with a particular scenario. */
 boost::multi_array< PolyhedralFunction , 2 > random_cuts;

 /// It indicates whether the random cuts have been initialized
 bool f_random_cuts_initialized = false;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// deserializes the i-th sub-Block out of the given group
 /** This auxiliary function deserializes the \p i-th sub-Block out of the
  * given \p group.
  *
  * @param group The netCDF::NcGroup containing the description of the
  *        sub-Block.
  *
  * @param i The index of the sub-Block to be deserialized.
  *
  * @return A pointer to the Block that was deserialized.
  */
 Block * deserialize_sub_Block( const netCDF::NcGroup & group , Index i ) {
  std::string sub_group_name = "StochasticBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   sub_group = group.getGroup( "StochasticBlock" );
   if( sub_group.isNull() )
    throw std::logic_error( "SDDPBlock::deserialize: neither group '" +
                            sub_group_name + "' nor 'StochasticBlock' "
                            "was found." );
   sub_group_name = "StochasticBlock";
  }

  auto type = sub_group.getAtt( "type" );
  if( type.isNull() )
   throw std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                           sub_group_name + "' must be present." );

  std::string type_name;
  type.getValues( type_name );

  if( type_name != "StochasticBlock" )
   throw std::logic_error( "SDDPBlock::deserialize: attribute 'type' of '" +
                           sub_group_name + "' must contain "
                           "'StochasticBlock'." );

  auto sub_Block = new_Block( sub_group , this );

  if( ! sub_Block )
   throw std::logic_error( "SDDPBlock::deserialize: sub-group '" +
                           sub_group_name + "' is incomplete." );

  if( sub_group_name != "StochasticBlock" ) {

   // If StochasticBlock_i does not have sub-group "Block" then
   // "StochasticBlock" must have one.
   if( sub_group.getGroup( "Block" ).isNull() ) {

    auto StochasticBlock_group = group.getGroup( "StochasticBlock" );
    if( StochasticBlock_group.isNull() )
     throw std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was not "
                             "provided neither in '" + sub_group_name +
                             "' nor in 'StochasticBlock'" );


    auto Block_group = StochasticBlock_group.getGroup( "Block" );
    if( Block_group.isNull() )
     throw std::logic_error( "SDDPBlock::deserialize: sub-group 'Block' was not "
                             "provided neither in '" + sub_group_name +
                             "' nor in 'StochasticBlock'" );

    auto inner_block = new_Block( Block_group, this );
    if( ! inner_block )
     throw std::logic_error( "SDDPBlock::deserialize: the 'Block' sub-group of "
                             "the 'StochasticBlock' group has an invalid or "
                             "incomplete description." );

    static_cast< StochasticBlock * >( sub_Block )->
     set_inner_block( inner_block );
   }

   // If StochasticBlock_i does not have the description of vector of
   // "DataMapping" then, if "StochasticBlock" has one, we use it.

   Index num_data_mappings;
   if( ! ::SMSpp_di_unipi_it::deserialize_dim( sub_group , "NumberDataMappings" ,
                                               num_data_mappings , true ) ) {

    auto StochasticBlock_group = group.getGroup( "StochasticBlock" );
    if( ! StochasticBlock_group.isNull() ) {

     if( ::SMSpp_di_unipi_it::deserialize_dim( StochasticBlock_group ,
                                               "NumberDataMappings" ,
                                               num_data_mappings , true ) ) {

      std::vector< std::unique_ptr< SimpleDataMappingBase > > data_mappings;
      data_mappings.reserve( num_data_mappings );
      SimpleDataMappingBase::deserialize
       ( group , data_mappings ,
         static_cast< StochasticBlock *>( sub_Block )->get_inner_block() );

      static_cast< StochasticBlock * >( sub_Block )->
       set_data_mappings( std::move( data_mappings ) );
     }
    }
   }
  }

  return sub_Block;
 }

 /*--------------------------------------------------------------------------*/

 /// initializes the structure that stores the random cuts
 void initialize_random_cuts();

/*--------------------------------------------------------------------------*/

};   // end( class SDDPBlock )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* SDDPBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SDDPBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
