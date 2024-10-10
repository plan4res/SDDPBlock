/*--------------------------------------------------------------------------*/
/*-------------------------- File SDDPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPSolver class.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright Copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BendersBlock.h"
#include "FRealObjective.h"
#include "Objective.h"
#include "SDDPBlock.h"
#include "SDDPSolver.h"
#include "StochasticBlock.h"

#include <chrono>
#include <Eigen/Core>

#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"

#include "StOpt/sddp/backwardForwardSDDP.h"
#include "StOpt/sddp/LocalConstRegressionForSDDP.h"
#include "StOpt/sddp/LocalLinearRegressionForSDDP.h"

#include "StOpt/sddp/LocalConstRegressionForSDDPGeners.h"
#include "StOpt/sddp/LocalLinearRegressionForSDDPGeners.h"

#define BENDERSBFUNCTION_DEBUG

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SDDPSolver );
SMSpp_insert_in_factory_cpp_0( SDDPSolverState );

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS of SDDPSolver -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPSolver::set_ComputeConfig( ComputeConfig * scfg ) {

 ThinComputeInterface::set_ComputeConfig( scfg );

 if( ! scfg ) { // factory reset
  delete f_inner_block_solver_config;
  f_inner_block_solver_config = nullptr;

  delete f_inner_block_config;
  f_inner_block_config = nullptr;

  delete f_get_var_solution_config;
  f_get_var_solution_config = nullptr;

  return;
 }

 if( ! scfg->f_extra_Configuration )
  // No extra Configuration has been provided. There is nothing else to do.
  return;

 BlockConfig * block_config = nullptr;
 BlockSolverConfig * block_solver_config = nullptr;

 // First, we try to extract a BlockConfig and/or a BlockSolverConfig from the
 // extra Configuration, as well as the Configuration to be passed to
 // get_var_solution() of the inner Solver.

 if( auto config = dynamic_cast< SimpleConfiguration<
     std::vector< Configuration * > > * >( scfg->f_extra_Configuration ) ) {

  // The extra Configuration is a vector. The first element of this vector, if
  // present and not nullptr, must be a BlockConfig for the inner Blocks of
  // the BendersBFunctions. The second element, if present and not nullptr,
  // must be a BlockSolverConfig for the inner Blocks of the
  // BendersBFunctions. Finally, the third element, if present and not
  // nullptr, must be a Configuration to be passed to get_var_solution() when
  // retrieving the Solutions to the inner Blocks of the BendersBFunctions.

  if( ( ! config->f_value.empty() ) && config->f_value.front() ) {
   // A BlockConfig must have been provided.
   if( ! ( block_config =
           dynamic_cast< BlockConfig * >( config->f_value.front() ) ) )
    throw( std::invalid_argument( "SDDPSolver::set_ComputeConfig: The first "
                                  "element of the extra Configuration is "
                                  "not a BlockConfig." ) );
  }

  if( config->f_value.size() >= 2 && config->f_value[ 1 ] ) {
   // A BlockSolverConfig must have been provided.
   if( ! ( block_solver_config =
           dynamic_cast< BlockSolverConfig * >( config->f_value[ 1 ] ) ) )
    throw( std::invalid_argument( "SDDPSolver::set_ComputeConfig: The second "
                                  "element of the extra Configuration is "
                                  "not a BlockSolverConfig." ) );
  }

  if( config->f_value.size() >= 3 && config->f_value[ 2 ] ) {
   // A Configuration for get_var_solution() of the Solver attached to the
   // inner Blocks.
   f_get_var_solution_config = config->f_value[ 2 ]->clone();
  }
 }
 else {
  // The extra Configuration must be either a BlockConfig or a
  // BlockSolverConfig.
  if( auto bc = dynamic_cast< BlockConfig * >( scfg->f_extra_Configuration ) )
   block_config = bc;
  else if( auto bsc =
           dynamic_cast< BlockSolverConfig * >( scfg->f_extra_Configuration ) )
   block_solver_config = bsc;
  else
   throw( std::invalid_argument( "SDDPSolver::set_ComputeConfig: The extra "
                                 "Configuration is invalid." ) );
 }

 // Now, replace the old Configurations if new ones have been provided.

 if( block_config ) {
  // A BlockConfig has been provided. Delete the old BlockConfig and clone the
  // given one.
  delete f_inner_block_config;
  f_inner_block_config = block_config->clone();
 }

 if( block_solver_config ) {
  // A BlockSolverConfig has been provided. Delete the old BlockSolverConfig
  // and clone the given one.
  delete f_inner_block_solver_config;
  f_inner_block_solver_config = block_solver_config->clone();
 }
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::set_Block( Block * block ) {
 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 Solver::set_Block( block );

 if( ! block )
  return;

 auto sddp_block = dynamic_cast< SDDPBlock * >( block );

 if( ! sddp_block )
  throw( std::invalid_argument( "SDDPSolver::set_Block: An SDDPSolver can "
                                "only be attached to an SDDPBlock." ) );

 const auto & scenario_set = sddp_block->get_scenario_set();
 sddp_optimizer->set_scenarios( scenario_set );

 // BlockConfig for the inner Blocks
 if( ( ! f_inner_block_config ) &&
     ( ! f_inner_block_config_filename.empty() ) ) {
  auto c = Configuration::deserialize( f_inner_block_config_filename );
  if( ! ( f_inner_block_config = dynamic_cast< BlockConfig * >( c ) ) ) {
   delete c;
   throw( std::invalid_argument
          ( "SDDPSolver::configure_inner_block: file " +
            f_inner_block_config_filename + " is not a BlockConfig." ) );
  }
 }

 // BlockSolverConfig for the inner Blocks
 if( ( ! f_inner_block_solver_config ) &&
     ( ! f_inner_block_solver_config_filename.empty() ) ) {
  auto c = Configuration::deserialize( f_inner_block_solver_config_filename );
  if( ! ( f_inner_block_solver_config =
          dynamic_cast< BlockSolverConfig * >( c ) ) ) {
   delete c;
   throw( std::invalid_argument
          ( "SDDPSolver::configure_inner_block: file " +
            f_inner_block_solver_config_filename +
            " is not a BlockSolverConfig." ) );
  }
 }

 // Configure the inner Blocks
 if( f_inner_block_config || f_inner_block_solver_config ) {

  for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {

   for( Index sub_block_index = 0 ;
        sub_block_index < sddp_block->get_num_sub_blocks_per_stage() ;
        ++sub_block_index ) {

    auto benders_function = get_benders_function( stage , sub_block_index );

    if( ! benders_function )
     throw( std::invalid_argument
            ( "SDDPSolver::set_Block: The BendersBFunction at stage " +
              std::to_string( stage ) + " is not present." ) );

    auto inner_block = benders_function->get_inner_block();

    if( ! inner_block )
     throw( std::invalid_argument
            ( "SDDPSolver::set_Block: The inner Block of the BendersBFunction "
              " at stage " + std::to_string( stage ) + " is not present." ) );

    if( f_inner_block_config )
     f_inner_block_config->apply( inner_block );

    if( f_inner_block_solver_config )
     f_inner_block_solver_config->apply( inner_block );
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int SDDPSolver::compute( bool changedvars ) {

 if( ! f_Block )
  return( kBlockLocked );

 // Possibly lock the SDDPBlock

 auto owned = f_Block->is_owned_by( f_id );        // check if already locked
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )  // if not try to lock
  return( kBlockLocked );                          // return error on failure

 process_outstanding_Modification();

 const auto time_horizon = get_time_horizon();

 // ostream for StOpt output
 boost::iostreams::stream< boost::iostreams::null_sink >
  null_sink( ( boost::iostreams::null_sink() ) );
 std::ostream * output_stream = & null_sink;
 if( f_log && log_verbosity >= 2 )
  output_stream = f_log;

 // log of the sub-Solvers

 std::vector< std::ofstream > sub_solvers_logfiles;

 if( ! f_sub_solver_filename_prefix.empty() ) {

  const auto num_sub_blocks_per_stage =
   static_cast< SDDPBlock * >( f_Block )->get_num_sub_blocks_per_stage();

  sub_solvers_logfiles.reserve( time_horizon * num_sub_blocks_per_stage );

  for( Index t = 0 ; t < time_horizon ; ++t ) {
   for( Index i = 0 ; i < num_sub_blocks_per_stage ; ++i ) {
    auto suffix = std::to_string( t );
    if( num_sub_blocks_per_stage > 1 )
     suffix += "-" + std::to_string( i );
    const auto filename = f_sub_solver_filename_prefix + suffix;

    sub_solvers_logfiles.emplace_back
     ( std::ofstream{ filename , std::ofstream::out | std::ofstream::app } );

    auto benders_function = get_benders_function( t , i );
    auto solver = benders_function->get_solver();
    solver->set_log( & sub_solvers_logfiles.back() );
   }
  }
 }

 /* "dates" must be an array with size T + 1, where T is the time_horizon,
  * such that dates[ t ] contains the t-th time step (in our case it is simply
  * t) for each t in {0, ..., T-1}. The element in this array, with index T,
  * is associated with the cut to be used at the last time step. */

 Eigen::ArrayXd dates =
  Eigen::ArrayXd::LinSpaced( time_horizon + 1 , 0 , time_horizon );

 /* As input to the StOpt SDDP solver, it contains the maximum number of
  * iterations that the solver should perform. As output, it contains the
  * number of iterations performed by the StOpt SDDP solver. */
 number_iterations_performed = maximum_number_iterations;

 /* As input to the StOpt SDDP solver, it contains the desired accuracy that
  * the method should seek. As output, it contains the accuracy achieved by
  * the StOpt SDDP solver, which is given by
  *
  * | backwardValue - forwardValue | / forwardValue
  *
  * where backwardValue is the value of the last backward pass and
  * forwardValueForConv is the value obtained during the forward pass when
  * checking for convergence. */

 auto accuracy_achieved_stopt = accuracy;

 /*****************/
 /* INITIAL STATE */
 /*****************/

 const auto number_state_variables = sddp_optimizer->getStateSize();

 if( ! initial_state.empty() ) {
  /* If an initial state for the first stage is provided, then the initial
   * state for the first stage is updated to the given one. */

  if( initial_state.size() !=
      decltype( initial_state )::size_type( number_state_variables ) ) {
   throw( std::logic_error
          ( "SDDPSolver::compute: the size of the given initial state (" +
            std::to_string( initial_state.size() ) + ") is different from the "
            "size of the admissible state ("
            + std::to_string( number_state_variables ) + ")" ) );
  }

  // Set the initial state for the first stage
  static_cast< SDDPBlock * >( f_Block )->set_state( initial_state , 0 );
 }
 else {
  const auto & state =
   static_cast< SDDPBlock * >( f_Block )->get_initial_state();
  if( ! state.empty() )
   // Use the initial state given by SDDPBlock
   static_cast< SDDPBlock * >( f_Block )->set_state( state , 0 );
 }

 /* The "initial state" that is passed to StOpt is the admissible state of the
  * first stage. This is because StOpt will use this state as the initial
  * state for the second stage problem (during the first backward
  * pass). However, StOpt also uses this state as the initial state for the
  * first stage problem. But this is not an issue because the initial state
  * for the first stage problem is set at most once here in compute() and it
  * is ignored in the oneStep* functions. */

 Eigen::ArrayXd initial_state_array = sddp_optimizer->oneAdmissibleState( 0 );

 /***************************/
 /* CUTS FOR THE LAST STAGE */
 /***************************/

 /* The cuts to be used at the last time instant. If any cuts are provided,
  * they are added to the last stage right now and do not need to be passed to
  * StOpt. Notice that any cut that is currently at the last stage are kept
  * there. */

 PolyhedralFunction::MultiVector A;
 PolyhedralFunction::RealVector b;

 if( ! last_stage_cuts.empty() ) {

  if( last_stage_cuts.size() % ( number_state_variables + 1 ) != 0 )
   throw( std::logic_error( "SDDPSolver::compute: Invalid size of the cuts "
                            "for the last stage." ) );

  const auto number_of_cuts =
   last_stage_cuts.size() / ( number_state_variables + 1 );

  // Pack the cuts

  A.resize( number_of_cuts );
  b.resize( number_of_cuts );

  for( Index i = 0 ; i < number_of_cuts ; ++i ) {
   const auto begin = i * ( number_state_variables + 1 );
   b[ i ] = last_stage_cuts[ begin + number_state_variables ];
   A[ i ].resize( number_state_variables );
   for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j )
    A[ i ][ j ] = last_stage_cuts[ begin + j ];
  }
 }
 else {
  /* No cut for the last stage has been provided. We check whether the
   * PolyhedralFunction at the last stage contains a finite bound or at least
   * one cut (row). */

  const auto polyhedral_function =
   static_cast< SDDPBlock * >( f_Block )->get_polyhedral_functions().back();

  if( ( ! polyhedral_function->is_bound_set() ) &&
      ( polyhedral_function->get_nrows() == 0 ) ) {
   if( f_log )
#ifdef USE_MPI
    if( ! mpi_communicator.rank() )
#endif
     *f_log << "Warning: SDDPSolver::compute: No cut for the last stage has "
            << "been provided and\nthe PolyhedralFunction at the last stage "
            << "has no bound and no row (cut). By\ndefault, the all-zero cut"
            << " will then be used for the last stage." << std::endl;
   b.resize( 1 , 0 );
   A.resize( 1 );
   A.front().resize( number_state_variables , 0 );
  }
 }

 if( ! A.empty() ) {
  // Add the cuts to the last stage
  static_cast< SDDPBlock * >( f_Block )->add_cuts
   ( std::move( A ) , std::move( b ) , get_time_horizon() - 1 );
 }

 /* The cuts for the subproblem at the last stage have just been added. StOpt
  * requires a StOpt::SDDPFinalCut as argument, representing the cuts for the
  * last stage. So, we create a dummy one, which will be ignored within
  * oneStepBackward() and oneStepForward() when the subproblem at the last
  * stage is being solved. */

 StOpt::SDDPFinalCut final_cut
  ( Eigen::ArrayXXd::Zero( number_state_variables + 1 , 1 ) );

 /* In some situations, we can avoid adding the cuts every time in
  * oneStepBackward() and oneStepForward(). One of these situations, for
  * instance, happens when no mesh discretization is provided. In this case,
  * the set of cuts to be considered in oneStepBackward() at any given time is
  * always a superset of that provided at all previous iterations for that
  * same time. So, we can add only the new cuts. The following vector stores
  * the number of cuts currently present at each stage and help determine
  * which cuts were generated by StOpt. */

 number_initial_cuts.resize( get_time_horizon() );
 for( Index stage = 0 ; stage < get_time_horizon() ; ++stage )
  number_initial_cuts[ stage ] =
   static_cast< SDDPBlock * >( f_Block )->get_number_cuts( stage , 0 );

 /***********************/
 /* MESH DISCRETIZATION */
 /***********************/

 // Meshes for regression
 Eigen::ArrayXi mesh_discretization_array;
 if( mesh_provided() ) {
  auto simulator = std::static_pointer_cast< ScenarioSimulator >
   ( sddp_optimizer->getSimulatorBackward() );
  const auto particle_length = simulator->get_particle_length();
  if( mesh_discretization.size() != particle_length )
   throw( std::logic_error
         ( "SDDPSolver::compute: simulation particle has dimension " +
           std::to_string( particle_length ) + " but given mesh discretization "
           "has dimension " + std::to_string( mesh_discretization.size() ) ) );
  mesh_discretization_array.resize( mesh_discretization.size() );
  for( Index i = 0 ; i < mesh_discretization.size() ; ++i )
   mesh_discretization_array( i ) = mesh_discretization[ i ];
 }

 /***************************/
 /* RESET THE SDDPOptimizer */
 /***************************/

 sddp_optimizer->reset();

 /***********************/
 /* SOLVING THE PROBLEM */
 /***********************/

 // Invoke the StOpt SDDP solver
 auto backward_forward_values =
  StOpt::backwardForwardSDDP< StOpt::LocalLinearRegressionForSDDP >
  ( sddp_optimizer , number_simulations_for_convergence , initial_state_array ,
    final_cut , dates , mesh_discretization_array , regressors_filename ,
    cuts_filename , visited_states_filename , number_iterations_performed ,
    accuracy_achieved_stopt , convergence_frequency , *output_stream ,
#ifdef USE_MPI
    mpi_communicator ,
#endif
    print_cpu_time );

 // Possibly output the future cost functions and/or save the State

 if( output_frequency > 0 ) {
  file_output();
 }

 // Unlock the SDDPBlock

 if( ! owned )              // if the Block was actually locked
  f_Block->unlock( f_id );  // unlock it

 // Retrieve the backward and forward values

 backward_value = backward_forward_values.first;
 forward_value = backward_forward_values.second;

 // Log

 if( f_log && log_verbosity > 0 ) {
#ifdef USE_MPI
  if( ! mpi_communicator.rank() ) {
#endif
   *f_log << "Backward value: " << std::setprecision( 20 )
          << backward_value << std::endl;
   *f_log << "Forward value:  " << std::setprecision( 20 )
          << forward_value << std::endl;
#ifdef USE_MPI
  }
#endif
 }

 // Close the log files of the sub-Solvers

 for( auto & logfile : sub_solvers_logfiles )
  logfile.close();

 // Compute the accuracy achieved

 if( forward_value != 0.0 )
  accuracy_achieved = std::abs( ( backward_value - forward_value ) /
                                forward_value );
 else
  accuracy_achieved = std::abs( backward_value );

 // Determine the status of SDDPSolver

 if( accuracy_achieved_stopt == 0.0 && accuracy_achieved != 0.0 )
  status = kCurveCross;
 else if( accuracy_achieved_stopt <= accuracy )
  status = kOK;
 else if( number_iterations_performed == maximum_number_iterations )
  status = kStopIter;
 else
  status = kError; // TODO

 return status;
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::get_lb( void ) {
 if( ! f_Block )
  return -Inf< double >();

 const auto minimization =
  ( f_Block->get_objective_sense() == Objective::eMin );

 if( status == kOK || status == kLowPrecision ) {
  if( minimization )
   return backward_value;
  return forward_value;
 }

 if( status == kInfeasible ) {
  if( minimization )
   return Inf< double >();
  return -Inf< double >();
 }

 return -Inf< double >();
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::get_ub( void ) {
 if( ! f_Block )
  return Inf< double >();

 const auto minimization =
  ( f_Block->get_objective_sense() == Objective::eMin );

 if( status == kOK || status == kLowPrecision ) {
  if( minimization )
   return forward_value;
  return backward_value;
 }

 if( status == kInfeasible ) {
  if( minimization )
   return Inf< double >();
  return -Inf< double >();
 }

 return Inf< double >();
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::process_outstanding_Modification() {
 v_mod.clear();
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/

void SDDPSolver::reset_event_handler( int type , EventID id ) {
 if( type != eEverykIteration )
  throw( std::invalid_argument( "SDDPSolver::reset_event_handler: unsupported"
                                " event type " + std::to_string( type ) ) );

 if( id >= v_events[ type ].size() )
  throw( std::invalid_argument( "SDDPSolver::reset_event_handler: incorrect "
                                "event id " + std::to_string( id ) +
                                " for type " + std::to_string( type ) ) );

 static auto do_nothing = []() -> int {
  return( ThinComputeInterface::eContinue ); };

 if( id == v_events[ type ].size() - 1 ) {
  // if the event is the last of its type, shorten the vector; moreover, if
  // any of the previous events is a do_nothing, keep shortening
  do
   v_events[ type ].pop_back();
  while( ( ! v_events[ type ].empty() ) &&
         ( *( v_events[ type ].back().target < int( * )() > ( ) ) ==
           do_nothing ) );
 }
 else
  // the event is not the last of its type: replace it with a do_nothing to
  // avoid messing up with the id-s, which are positions in the vector
  v_events[ type ][ id ] = do_nothing;
}

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR HANDLING THE State OF THE SDDPSolver ------------*/
/*--------------------------------------------------------------------------*/

State * SDDPSolver::get_State( void ) const {
 return new SDDPSolverState( this );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::put_State( const State & state ) {

 // If this SDDPSolver is not currently attached to an SDDPBlock, nothing is
 // done
 auto sddp_block = static_cast< SDDPBlock * >( f_Block );
 if( ! sddp_block )
  return;

 auto s = dynamic_cast< const SDDPSolverState & >( state );

 const auto time_horizon = get_time_horizon();

 const auto num_polyhedral_per_sub_block =
  sddp_block->get_num_polyhedral_function_per_sub_block();

 const auto num_sub_blocks_per_stage =
  sddp_block->get_num_sub_blocks_per_stage();

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
   for( Index sub_block_index = 0 ;
        sub_block_index < num_sub_blocks_per_stage ; ++sub_block_index ) {

    auto polyhedral_function =
     sddp_block->get_polyhedral_function( t , i , sub_block_index );

    assert( polyhedral_function );

    auto A = s.v_A[ t ];
    auto b = s.v_b[ t ];

    polyhedral_function->set_PolyhedralFunction
     ( std::move( A ) , std::move( b ) , s.v_bound[ t ] , s.v_is_convex[ t ] );
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::put_State( State && state ) {

 // If this SDDPSolver is not currently attached to an SDDPBlock, nothing is
 // done
 auto sddp_block = static_cast< SDDPBlock * >( f_Block );
 if( ! sddp_block )
  return;

 auto s = dynamic_cast< SDDPSolverState && >( state );

 const auto time_horizon = get_time_horizon();

 const auto num_polyhedral_per_sub_block =
  sddp_block->get_num_polyhedral_function_per_sub_block();

 const auto num_sub_blocks_per_stage =
  sddp_block->get_num_sub_blocks_per_stage();

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
   for( Index sub_block_index = 0 ;
        sub_block_index < num_sub_blocks_per_stage ; ++sub_block_index ) {

    auto polyhedral_function =
     sddp_block->get_polyhedral_function( t , i , sub_block_index );

    assert( polyhedral_function );

    polyhedral_function->set_PolyhedralFunction
     ( std::move( s.v_A[ t ] ) , std::move( s.v_b[ t ] ) , s.v_bound[ t ] ,
       s.v_is_convex[ t ] );
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::serialize_State
( netCDF::NcGroup & group , const std::string & sub_group_name ) const {

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );
 if( ! sddp_block )
  return;

 if( ! sub_group_name.empty() ) {
  auto sub_group = group.addGroup( sub_group_name );
  if( sub_group.isNull() )
   throw( std::invalid_argument( "SDDPSolver::serialize_State: it was not "
                                 "possible to add group " + sub_group_name ) );
  serialize_State( sub_group );
  return;
 }

 group.putAtt( "type" , "SDDPSolverState" );

 const auto time_horizon = get_time_horizon();
 group.addDim( "TimeHorizon" , get_time_horizon() );

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  const auto polyhedral_function = sddp_block->get_polyhedral_function( t );
  assert( polyhedral_function );
  auto num_var = polyhedral_function->get_num_active_var();
  auto is_convex = polyhedral_function->is_convex();
  auto bound = polyhedral_function->get_global_bound();
  const auto & A = polyhedral_function->get_A();
  const auto & b = polyhedral_function->get_b();
  SDDPSolverState::serialize( group , t , num_var , is_convex , bound , A , b );
 }
}  // end( SDDPSolver::serialize_State )

/*--------------------------------------------------------------------------*/

SDDPBlock::Index SDDPSolver::get_time_horizon( void ) const {
 return static_cast< SDDPBlock * >( f_Block )->get_time_horizon();
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::add_cuts( const Eigen::ArrayXXd & cuts ,
                           SDDPBlock::Index stage ,
                           SDDPBlock::Index sub_block_index ) const {

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::add_cuts: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 // Number of cuts that will be added. In principle, all given cuts are added.

 auto number_cuts_to_be_added = cuts.cols();

 /* If a mesh has been provided, all cuts are replaced by the given ones, but
  * the initial cuts (those there were present at the time StOpt was called)
  * are kept. */
 auto number_cuts_to_keep = number_initial_cuts[ stage ];

 if( ! mesh_provided() ) {

  // Total number of cuts previously added by StOpt
  auto number_cuts_previously_added =
   static_cast< SDDPBlock * >( f_Block )->
   get_number_cuts( stage , sub_block_index ) - number_initial_cuts[ stage ];

  assert( cuts.cols() >= number_cuts_previously_added );

  /* StOpt provides all cuts that were ever generated. Since no mesh has been
   * provided, we can consider only the most recent cuts (since the previous
   * ones have already been added before). We are assuming that the new cuts
   * provided by StOpt appear in the last columns of the matrix "cuts". */

  if( cuts.cols() == number_cuts_previously_added )
   return; // all cuts are already there

  // Add only the most recent cuts
  number_cuts_to_be_added = cuts.cols() - number_cuts_previously_added;

  // and keep the current ones
  number_cuts_to_keep = Inf< Index >();
 }

 // Store the given cuts in A and b

 PolyhedralFunction::MultiVector A( number_cuts_to_be_added );
 PolyhedralFunction::RealVector b( number_cuts_to_be_added );

 /* Each column in "cuts" contains a cut. The first element is the constant
  * term of the cut, and all the other elements are the coefficients. */

 for( Index k = 0 ; k < number_cuts_to_be_added ; ++k ) {
  const auto i = number_cuts_to_be_added - k - 1;
  const auto col = cuts.cols() - k - 1;
  A[ i ].resize( cuts.rows() - 1 );
  b[ i ] = cuts( 0 , col );
  for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j ) {
   A[ i ][ j ] = cuts( j + 1 , col );
  }
 }

 // Log

 if( f_log && log_verbosity >= 30 ) {
  *f_log << "  Adding the following cuts:" << std::endl;
  for( decltype( b )::size_type i = 0 ; i < b.size() ; ++i ) {
   *f_log << "    (" << b[ i ];
   for( decltype( A[ i ].size() ) j = 0 ; j < A[ i ].size() ; ++j )
    *f_log << ", " << A[ i ][ j ];
   *f_log << ")" << std::endl;
  }
 }

 // Finally add the cuts

 static_cast< SDDPBlock * >( f_Block )->add_cuts
  ( std::move( A ) , std::move( b ) , stage , sub_block_index ,
    number_cuts_to_keep );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::set_state( const Eigen::ArrayXd & state ,
                            SDDPBlock::Index stage ,
                            SDDPBlock::Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::set_state: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_state( state , stage ,
                                                   sub_block_index );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::set_scenario( SDDPBlock::Index scenario_id ,
                               SDDPBlock::Index stage ,
                               SDDPBlock::Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::set_scenario: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_scenario( scenario_id , stage ,
                                                      sub_block_index );
}

/*--------------------------------------------------------------------------*/

BendersBFunction *
SDDPSolver::get_benders_function( SDDPBlock::Index stage ,
                                  SDDPBlock::Index sub_block_index ) const {
 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< SDDPBlock * >( f_Block )->get_sub_Block
    ( stage , sub_block_index )->get_nested_Blocks().front() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 return static_cast< BendersBFunction * >( objective->get_function() );
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::solve( SDDPBlock::Index stage ,
                          SDDPBlock::Index sub_block_index ) {

 /* Solving the subproblem consists in evaluating the Objective of the
  * BendersBFunction associated with the subproblem of the given stage. */

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::solve: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 auto benders_block = static_cast< BendersBlock * >
  ( static_cast< SDDPBlock * >( f_Block )->get_sub_Block
    ( stage , sub_block_index )->get_nested_Blocks().front() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 auto benders_function = static_cast< BendersBFunction * >
  ( objective->get_function() );

 auto status = benders_function->compute();

 if( status != kOK && status != kLowPrecision ) {
  // No feasible solution has been found

  std::string message;
  if( status == kUnbounded )
   message = " The subproblem is unbounded.";
  else if( status == kInfeasible )
   message = " The subproblem is infeasible.";
  else if( status == kError )
   message = " An error occurred while solving the subproblem.";
  else if( status == kStopTime )
   message = " A time limit has been reached while solving the subproblem.";
  else if( status == kStopIter )
   message = " A maximum number of iterations has been reached while solving "
    "the subproblem.";

  throw( std::logic_error( "SDDPSolver::solve: the subproblem at stage " +
                           std::to_string( stage ) + " was not solved." +
                           message ) );
 }

 auto solver = benders_function->get_solver();
 if( ! solver->has_var_solution() ) {
  throw( std::logic_error( "SDDPSolver::solve: the subproblem at stage " +
                           std::to_string( stage ) + " has no solution." ) );
 }

 solver->get_var_solution( f_get_var_solution_config );

 return benders_function->get_value();
}

/*--------------------------------------------------------------------------*/

template< class T >
T SDDPSolver::get_solution( SDDPBlock::Index stage ,
                            SDDPBlock::Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPSolver::get_solution: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 const auto polyhedral_function =
  static_cast< SDDPBlock * >( f_Block )->get_polyhedral_function
  ( stage , 0 , sub_block_index );

 T solution( polyhedral_function->get_num_active_var() );

 auto data = solution.data();
 for( const auto & variable : * polyhedral_function ) {
  *data = static_cast< const ColVariable & >( variable ).get_value();
  data++;
 }
 return solution;
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::file_output() const {

#ifdef USE_MPI
 if( mpi_communicator.rank() )
  return;
#endif

 if( ! f_output_filename.empty() ) {
  // Output the future cost functions
  std::string cuts_filename = f_output_filename;
  if( f_add_suffix )
   cuts_filename += f_filename_suffix;
  output_future_cost_functions( cuts_filename );
 }

 if( ! f_state_filename.empty() ) {
  // Serialize the State
  std::string state_filename = f_state_filename;
  if( f_add_suffix )
   state_filename += f_filename_suffix;
  serialize_State( state_filename );
 }

 if( ! f_random_cuts_filename.empty() ) {
  // Output the random cuts
  std::string random_cuts_filename = f_random_cuts_filename;
  if( f_add_suffix )
   random_cuts_filename += f_filename_suffix;
  auto sddp_block = static_cast< SDDPBlock *>( f_Block );
  sddp_block->serialize_random_cuts( random_cuts_filename );
 }

 f_add_suffix = ! f_add_suffix;
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::output_future_cost_functions( const std::string & filename )
 const {

 if( filename.empty() )
  return;

 auto sddp_block = static_cast< SDDPBlock *>( f_Block );

 const auto & functions = sddp_block->get_polyhedral_functions();
 if( functions.empty() )
  return;

 std::ofstream output( filename , std::ios::out );

 const char separator_character = ',';
 const auto num_var = functions.front()->get_num_active_var();

 output << "Timestep";
 for( Index i = 0 ; i < num_var ; ++i ) {
  output << separator_character << "a_" << std::to_string( i );
 }
 output << separator_character << "b" << std::endl;

 for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {

  auto function = sddp_block->get_polyhedral_function( stage , 0 , 0 );

  const auto & b = function->get_b();
  const auto & A = function->get_A();

  assert( b.size() == A.size() );

  for( Index i = 0 ; i < b.size() ; ++i ) {
   output << stage;
   for( Index j = 0 ; j < A[ i ].size() ; ++j )
    output << separator_character << std::setprecision( 20 ) << A[ i ][ j ];
   output << separator_character << std::setprecision( 20 ) << b[ i ]
          << std::endl;
  }
 }

 output.close();

}

/*--------------------------------------------------------------------------*/

void SDDPSolver::serialize_State( const std::string & filename ) const {

 if( filename.empty() )
  return;

 netCDF::NcFile file( filename , netCDF::NcFile::replace );

 auto group = file.addGroup( "SDDPSolverState" );

 serialize_State( group );

 file.close();
}

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of SDDPOptimizer -----------------------*/
/*--------------------------------------------------------------------------*/

Eigen::ArrayXd SDDPSolver::SDDPOptimizer::oneStepBackward
( const StOpt::SDDPCutOptBase & sddp_cut ,
  const std::tuple< std::shared_ptr< Eigen::ArrayXd > , int , int > & state ,
  const Eigen::ArrayXd & particle , const int & simulation_id ) const {

 /* The last argument is the simulation id indicating in which scenario the
  * resolution will be done. */

 const auto scenario_index = get_backward_scenario_index( simulation_id );

 // TODO Work with multiple sub-Blocks
 const Index sub_block_index = 0;
 const bool scenario_must_be_set = true;

 const auto backward = SDDPSolver::SDDPOptimizer::oneStepBackward
  ( sddp_cut , state , particle , simulation_id , scenario_index ,
    scenario_must_be_set , sub_block_index );

 previous_pass_was_backward = true;

 return backward;
}

/*--------------------------------------------------------------------------*/

Eigen::ArrayXd SDDPSolver::SDDPOptimizer::oneStepBackward
( const StOpt::SDDPCutOptBase & sddp_cut ,
  const std::tuple< std::shared_ptr< Eigen::ArrayXd > , int , int > & state ,
  const Eigen::ArrayXd & particle , const int & simulation_id ,
  const Index scenario_index , const bool scenario_must_be_set ,
  const Index sub_block_index ) const {

 const auto start_time = std::chrono::system_clock::now();

 const auto current_stage = get_current_backward_stage();

 /* In the first stage, the scenario will not be set if the given scenario
  * index is negative. */

 const bool scenario_provided = ( current_stage > 0 ) ||
  ( scenario_index < Inf< Index >() );

 // Log

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepBackward *****" << std::endl;
  *log << "  Stage:          " << current_stage << std::endl;
  *log << "  Block index:    " << sub_block_index << std::endl;
  *log << "  Scenario index: ";
  if( scenario_provided ) *log << scenario_index << std::endl;
  else *log << "none" << std::endl;
  if( current_stage > 0 ) {
   *log << "  Simulation id:  " << simulation_id << std::endl;
   if( sddp_solver->log_verbosity >= 4 ) {
    *log << "  Particle:       (";
    for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
     if( i > 0 ) *log << ", ";
     *log << particle( i );
    }
    *log << ")" << std::endl;
   }
  }

  if( current_stage > 0 && sddp_solver->log_verbosity >= 10 ) {
   *log << "  State:          (";
   const auto & state_variables = * std::get<0>( state );
   for( decltype( state_variables.size() ) i = 0 ;
        i < state_variables.size() ; ++i ) {
    if( i > 0 ) *log << ", ";
    *log << state_variables( i );
   }
   *log << ")" << std::endl;
  }
 }

 /***************/
 /* ADDING CUTS */
 /***************/

 if( current_stage < sddp_solver->get_time_horizon() - 1 ) {
  /* The cuts for the last stage are added only once in the beginning of
   * compute(). */

  const auto cuts = sddp_cut.getCutsAssociatedToTheParticle
   ( std::get<1>( state ) );
  sddp_solver->add_cuts( cuts , current_stage , sub_block_index );
 }

 /*******************/
 /* STATE VARIABLES */
 /*******************/

 if( current_stage > 0 ) {
  /* The initial state for the first stage problem is set only once, in the
   * beggining of compute(). */
  sddp_solver->set_state( * std::get<0>( state ).get() , current_stage ,
                          sub_block_index );
 }

 /***************/
 /* RANDOM DATA */
 /***************/

 /* In the first stage, the scenario is not set if the given scenario index
  * (for the first stage) is negative. */
 if( scenario_must_be_set && scenario_provided )
  sddp_solver->set_scenario( scenario_index , current_stage , sub_block_index );

 /**************************/
 /* SOLVING THE SUBPROBLEM */
 /**************************/

 auto objective_value = sddp_solver->solve( current_stage , sub_block_index );

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  *( sddp_solver->f_log ) << "  Objective:      " << objective_value
                          << std::endl;
 }

 /**********************************/
 /* CONSTRUCTING THE LINEARIZATION */
 /**********************************/

 /* The oneStepBackard function returns a one-dimensional array (called
  * "linearization" and declared below) whose size is the number of state
  * variables plus one and that contains a linearization of the
  * BendersBFunction. The first component contains a value that will be used
  * by StOpt to compute the linearization constant and the remaining
  * components contain the coefficients of the linearization of the
  * BendersBFunction. For i in {1, ..., number_state_variables},
  * linearization( i ) contains the coefficient of the linearization of the
  * BendersBFunction associated with the i-th state variable.
  *
  * StOpt will compute the linearization constant based on the value that is
  * returned at the first position of the "linearization" array. By letting z
  * the value that we return in linearization(0), StOpt will compute the
  * linearization constant as z - g'y, where g is the vector with the
  * linearization coefficients and y is the vector with the values of the
  * state variables. Therefore, we return in linearization(0) the value alpha
  * + g'y, where alpha is the linearization constant provided by the
  * BendersBFunction, i.e., we return a lower bound (upper bound if the
  * subproblem is a maximization) on the value of the solution. */

 const auto number_state_variables = std::get<0>( state )->size();
 Eigen::ArrayXd linearization( number_state_variables + 1 );

 auto benders_function =
  sddp_solver->get_benders_function( current_stage , sub_block_index );

 if( benders_function->has_linearization( true ) ) {

  // Retrieve the linearization coefficients
  benders_function->get_linearization_coefficients( linearization.data() + 1 );

  // Compute g'y
  double gy = 0;

  if( current_stage == 0 ) {
   /* Retrieve the initial state from the SDDPBlock, as the "state" parameter
    * does not contain the state for the first stage problem. */
   const auto state_var = static_cast< SDDPBlock * >
    ( sddp_solver->f_Block )->get_state( current_stage , sub_block_index );
   for( decltype( state_var.size() ) j = 0 ; j < state_var.size() ; ++j )
    gy += linearization( j + 1 ) * state_var[ j ];
  }
  else {
   // Use the state variables given as argument
   const auto state_var = std::get<0>( state ).get();
   for( Index j = 0 ; j < state_var->size() ; ++j )
    gy += linearization( j + 1 ) * ( *state_var )( j );
  }

  // Retrieve the linearization constant
  const auto alpha = benders_function->get_linearization_constant();

  linearization( 0 ) = alpha + gy;

  // Store the random cut if required
  if( sddp_solver->store_random_cuts() ) {
   auto sddp_block = static_cast< SDDPBlock * >( sddp_solver->f_Block );
   std::vector< double > coefficients
    ( linearization.data() + 1 , linearization.data() + linearization.size() );

   auto actual_scenario_index = scenario_index;
   if( actual_scenario_index == Inf< Index >() ) {
    // The scenario index is infinity. Thus, this must be the first stage.
    assert( current_stage == 0 );

    // Since there is no scenario associated with the first stage, we store
    // the cut as if it were associated with the first scenario.
    actual_scenario_index = 0;
   }

   sddp_block->store_random_cut( std::move( coefficients ) , alpha ,
                                 current_stage , actual_scenario_index );
  }

 /*************/
 /* DEBUGGING */
 /*************/

  // Debugging the BendersBFunction

#ifdef BENDERSBFUNCTION_DEBUG
  check_linearization( objective_value , alpha , gy, linearization ,
                       sub_block_index );
#endif

 }
 else {
  // No diagonal linearization is available
  throw( std::logic_error( "SDDPOptimizer::oneStepBackward: no "
                           "linearization is available." ) );
 }

 /********************/
 /* FINAL LOG OUTPUT */
 /********************/

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {

  const auto end_time = std::chrono::system_clock::now();
  const std::chrono::duration< double > duration = end_time - start_time;
  const auto time = duration.count();
  const auto log = sddp_solver->f_log;
  *log << "  Time (s):       " << time << std::endl;

  if( sddp_solver->log_verbosity >= 10 ) {
   auto solution = sddp_solver->get_solution( current_stage , sub_block_index );
   *log << "  Solution:       (";
   for( decltype( solution.size() ) i = 0 ; i < solution.size() ; ++i ) {
    if( i > 0 ) *( sddp_solver->f_log ) << ", ";
    *log << solution( i );
   }
   *log << ")" << std::endl;
  }
 }

 return linearization;
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::SDDPOptimizer::oneStepForward
( const Eigen::ArrayXd & particle , Eigen::ArrayXd & state ,
  Eigen::ArrayXd & state_to_store , const StOpt::SDDPCutOptBase & sddp_cut ,
  const int & simulation_id ) const {

 // Update control variables and output

 if( previous_pass_was_backward ) {
  if( sddp_solver->output_frequency > 0 &&
      ( current_iteration % sddp_solver->output_frequency == 0 ) ) {
   sddp_solver->file_output();
  }

  if( sddp_solver->f_handle_events_every_k_iter &&
      ( current_iteration % sddp_solver->f_handle_events_every_k_iter == 0 ) )
   for( auto & event : sddp_solver->v_events[ eEverykIteration ] )
    event();

  current_iteration++;
 }

 const auto scenario_index = get_forward_scenario_index( simulation_id );

 // TODO Work with multiple sub-Blocks
 const Index sub_block_index = 0;
 const auto scenario_must_be_set = true;

 const auto forward = SDDPSolver::SDDPOptimizer::oneStepForward
  ( particle , state , state_to_store , sddp_cut , simulation_id ,
    scenario_index , scenario_must_be_set , sub_block_index );

 // Update control variable
 previous_pass_was_backward = false;

 return forward;
}

/*--------------------------------------------------------------------------*/

double SDDPSolver::SDDPOptimizer::oneStepForward
( const Eigen::ArrayXd & particle , Eigen::ArrayXd & state ,
  Eigen::ArrayXd & state_to_store , const StOpt::SDDPCutOptBase & sddp_cut ,
  const int & simulation_id , const Index scenario_index ,
  const bool scenario_must_be_set , const Index sub_block_index ) const {

 const auto start_time = std::chrono::system_clock::now();

 const auto current_stage = get_current_forward_stage();

 /* In the first stage, the scenario will not be set if the given scenario
  * index is negative. */

 const bool scenario_provided = ( current_stage > 0 ) ||
  ( scenario_index < Inf< Index >() );

 // Log

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  auto log = sddp_solver->f_log;
  *log << "***** SDDPSolver::SDDPOptimizer::oneStepForward *****"
            << std::endl;
  *log << "  Stage:          " << current_stage << std::endl;
  *log << "  Block index:    " << sub_block_index << std::endl;
  *log << "  Scenario index: ";
  if( scenario_provided ) *log << scenario_index << std::endl;
  else *log << "none" << std::endl;
  if( current_stage > 0 ) {
   *log << "  Simulation id:  " << simulation_id << std::endl;
   if( sddp_solver->log_verbosity >= 4 ) {
    *log << "  Particle:       (";
    for( decltype( particle.size() ) i = 0 ; i < particle.size() ; ++i ) {
     if( i > 0 ) *log << ", ";
     *log << particle( i );
    }
    *log << ")" << std::endl;
   }
  }

  if( current_stage > 0 && sddp_solver->log_verbosity >= 10 ) {
   *log << "  State:          (";
   for( decltype( state.size() ) i = 0 ; i < state.size() ; ++i ) {
    if( i > 0 ) *log << ", ";
    *log << state( i );
   }
   *log << ")" << std::endl;
  }
 }

 /***************/
 /* ADDING CUTS */
 /***************/

 if( sddp_solver->mesh_provided() && ( current_stage > 0 )  &&
     ( current_stage < sddp_solver->get_time_horizon() - 1 ) ) {
  /* Cuts are added if and only if meshes were provided and the current stage
   * is not the first or the last one. There is a single mesh associated with
   * the first stage and the cuts have already been added during the backward
   * pass. The cuts for the last stage are added once in the beginning of
   * compute(). The array "particle" contains the random quantities on which
   * the regression over the expectation of the value function will be
   * based. */

  const auto cuts = sddp_cut.getCutsAssociatedToAParticle( particle );
  sddp_solver->add_cuts( cuts , current_stage , sub_block_index );
 }

 /***************/
 /* RANDOM DATA */
 /***************/

 if( scenario_must_be_set && scenario_provided )
  sddp_solver->set_scenario( scenario_index , current_stage , sub_block_index );

 /*********************************/
 /* VARIABLES FROM PREVIOUS STAGE */
 /*********************************/

 if( current_stage > 0 ) {
  /* The initial state for the first stage problem is set only once, in the
   * beggining of compute(). */
  sddp_solver->set_state( state , current_stage , sub_block_index );
 }

 /**************************/
 /* SOLVING THE SUBPROBLEM */
 /**************************/

 auto objective_value = sddp_solver->solve( current_stage , sub_block_index );

 /* The objective_value takes into account the value of the future cost
  * function. For all stages other than the last one, we subtract the value of
  * the future cost function from objective_value. The subproblem at the last
  * stage, however, has a fixed future cost and it is kept as it is considered
  * part of the cost of that stage. */

 if( current_stage < sddp_solver->get_time_horizon() - 1 )
  objective_value -= static_cast< SDDPBlock * >( sddp_solver->f_Block )->
   get_future_cost( current_stage , sub_block_index );

 /**************************/
 /* RETRIVING THE SOLUTION */
 /**************************/

 // Retrieve the solution x_t of the Block associated with the current stage.

 auto solution = sddp_solver->get_solution( current_stage , sub_block_index );

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  *( sddp_solver->f_log ) << "  Objective:      " << objective_value
                          << std::endl;
  if( sddp_solver->log_verbosity >= 10 ) {
   *( sddp_solver->f_log ) << "  Solution:       (";
   for( decltype( solution.size() ) i = 0 ; i < solution.size() ; ++i ) {
    if( i > 0 ) *( sddp_solver->f_log ) << ", ";
    *( sddp_solver->f_log ) << solution( i );
   }
   *( sddp_solver->f_log ) << ")" << std::endl;
  }
 }

 // Store in state the current state, i.e, ( x_t, w_t^{dep} ).

 state.resize( solution.size() );
 state << solution;

 // Store in state_to_store the vector ( x_t, w_{t-1}^{dep} ).

 state_to_store.resize( solution.size() );
 state_to_store << solution;

 /********************/
 /* FINAL LOG OUTPUT */
 /********************/

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 3 ) {
  const auto end_time = std::chrono::system_clock::now();
  const std::chrono::duration< double > duration = end_time - start_time;
  const auto time = duration.count();
  *( sddp_solver->f_log ) << "  Time (s):       " << time << std::endl;
 }

 /*************************/
 /* RETURN SOLUTION VALUE */
 /*************************/

 return objective_value;
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::SDDPOptimizer::updateDates
( const double & date , const double & date_next ) {

  // We assume that the given arguments correspond to consecutive stages
  assert( date_next == date + 1.0 );

  // We assume that -1 <= date < T
  assert( - 1.0 <= date && date < sddp_solver->get_time_horizon() );

  this->date = date;
  this->date_next = date_next;
}

/*--------------------------------------------------------------------------*/

Eigen::ArrayXd
SDDPSolver::SDDPOptimizer::oneAdmissibleState( const double & stage ) {

 const auto sddp_block = static_cast< SDDPBlock * >( sddp_solver->f_Block );

 auto state_size = sddp_block->get_admissible_state_size( stage );
 Eigen::ArrayXd state( state_size );

 auto state_iterator = sddp_block->get_admissible_state( stage );

 auto data = state.data();
 for( decltype( state_size ) i = 0 ; i < state_size ;
      ++i , ++data , ++state_iterator )
  *data = *state_iterator;

 return state;
}

/*--------------------------------------------------------------------------*/

int SDDPSolver::SDDPOptimizer::getStateSize() const {
 return static_cast< SDDPBlock * >( sddp_solver->f_Block )->
  get_admissible_state_size( 0 );
}

/*--------------------------------------------------------------------------*/

void SDDPSolver::SDDPOptimizer::check_linearization
( double objective_value , double alpha , double gy ,
  const Eigen::ArrayXd & linearization , Index sub_block_index ) const {

 if( sddp_solver->f_log && sddp_solver->log_verbosity >= 20 ) {
  *( sddp_solver->f_log ) << "  Linearization: " << std::endl;
  *( sddp_solver->f_log ) << "    alpha:        " << alpha << std::endl;
  if( sddp_solver->log_verbosity >= 30 ) {
   *( sddp_solver->f_log ) << "    coefficients: (";
   for( decltype( linearization.size() ) i = 1 ;
        i < linearization.size() ; ++i ) {
    if( i > 1 )
     *( sddp_solver->f_log ) << ", ";
    *( sddp_solver->f_log ) << linearization( i );
   }
   *( sddp_solver->f_log ) << ")" << std::endl;
  }
 }

 const double epsilon = 1.0e-4;
 const auto scale =
  std::max( 1.0 , std::min( abs( objective_value ) , abs( alpha + gy ) ) );
 const auto diff = std::abs( objective_value - ( alpha + gy ) );
 if( diff > epsilon * scale ) {
  auto log = sddp_solver->f_log;
  if( ! log )
   log = &std::cerr;

  *log << "  SDDPOptimizer::oneStepBackward: linearization precision "
       << "was not achieved:" << std::endl;
  *log << "    precision required: " << std::setprecision( 20 )
       << epsilon << std::endl;
  *log << "    precision achieved: " << std::setprecision( 20 )
       << ( diff / scale ) << std::endl;
  *log << "    objective: " << std::setprecision( 20 )
       << objective_value << std::endl;
  *log << "    alpha:     " << std::setprecision( 20 )
       << alpha << std::endl;
  *log << "    g'y:       " << std::setprecision( 20 ) << gy << std::endl;
 }
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of SDDPSolverState -----------------------*/
/*--------------------------------------------------------------------------*/

SDDPSolverState::SDDPSolverState( const SDDPSolver * solver ) {

 if( ! solver )
  return;

 auto sddp_block = static_cast< SDDPBlock * >( solver->f_Block );

 if( ! sddp_block )
  return;

 const auto time_horizon = solver->get_time_horizon();

 v_is_convex.reserve( time_horizon );
 v_num_var.reserve( time_horizon );
 v_A.reserve( time_horizon );
 v_b.reserve( time_horizon );
 v_bound.reserve( time_horizon );

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  auto polyhedral_function = sddp_block->get_polyhedral_function( t );
  assert( polyhedral_function );
  v_is_convex.push_back( polyhedral_function->is_convex() );
  v_num_var.push_back( polyhedral_function->get_num_active_var() );
  v_A.push_back( polyhedral_function->get_A() );
  v_b.push_back( polyhedral_function->get_b() );
  v_bound.push_back( polyhedral_function->get_global_bound() );
 }
}

/*--------------------------------------------------------------------------*/

void SDDPSolverState::serialize
( netCDF::NcGroup & group , Index t , Index num_var , bool is_convex ,
  PolyhedralFunction::FunctionValue bound ,
  const PolyhedralFunction::MultiVector & A ,
  const PolyhedralFunction::RealVector & b ) {

 if( is_convex )
  group.addDim( "PolyFunction_sign_" + std::to_string( t ) , 1 );
 else
  group.addDim( "PolyFunction_sign_" + std::to_string( t ) , 0 );

 ( group.addVar( "PolyFunction_lb_" + std::to_string( t ) ,
                 netCDF::NcDouble() ) ).putVar( &bound );

 auto nv = group.addDim( "PolyFunction_NumVar_" + std::to_string( t ) ,
                         num_var );

 auto num_rows = b.size();

 if( num_rows ) {

  auto nr = group.addDim( "PolyFunction_NumRow_" + std::to_string( t ) ,
                          num_rows );

  auto ncdA = group.addVar( "PolyFunction_A_" + std::to_string( t ) ,
                            netCDF::NcDouble() , { nr , nv } );

  for( Index i = 0 ; i < num_rows ; ++i )
   ncdA.putVar( { i , 0 } , { 1 , num_var } , A[ i ].data() );

  ( group.addVar( "PolyFunction_b_" + std::to_string( t ) ,
                  netCDF::NcDouble() , nr ) ).
   putVar( { 0 } , { num_rows } , b.data() );
 }
}

/*--------------------------------------------------------------------------*/

void SDDPSolverState::serialize( netCDF::NcGroup & group ) const {

 State::serialize( group );

 const auto time_horizon = v_num_var.size();

 group.addDim( "TimeHorizon" , time_horizon );

 for( Index t = 0 ; t < time_horizon ; ++t )
  serialize( group , t , v_num_var[ t ] , v_is_convex[ t ] , v_bound[ t ] ,
             v_A[ t ] , v_b[ t ] );
}

/*--------------------------------------------------------------------------*/

void SDDPSolverState::deserialize( const netCDF::NcGroup & group ) {

 auto TimeHorizon = group.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() )
  throw( std::logic_error
         ( "SDDPSolverState::deserialize: TimeHorizon dimension is "
           "required, but it is not in the given group." ) );

 auto time_horizon = TimeHorizon.getSize();

 v_is_convex.resize( time_horizon );
 v_A.resize( time_horizon );
 v_b.resize( time_horizon );
 v_bound.resize( time_horizon );
 v_num_var.resize( time_horizon );

 for( decltype( time_horizon ) t = 0 ; t < time_horizon ; ++t ) {

  auto nv = group.getDim( "PolyFunction_NumVar_" + std::to_string( t ) );
  if( nv.isNull() )
   throw( std::logic_error( "SDDPSolverState::deserialize: PolyFunction_NumVar_"
                            + std::to_string( t ) + " dimension is required,"
                            " but it is not in the given group.") );

  v_num_var[ t ] = nv.getSize();

  auto nr = group.getDim( "PolyFunction_NumRow_" + std::to_string( t ) );
  if( ( ! nr.isNull() ) && ( nr.getSize() ) ) {
   auto ncdA = group.getVar( "PolyFunction_A_" + std::to_string( t ) );
   if( ncdA.isNull() )
    throw( std::logic_error( "SDDPSolverState::deserialize: PolyFunction_A_"
                             + std::to_string( t ) + " dimension is required,"
                             " but it is not in the given group.") );

   auto ncdb = group.getVar( "PolyFunction_b_" + std::to_string( t ) );
   if( ncdb.isNull() )
    throw( std::logic_error( "SDDPSolverState::deserialize: PolyFunction_b_"
                             + std::to_string( t ) + " dimension is required,"
                             " but it is not in the given group.") );

   v_A[ t ].resize( nr.getSize() );
   for( Index i = 0 ; i < v_A[ t ].size() ; ++i ) {
    v_A[ t ][ i ].resize( v_num_var[ t ] );
    ncdA.getVar( { i , 0 } , { 1 , v_num_var[ t ] } , v_A[ t ][ i ].data() );
   }

   v_b[ t ].resize( nr.getSize() );
   ncdb.getVar( v_b[ t ].data() );
  }

  v_is_convex[ t ] = true;
  auto sgn = group.getDim( "PolyFunction_sign_" + std::to_string( t ) );
  if( ! sgn.isNull() )
   v_is_convex[ t ] = sgn.getSize() > 0 ? true : false;

  auto nclb = group.getVar( "PolyFunction_lb_" + std::to_string( t ) );
  if( nclb.isNull() ) {
   if( v_is_convex[ t ] )
    v_bound[ t ] = -Inf< PolyhedralFunction::FunctionValue >();
   else
    v_bound[ t ] = Inf< PolyhedralFunction::FunctionValue >();
  }
  else
   nclb.getVar( & v_bound[ t ] );
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File SDDPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
