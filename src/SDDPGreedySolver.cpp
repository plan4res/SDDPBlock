/*--------------------------------------------------------------------------*/
/*----------------------- File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPGreedySolver class.
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
#include "BlockSolverConfig.h"
#include "CDASolver.h"
#include "FRealObjective.h"
#include "SDDPBlock.h"
#include "SDDPGreedySolver.h"
#include "StochasticBlock.h"

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SDDPGreedySolver );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of SDDPGreedySolver ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING SDDPGreedySolver --------------*/
/*--------------------------------------------------------------------------*/

SDDPGreedySolver::~SDDPGreedySolver() {
 delete f_inner_block_config;
 delete f_inner_block_solver_config;
 delete f_get_var_solution_config;
 delete f_get_dual_solution_config;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_ComputeConfig( ComputeConfig * scfg ) {

 ThinComputeInterface::set_ComputeConfig( scfg );

 if( ! scfg ) { // factory reset
  delete f_inner_block_solver_config;
  f_inner_block_solver_config = nullptr;

  delete f_inner_block_config;
  f_inner_block_config = nullptr;

  delete f_get_var_solution_config;
  f_get_var_solution_config = nullptr;

  delete f_get_dual_solution_config;
  f_get_dual_solution_config = nullptr;

  return;
 }

 if( ! scfg->f_extra_Configuration ) {
  // No extra Configuration has been provided.

  if( f_early_config && f_Block &&
      ( ! ( f_inner_block_config_filename.empty() &&
            f_inner_block_solver_config_filename.empty() ) ) ) {
   // If required, configure the SDDPBlock right away.
   const auto time_horizon = get_time_horizon();
   for( Index stage = 0 ; stage < time_horizon ; ++stage )
    configure_inner_block( stage );
  }

  // There is nothing else to do.
  return;
 }

 BlockConfig * block_config = nullptr;
 BlockSolverConfig * block_solver_config = nullptr;

 // First, we try to extract a BlockConfig and/or a BlockSolverConfig from the
 // extra Configuration, as well as the Configuration to be passed to
 // get_var_solution() and get_dual_solution() of the inner Solver.

 if( auto config = dynamic_cast< SimpleConfiguration<
     std::vector< Configuration * > > * >( scfg->f_extra_Configuration ) ) {

  // The extra Configuration is a vector. The first element of this
  // vector, if present and not nullptr, must be a BlockConfig for the inner
  // Blocks of the BendersBFunctions. The second element, if present and not
  // nullptr, must be a BlockSolverConfig for the inner Blocks of the
  // BendersBFunctions. The third element, if present and not nullptr, must be
  // a Configuration to be passed to get_var_solution() when retrieving the
  // Solutions to the inner Blocks of the BendersBFunctions. Finally, the
  // fourth element, if present and not nullptr, must be a Configuration to be
  // passed to get_dual_solution() when retrieving the Solutions to the inner
  // Blocks of the BendersBFunctions.

  if( ( ! config->f_value.empty() ) && config->f_value.front() ) {
   // A BlockConfig must have been provided.
   if( ! ( block_config =
           dynamic_cast< BlockConfig * >( config->f_value.front() ) ) )
    throw( std::invalid_argument( "SDDPGreedySolver::set_ComputeConfig: The "
                                  "first element of the extra Configuration "
                                  "is not a BlockConfig." ) );
  }

  if( config->f_value.size() >= 2 && config->f_value[ 1 ] ) {
   // A BlockSolverConfig must have been provided.
   if( ! ( block_solver_config =
           dynamic_cast< BlockSolverConfig * >( config->f_value[ 1 ] ) ) )
    throw( std::invalid_argument( "SDDGreedyPSolver::set_ComputeConfig: The "
                                  "second element of the extra Configuration "
                                  "is not a BlockSolverConfig." ) );
  }

  if( config->f_value.size() >= 3 && config->f_value[ 2 ] ) {
   // A Configuration for get_var_solution() of the Solver attached to the
   // inner Blocks.
   f_get_var_solution_config = config->f_value[ 2 ]->clone();
  }

  if( config->f_value.size() >= 4 && config->f_value[ 3 ] ) {
   // A Configuration for get_dual_solution() of the Solver attached to the
   // inner Blocks.
   f_get_dual_solution_config = config->f_value[ 3 ]->clone();
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
   throw( std::invalid_argument( "SDDPGreedySolver::set_ComputeConfig: The "
                                 "extra Configuration is invalid." ) );
 }

 // Now, replace the old Configurations if new ones have been provided.

 if( block_config ) {
  // A BlockConfig has been provided. Delete the old BlockConfig and clone the
  // given one.

  if( f_inner_block_config &&
      std::any_of( v_inner_block_configured.cbegin() ,
                   v_inner_block_configured.cend() ,
                   []( auto b ) { return b; } ) ) {

   f_inner_block_config ->clear();

   const auto time_horizon = get_time_horizon();
   for( Index stage = 0 ; stage < time_horizon ; ++stage ) {
    auto benders_function = get_benders_function( stage );
    auto inner_block = benders_function->get_inner_block();
    f_inner_block_config->apply( inner_block );
   }
  }

  delete f_inner_block_config;
  f_inner_block_config = block_config->clone();
  v_inner_block_configured.assign( get_time_horizon() , false );
 }

 if( block_solver_config ) {
  // A BlockSolverConfig has been provided. Delete the old BlockSolverConfig
  // and clone the given one.

  if( f_inner_block_solver_config &&
      std::any_of( v_inner_solver_configured.cbegin() ,
                   v_inner_solver_configured.cend() ,
                   []( auto b ) { return b; } ) ) {

   f_inner_block_solver_config ->clear();

   const auto time_horizon = get_time_horizon();
   for( Index stage = 0 ; stage < time_horizon ; ++stage ) {
    auto benders_function = get_benders_function( stage );
    auto inner_block = benders_function->get_inner_block();
    f_inner_block_solver_config->apply( inner_block );
   }
  }

  delete f_inner_block_solver_config;
  f_inner_block_solver_config = block_solver_config->clone();
  v_inner_solver_configured.assign( get_time_horizon() , false );
 }

 if( f_early_config && f_Block ) {
  // If required, configure the SDDPBlock right away.
  const auto time_horizon = get_time_horizon();
  for( Index stage = 0 ; stage < time_horizon ; ++stage )
   configure_inner_block( stage );
 }

}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_Block( Block * block ) {

 if( f_Block == block )
  return;

 if( f_Block ) {
  // TODO clean
  v_inner_block_configured.clear();
  v_inner_solver_configured.clear();
 }

 Solver::set_Block( block );

 if( ! f_Block )
  return;

 SDDPBlock * sddp_block;
 if( ! ( sddp_block = dynamic_cast< SDDPBlock * >( block ) ) )
  throw( std::invalid_argument( "SDDPGreedySolver::set_Block: given Block "
                                "is not an SDDPBlock." ) );

 v_inner_block_configured.assign( sddp_block->get_time_horizon() , false );
 v_inner_solver_configured.assign( sddp_block->get_time_horizon() , false );

}  // end( SDDPGreedySolver::set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int SDDPGreedySolver::compute( bool changedvars ) {

 reset_compute_time();

 if( ! f_Block )
  return kBlockLocked;

 // Possibly lock the SDDPBlock

 auto owned = f_Block->is_owned_by( f_id );        // check if already locked
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )  // if not try to lock
  return( kBlockLocked );                          // return error on failure

 process_outstanding_Modification();

 const auto time_horizon = get_time_horizon();
 status_compute = Solver::kOK;
 fault_stage = Inf< Index >();
 solution_value = 0.0;
 f_has_var_solution = false;
 f_has_dual_solution = true;

 // Initial state for the first stage problem
 if( time_horizon > 0 ) {
  if( ! v_initial_state.empty() )
   // Use the initial state given as parameter to the SDDPGreedySolver
   set_state( v_initial_state , 0 );
  else {
   const auto & state =
    static_cast< SDDPBlock * >( f_Block )->get_initial_state();
   if( ! state.empty() )
    // Use the initial state given by SDDPBlock
    set_state( state , 0 );
  }
 }

 // Possibly set the scenario of the first stage
 if( f_first_stage_scenario_id == -1 )
  set_scenario( f_scenario_id , 0 );
 else if( f_first_stage_scenario_id >= 0 )
  set_scenario( f_first_stage_scenario_id , 0 );

 // If required, load the random cuts
 if( ! f_random_cuts_filename.empty() )
  static_cast< SDDPBlock * >( f_Block )->
   deserialize_random_cuts( f_random_cuts_filename );

 // Clear the data from previous call to compute()
 f_simulation_data.clear();

 // Header of the log
 Logger logger( this , f_log , log_verbosity );
 logger.log_header();

 for( Index stage = 0 ; stage < time_horizon ; ++stage ) {

  reset_subproblem_time();

  // If required, sample a scenario for this stage.
  if( f_seed < Inf< Index >() )
   sample_scenario( stage );

  logger.log( stage , get_scenario_id( stage ) );

  if( stage > 0 ) {
   // Set the state of the subproblem as that given by the solution of the
   // subproblem at the previous stage.
   set_state( get_solution( stage - 1 ) , stage );

   // Set the scenario.
   set_scenario( get_scenario_id( stage ) , stage );
  }

  if( callback ) callback( stage );

  configure_inner_block( stage );

  if( ( ! f_load_cuts_once ) ||
      ( stage >= v_cuts_loaded.size() ) || ( ! v_cuts_loaded[ stage ] ) )
   load_cuts( stage );

  auto sub_status = solve( stage , true );

  if( sub_status == Solver::kInfeasible ) {
   const auto obj_sign = get_benders_function( stage )->is_convex() ? - 1 : 1;
   logger.log( - obj_sign * Inf< double >() );

   fault_stage = stage;
   if( stage == 0 ) status_compute = kInfeasible;
   else status_compute = kSubproblemInfeasible;
   break;
  }
  else if( sub_status == Solver::kUnbounded ) {
   const auto obj_sign = get_benders_function( stage )->is_convex() ? - 1 : 1;
   logger.log( obj_sign * Inf< double >() );

   fault_stage = stage;
   status_compute = Solver::kUnbounded;
   break;
  }
  else if( sub_status >= Solver::kError ) {
   logger.log();
   fault_stage = stage;
   status_compute = kError;
   break;
  }
  else if( ! get_sub_solver( stage )->has_var_solution() ) {
   logger.log();
   fault_stage = stage;
   status_compute = kSolutionNotFound;
   break;
  }
  else if( sub_status == Solver::kStopTime || sub_status == Solver::kStopIter ) {
   const auto sub_solution_value = get_sub_solution_value( stage );
   solution_value += sub_solution_value;
   if( fault_stage == Inf< Index >() ) {
    fault_stage = stage;
    status_compute = sub_status;
   }

   logger.log( sub_solution_value , get_future_value( stage ) );
  }
  else {
   const auto sub_solution_value = get_sub_solution_value( stage );
   solution_value += sub_solution_value;

   logger.log( sub_solution_value , get_future_value( stage ) );

   if( ! f_simulation_data_filename.empty() ) {
    store_simulation_data( stage , get_scenario_id( stage ) );
   }
  }

  if( f_unregister_solver )
   unregister_solver_inner_block( stage );
 }

 // Unlock the SDDPBlock

 if( ! owned )              // if the Block was actually locked
  f_Block->unlock( f_id );  // unlock it

 // Has a feasible solution been found?

 f_has_var_solution = ( status_compute == Solver::kOK )
  || ( status_compute == Solver::kLowPrecision )
  || ( status_compute == Solver::kStopIter )
  || ( status_compute == Solver::kStopTime );

 // Output the data obtained during the simulation if required.

 output_simulation_data( f_simulation_data_filename );

 // Final log

 logger.show_status();

 return status_compute;
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING THE DATA ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::get_var_solution( Configuration *solc ) {
 /* During the call to compute(), every Block associated with a stage in {0,
  * ..., T-2} has its solutions written in it. Therefore, we only need to
  * write the solution associated with the last stage. */

 auto solver = get_sub_solver( get_time_horizon() - 1 );

 if( ! solver )
  return; // The Solver must have been unregistered (but the Solution should
          // have already been written into the Block)

 if( ! solver->has_var_solution() )
  throw( std::logic_error( "SDDPGreedySolver::get_var_solution: subproblem "
                           "at the last stage does not have a solution." ) );
 else {
  solver->get_var_solution( f_get_var_solution_config );
 }
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::has_dual_solution() {
 return f_has_dual_solution;
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::is_dual_feasible() {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t ) {
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) ) {
   if( ! solver->is_dual_feasible() )
    return false;
  }
  else
   return false; // The sub-Solver is not a CDASolver
 }
 return true;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::get_dual_solution( Configuration * solc ) {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t )
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) )
   solver->get_dual_solution( solc );
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::new_dual_solution() {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t ) {
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) )
   if( solver->new_dual_solution() )
    return true;
 }
 return false;
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::has_dual_direction() {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t ) {
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) ) {
   if( ! solver->has_dual_direction() )
    return false;
  }
  else
   return false; // The sub-Solver is not a CDASolver
 }
 return true;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::get_dual_direction( Configuration * dirc ) {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t )
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) )
   solver->get_dual_direction( dirc );
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::new_dual_direction() {
 const auto time_horizon = get_time_horizon();
 for( Index t = 0 ; t < time_horizon ; ++t ) {
  if( auto solver = dynamic_cast< CDASolver * >( get_sub_solver( t ) ) )
   if( solver->new_dual_direction() )
    return true;
 }
 return false;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::configure_inner_block( Index stage ) {

 if( v_inner_block_configured[ stage ] && v_inner_solver_configured[ stage ] )
  return;

 auto benders_function = get_benders_function( stage );
 auto inner_block = benders_function->get_inner_block();

 // BlockConfig

 if( ! v_inner_block_configured[ stage ] ) {

  if( ( ! f_inner_block_config ) &&
      ( ! f_inner_block_config_filename.empty() ) ) {
   auto c = Configuration::deserialize( f_inner_block_config_filename );
   if( ! ( f_inner_block_config = dynamic_cast< BlockConfig * >( c ) ) ) {
    delete c;
    throw( std::invalid_argument
           ( "SDDPGreedySolver::configure_inner_block: file " +
             f_inner_block_config_filename + " is not a BlockConfig." ) );
   }
  }

  if( f_inner_block_config ) {
   f_inner_block_config->apply( inner_block );
   v_inner_block_configured[ stage ] = true;
  }
 }

 // BlockSolverConfig

 if( ! v_inner_solver_configured[ stage ] ) {

  if( ( ! f_inner_block_solver_config ) &&
      ( ! f_inner_block_solver_config_filename.empty() ) ) {
   auto c = Configuration::deserialize( f_inner_block_solver_config_filename );
   if( ! ( f_inner_block_solver_config =
           dynamic_cast< BlockSolverConfig * >( c ) ) ) {
    delete c;
    throw( std::invalid_argument
           ( "SDDPGreedySolver::configure_inner_block: file " +
             f_inner_block_solver_config_filename +
             " is not a BlockSolverConfig." ) );
   }
  }

  if( f_inner_block_solver_config ) {
   f_inner_block_solver_config->apply( inner_block );
   v_inner_solver_configured[ stage ] = true;
  }
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::unregister_solver_inner_block( Index stage ) {

 auto benders_function = get_benders_function( stage );
 auto inner_block = benders_function->get_inner_block();

 // BlockSolverConfig

 BlockSolverConfig * inner_block_solver_config = nullptr;

 if( v_BSC.size() > stage && v_BSC[ stage ] )
  inner_block_solver_config = v_BSC[ stage ]->clone();
 else if( f_inner_block_solver_config )
  inner_block_solver_config = f_inner_block_solver_config->clone();

 if( inner_block_solver_config ) {
  inner_block_solver_config->clear();
  inner_block_solver_config->apply( inner_block );
  delete inner_block_solver_config;
 }
 else {
  inner_block->unregister_Solvers();
 }

 v_inner_solver_configured[ stage ] = false;
}

/*--------------------------------------------------------------------------*/

int SDDPGreedySolver::solve( Index stage , bool write_solution ) {

 auto benders_function = get_benders_function( stage );

 auto status = benders_function->compute();

 auto solver = benders_function->get_solver();

 auto solver_has_dual_solution = false;

 if( write_solution ) {
  if( solver->has_var_solution() )
   solver->get_var_solution( f_get_var_solution_config );
  if( auto cda_solver = dynamic_cast< CDASolver * >( solver ) ) {
   if( cda_solver->has_dual_solution() ) {
    cda_solver->get_dual_solution( f_get_dual_solution_config );
    solver_has_dual_solution = true;
   }
  }
 }

 f_has_dual_solution = f_has_dual_solution && solver_has_dual_solution;

 return status;
}

/*--------------------------------------------------------------------------*/

Solver * SDDPGreedySolver::get_sub_solver( Index stage ) const {
 auto benders_function = get_benders_function( stage );
 return benders_function->get_solver();
}

/*--------------------------------------------------------------------------*/

std::vector< double > SDDPGreedySolver::get_solution
( SDDPBlock::Index stage ) const {

 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::get_solution: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 Index solution_size = 0;
 for( Index i = 0 ;
      i < sddp_block->get_num_polyhedral_function_per_sub_block() ; ++i ) {
  solution_size +=
   sddp_block->get_polyhedral_function( stage , i )->get_num_active_var();
 }

 std::vector< double > solution;
 solution.reserve( solution_size );

 for( Index i = 0 ;
      i < sddp_block->get_num_polyhedral_function_per_sub_block() ; ++i ) {
  const auto polyhedral_function =
   sddp_block->get_polyhedral_function( stage , i );

  for( const auto & variable : * polyhedral_function ) {
   solution.push_back
    ( static_cast< const ColVariable & >( variable ).get_value() );
  }
 }
 return solution;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_state( const std::vector< double > & state ,
                                  Index stage ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPGreedySolver::set_state: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 static_cast< SDDPBlock * >( f_Block )->set_state( state , stage , 0 );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::process_outstanding_Modification() {
 while( ! v_mod.empty() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification
  v_mod.pop_front();
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::set_scenario( Index scenario_id , Index stage ) {
 static_cast< SDDPBlock * >( f_Block )->set_scenario( scenario_id , stage );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_simulation_data( Index stage ,
                                              Index scenario_index ) {

 // First, we store the initial state.

 std::vector< double > initial_state;

 if( stage == 0 ) {
  // Store the initial state of the first stage problem.

  if( ! v_initial_state.empty() )
   // Use the initial state given as parameter to the SDDPGreedySolver.
   initial_state = v_initial_state;
  else {
   // Use the initial state given by SDDPBlock.
   initial_state = static_cast< SDDPBlock * >( f_Block )->get_initial_state();
  }

  if( ! initial_state.empty() )
   f_simulation_data.store_initial_state( initial_state , 0 );
 }

 // Store the solution at the given stage as the initial state of the next
 // stage.

 f_simulation_data.store_initial_state( get_solution( stage ) , stage + 1 );

 // Store the subgradient with respect to the final state.
 store_subgradient_final_state( stage );

 if( stage > 0 ) {
  // Store the subgradient with respect to the initial state.
  store_subgradient_initial_state( stage , scenario_index ,
                                   get_solution( stage - 1 ) );
 }

 // Store the objective value disregarding the future value

 f_simulation_data.store_objective_value( get_sub_solution_value( stage ) ,
                                          stage );
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_subgradient_final_state( Index stage ) {

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 // There must be a single PolyhedralFunction per sub-Block.
 assert( sddp_block->get_num_polyhedral_function_per_sub_block() == 1 );

 // The solution is already written into the Block, so there is no need to set
 // the values of the active Variables of the PolyhedralFunction.

 // Find the index of the active row.

 auto polyhedral_function = sddp_block->get_polyhedral_function( stage );
 polyhedral_function->compute();

 std::vector< double > subgradient;

 if( polyhedral_function->has_linearization() ) {
  // The PolyhedralFunction has a linearization and, therefore, there is a row
  // that is active. By optimality conditions, a subgradient of the objective
  // function if given by *minus* the cut of the future cost function that is
  // active at the final state.
  subgradient.resize( polyhedral_function->get_num_active_var() );
  polyhedral_function->get_linearization_coefficients( subgradient.data() );
 }
 else if( polyhedral_function->get_value() ==
          polyhedral_function->get_global_bound() ) {
  // The bound is active, so the subgradient is zero.
  subgradient.assign( polyhedral_function->get_num_active_var() , 0.0 );
 }

 if( ! subgradient.empty() ) {
  // Store the subgradient.
  f_simulation_data.store_subgradient_final_state( std::move( subgradient ) ,
                                                   stage );
 }
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::store_subgradient_initial_state
( Index stage , Index scenario_index ,
  const std::vector< double > & initial_state ) {

 if( stage == 0 )
  return;

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );

 // Get the random cut associated with the previous stage and the given
 // scenario index.

 auto & polyhedral_function = sddp_block->get_random_cut
  ( stage - 1 , scenario_index );

 if( polyhedral_function.get_num_active_var() == 0 ) {
  // The PolyhedralFunction representing the random cut has no
  // Variables. Thus, collect the active Variables of the PolyhedralFunction
  // representing the future cost function at the given stage and make them
  // active Variables of the random cut.
  const auto future_cost_function =
   sddp_block->get_polyhedral_function( stage );
  PolyhedralFunction::VarVector active_variables
   ( future_cost_function->get_num_active_var() );
  for( Index i = 0 ; i < future_cost_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( future_cost_function->get_active_var( i ) );

  polyhedral_function.set_variables( std::move( active_variables ) );
 }

 /* In order to evaluate the PolyhedralFunction at the initial state, we must
  * set the values of its active Variables to be equal to the initial
  * state. Thus, we save the current values of the active Variables in order
  * to undo this modification at the end. */

 std::vector< double > original_variable_values;
 original_variable_values.reserve( polyhedral_function.get_num_active_var() );

 auto initial_state_it = std::cbegin( initial_state );

 for( auto & variable : polyhedral_function ) {
  auto & col_variable = static_cast< ColVariable & >( variable );

  // Save the current value of the active Variable.
  original_variable_values.push_back( col_variable.get_value() );

  // Change the value of the active Variable.
  col_variable.set_value( * ( initial_state_it++ ) );
 }

 // Compute the PolyhedralFunction and obtain the index of the active row.

 polyhedral_function.compute();

 std::vector< double > subgradient;

 if( polyhedral_function.has_linearization() ) {
  // The PolyhedralFunction has a linearization and, therefore, there is a row
  // that is active. By optimality conditions, a subgradient of the objective
  // function if given by *minus* the cut of the future cost function that is
  // active at the final state.
  subgradient.resize( polyhedral_function.get_num_active_var() );
  polyhedral_function.get_linearization_coefficients( subgradient.data() );
 }
 else if( polyhedral_function.get_value() ==
          polyhedral_function.get_global_bound() ) {
  // The bound is active, so the subgradient is zero.
  subgradient.assign( polyhedral_function.get_num_active_var() , 0.0 );
 }

 if( ! subgradient.empty() ) {
  // Store the subgradient.
  f_simulation_data.store_subgradient_initial_state( std::move( subgradient ) ,
                                                     stage );
 }

 // Put back the original values of the active Variables.

 auto original_variable_values_it = std::cbegin( original_variable_values );
 for( auto & variable : polyhedral_function ) {
  auto & col_variable = static_cast< ColVariable & >( variable );
  col_variable.set_value( * ( original_variable_values_it++ ) );
 }

}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::load_cuts( Index stage ) {

 if( f_load_cuts_filename.empty() )
  return;

 const auto sddp_block = static_cast< SDDPBlock * >( f_Block );
 const auto time_horizon = sddp_block->get_time_horizon();

 if( stage >= time_horizon )
  throw( std::logic_error( "SDDPGreedySolver::load_cuts: invalid stage: " +
                           std::to_string( stage ) + "." ) );

 std::ifstream cuts_file( f_load_cuts_filename );

 // Make sure the file is open.
 if( ! cuts_file.is_open() )
  throw( std::runtime_error( "SDDPGreedySolver::load_cuts: It was not possible "
                             "to open the file \"" + f_load_cuts_filename +
                             "\"." ) );

 PolyhedralFunction::MultiVector A;
 PolyhedralFunction::RealVector b;

 std::string line;

 if( cuts_file.good() )
  // Skip the first line containing the header.
  std::getline( cuts_file , line );

 auto polyhedral_function = sddp_block->get_polyhedral_function( stage );
 const auto num_active_var = polyhedral_function->get_num_active_var();

 int line_number = 0;

 // Read the cuts.

 while( std::getline( cuts_file , line ) ) {
  ++line_number;

  std::stringstream line_stream( line );

  // Try to read the stage.
  Index current_stage;
  if( ! ( line_stream >> current_stage ) )
   break;

  if( current_stage >= time_horizon )
   throw( std::logic_error
          ( "SDDPGreedySolver::load_cuts: File \"" + f_load_cuts_filename + "\""
            " contains an invalid stage: " + std::to_string( current_stage ) +
            "." ) );

  if( current_stage != stage )
   continue;

  if( line_stream.peek() != ',' )
   throw( std::logic_error( "SDDPGreedySolver::load_cuts: File \"" +
                            f_load_cuts_filename + "\" has an invalid "
                            "format." ) );
  line_stream.ignore();

  // Read the cut.

  PolyhedralFunction::RealVector a( num_active_var );

  Index i = 0;
  double value;
  while( line_stream >> value ) {
   if( i > num_active_var )
    throw( std::logic_error
           ( "SDDPGreedySolver::load_cuts: File \"" + f_load_cuts_filename +
             "\" contains an invalid cut at line " +
             std::to_string( line_number ) + "." ) );

   if( i < num_active_var )
    a[ i ] = value;
   else
    b.push_back( value );

   ++i;

   if( line_stream.peek() == ',' )
    line_stream.ignore();
  }

  if( i < num_active_var )
   throw( std::logic_error
          ( "SDDPGreedySolver::load_cuts: File \"" + f_load_cuts_filename +
            "\" contains an invalid cut at line " +
            std::to_string( line_number ) + "." ) );

  A.push_back( a );
 }

 cuts_file.close();

 // Now, add the cuts to the PolyhedralFunction. Notice that, even if the
 // SDDPBlock has multiple sub-Blocks per stage, we only add cuts to the first
 // sub-Block of a given stage. This is so because only the first sub-Block
 // associated with each stage is used during the simulation and it may also
 // be the only sub-Block that has been configured.

 // We alsso assume that there is only one PolyhedralFunction per sub-Block.

 assert( sddp_block->get_num_polyhedral_function_per_sub_block() == 1 );

 polyhedral_function->add_rows( std::move( A ) , b );

 // Mark that cuts have been loaded to this stage.

 if( stage >= v_cuts_loaded.size() )
  v_cuts_loaded.resize( get_time_horizon() , false );
 v_cuts_loaded[ stage ] = true;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::output_simulation_data
( const std::string & filename ) const {

 if( filename.empty() )
  return;

 std::ofstream file( filename );

 if( ! file.is_open() )
  throw( std::runtime_error( "SDDPGreedySolver::output_simulation_data: "
                             "it was not possible to open the file \"" +
                             filename + "\"." ) );

 const auto separator = ",";

 // Set the precision of the output

 file << std::setprecision( f_simulation_data_output_precision );

 // Output the time horizon and the number of initial states

 file << get_time_horizon() << separator
      << f_simulation_data.initial_states.size() << std::endl;

 // Initial states

 for( const auto & state : f_simulation_data.initial_states ) {
  file << state.first;
  for( const auto & component : state.second )
   file << separator << component;
  file << std::endl;
 }

 // Subgradients with respect to the initial state

 for( const auto & subgradient : f_simulation_data.subgradients_initial_state ) {
  file << subgradient.first << separator << "I";
  for( const auto & component : subgradient.second )
   file << separator << component;
  file << std::endl;
 }

 // Subgradients with respect to the final state

 for( const auto & subgradient : f_simulation_data.subgradients_final_state ) {
  file << subgradient.first << separator << "F";
  for( const auto & component : subgradient.second )
   file << separator << component;
  file << std::endl;
 }

 // Objective values

 for( const auto & objective_value : f_simulation_data.objective_values ) {
  file << objective_value.first << separator << objective_value.second
       << std::endl;
 }

 // Finally, we output the scenarios

 if( f_output_scenario < 0 ) {
  // Only output the ID of the scenario
  for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {
   auto scenario_id = get_scenario_id( stage );
   file << stage << separator << scenario_id << std::endl;
  }
 }
 else if( f_output_scenario > 0 ) {
  // Output the full scenario
  const auto sddp_block = static_cast< SDDPBlock * >( f_Block );
  const auto & scenario_set = sddp_block->get_scenario_set();

  for( Index stage = 0 ; stage < get_time_horizon() ; ++stage ) {

   auto scenario_id = get_scenario_id( stage );
   auto scenario_begin = scenario_set.sub_scenario_begin( scenario_id , stage );
   auto scenario_end = scenario_set.sub_scenario_end( scenario_id , stage );

   file << stage;

   for( auto scenario_it = scenario_begin ; scenario_it != scenario_end ;
        ++scenario_it )
    file << separator << *scenario_it;
   file << std::endl;
  }
 }

 file.close();
}

/*--------------------------------------------------------------------------*/

Index SDDPGreedySolver::get_scenario_id( Index stage ) const {
 if( stage == 0 ) {
  if( f_first_stage_scenario_id >= 0 )
   return f_first_stage_scenario_id;
  return f_scenario_id;
 }
 else if( f_seed < Inf< Index >() ) { // Random scenarios are being considered
  assert( stage < v_random_scenario_id.size() );
  return v_random_scenario_id[ stage ];
 }
 return f_scenario_id;
}

/*--------------------------------------------------------------------------*/

bool SDDPGreedySolver::should_sample( Index stage ) const {
 if( f_seed == Inf< Index >() )
  return false;

 if( ( f_scenario_sample_frequency > 0 ) &&
     ( stage % f_scenario_sample_frequency == 0 ) )
  return true;

 if( std::find( v_stages_to_sample.begin() , v_stages_to_sample.end() , stage )
     != v_stages_to_sample.end() )
  return true;

 return false;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::sample_scenario( Index stage ) {
 v_random_scenario_id.resize( get_time_horizon() );

 if( stage == 0 ) {
  // The ID of the scenario for the first stage subproblem is not random.
  v_random_scenario_id[ stage ] = get_scenario_id( stage );
  return;
 }

 if( should_sample( stage ) ) {
  using param_type = std::uniform_int_distribution< Index >::param_type;
  const auto num_scenarios =
   static_cast< SDDPBlock * >( f_Block )->get_scenario_set().size();
  v_random_scenario_id[ stage ] = scenario_distribution
   ( random_number_engine , param_type( 0 , num_scenarios - 1 ) );
 }
 else {
  v_random_scenario_id[ stage ] = v_random_scenario_id[ stage - 1 ];
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------------- METHODS of Logger ---------------------------*/
/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::log( double objective_value ,
                                    double future_value ) const {

 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;

 if( objective_value == Inf< double >() )
  *f_log << std::setw( width + 1 ) << "+Inf";
 else if( objective_value == -Inf< double >() )
  *f_log << std::setw( width + 1 ) << "-Inf";
 else {
  int offset = width + 1;
  if( objective_value < 0 )
   offset = width;
  *f_log << std::setw( offset ) << std::setprecision( precision )
         << std::scientific << objective_value;
 }

 *f_log << std::setw( 2 ) << "";

 if( future_value == Inf< double >() )
  *f_log << std::setw( width + 1 ) << "+Inf";
 else if( future_value == -Inf< double >() )
  *f_log << std::setw( width + 1 ) << "-Inf";
 else {
  int offset = width + 1;
  if( future_value < 0 )
   offset = width;
  *f_log << std::setw( offset ) << std::setprecision( precision )
         << std::scientific << future_value;
 }

 *f_log << std::setw( 3 ) << "";

 *f_log << solver->get_subproblem_time() << "   "
        << solver->get_compute_time() << std::endl;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::log( double objective_value ) const {

 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;

 if( objective_value == Inf< double >() )
  *f_log << std::setw( width + 1 ) << "+Inf";
 else if( objective_value == -Inf< double >() )
  *f_log << std::setw( width + 1 ) << "-Inf";
 else {
  int offset = width + 1;
  if( objective_value < 0 )
   offset = width;
  *f_log << std::setw( offset ) << std::setprecision( precision ) <<
   std::scientific << objective_value;
 }

 *f_log << std::setw( width + 6 ) << "-   "
        << solver->get_subproblem_time() << "   "
        << solver->get_compute_time() << std::endl;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::log() const {
 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;
 *f_log << std::setw( width + 3 ) << "-   "
        << std::setw( width + 3 ) << "-   "
        << solver->get_subproblem_time() << "   "
        << solver->get_compute_time() << std::endl;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::log( Index stage , Index scenario ) const {
 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;
 *f_log << std::setw( stage_width ) << stage << std::setw( 3 ) << "";
 *f_log << std::setw( scenario_width - 1 ) << scenario << std::setw( 2 ) << "";
 *f_log << std::flush;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::log_header() const {
 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;

 *f_log << std::setw( stage_width ) << "      "
        << " |" << std::setw( scenario_width ) << "         "
        << " |        Objective value        |               |    Total"
        << std::endl;

 *f_log << std::setw( stage_width ) << " Stage"
        << " |" << std::setw( scenario_width ) << " Scenario"
        << " |    Present    |    Future     |   Time (s)    |   time (s)"
        << std::endl;

 *f_log << std::string( stage_width + scenario_width - 1 , '-' )
        << "-------------------------------------------------------------------"
        << std::endl;
}

/*--------------------------------------------------------------------------*/

void SDDPGreedySolver::Logger::show_status() const {

 if( ( ! f_log ) || ( ! log_verbosity ) )
  return;

 const auto fault_stage = solver->get_fault_stage();

 switch( solver->get_status() ) {

  case( kError ):
   *f_log << "Error while solving the subproblem at stage "
          << fault_stage << std::endl;
   break;

  case( kUnbounded ):
   *f_log << "The subproblem at stage " << fault_stage
          << " is unbounded." << std::endl;
   break;

  case( kInfeasible ):
   *f_log << "The problem is infeasible." << std::endl;
   break;

  case( kStopTime ):
   *f_log << "A feasible solution has been found. The solution process "
          << "of subproblem at stage " << fault_stage
          << " terminated due a time limit." << std::endl;
   break;

  case( kStopIter ):
   *f_log << "A feasible solution has been found. The solution process "
          << "of subproblem at stage " << fault_stage
          << " terminated due to an iteration limit." << std::endl;
   break;

  case( kLowPrecision ):
   *f_log << "A feasible solution has been found." << std::endl;
   break;

  case( kSubproblemInfeasible ):
   *f_log << "The subproblem at stage " << fault_stage
          << " is infeasible." << std::endl;
   break;

  case( kSolutionNotFound ):
   *f_log << "A solution for the subproblem at stage "
          << fault_stage << " has not been found." << std::endl;
   break;
 }

 if( solver->has_var_solution() )
  *f_log << "Objective value: " << solver->get_var_value() << std::endl;

 *f_log << "Total time (s):  " << solver->get_compute_time() << std::endl;
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SDDPGreedySolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
