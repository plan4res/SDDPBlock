/*--------------------------------------------------------------------------*/
/*---------------------- File ParallelSDDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ParallelSDDPSolver class.
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
/*------------------------------- INCLUDES ---------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ParallelSDDPSolver.h"

#include <Eigen/Core>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( ParallelSDDPSolver );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS of ParallelSDDPOptimizer --------------------*/
/*--------------------------------------------------------------------------*/

Eigen::ArrayXd ParallelSDDPSolver::ParallelSDDPOptimizer::oneStepBackward
( const StOpt::SDDPCutOptBase & sddp_cut ,
  const std::tuple< std::shared_ptr< Eigen::ArrayXd > , int , int > & state ,
  const Eigen::ArrayXd & particle , const int & simulation_id ) const {

 const auto current_stage = get_current_backward_stage();
 const auto scenario_index = get_backward_scenario_index( simulation_id );
 const auto sub_block_index = lock( current_stage , simulation_id , true );

 const auto scenario_must_be_set =
  ! is_scenario_set( current_stage , sub_block_index , scenario_index );

 auto backward = SDDPSolver::SDDPOptimizer::oneStepBackward
  ( sddp_cut , state , particle , simulation_id , scenario_index ,
    scenario_must_be_set , sub_block_index );

 unlock( current_stage , sub_block_index , scenario_index , true );

 return backward;
}

/*--------------------------------------------------------------------------*/

double ParallelSDDPSolver::ParallelSDDPOptimizer::oneStepForward
( const Eigen::ArrayXd & particle , Eigen::ArrayXd & state ,
  Eigen::ArrayXd & state_to_store , const StOpt::SDDPCutOptBase & sddp_cut ,
  const int & simulation_id ) const {

 const auto current_stage = get_current_forward_stage();
 const auto scenario_index = get_forward_scenario_index( simulation_id );
 const auto sub_block_index = lock( current_stage , simulation_id , false );

 const auto scenario_must_be_set =
  ! is_scenario_set( current_stage , sub_block_index , scenario_index );

 auto forward = SDDPSolver::SDDPOptimizer::oneStepForward
  ( particle , state , state_to_store , sddp_cut , simulation_id ,
    scenario_index , scenario_must_be_set , sub_block_index );

 unlock( current_stage , sub_block_index , scenario_index , true );

 return forward;
}

/*--------------------------------------------------------------------------*/

bool ParallelSDDPSolver::ParallelSDDPOptimizer::is_locked
( Index sub_block_index ) const {
 if( locked.empty() )
  return false;
 return locked[ sub_block_index ];
}

/*--------------------------------------------------------------------------*/

void ParallelSDDPSolver::ParallelSDDPOptimizer::lock_sub_block
( Index sub_block_index ) const {
 if( locked.empty() ) {
  locked.resize( get_num_sub_blocks_per_stage() , false );
 }
 locked[ sub_block_index ] = true;
}

/*--------------------------------------------------------------------------*/

Index ParallelSDDPSolver::ParallelSDDPOptimizer::get_sub_block_with_scenario
( Index stage , Index scenario_index ) const {
 if( scenario_currently_set.empty() ||
     scenario_currently_set[ stage ].empty() )
  return Inf< Index >();
 auto iter = std::find( scenario_currently_set[ stage ].cbegin() ,
                        scenario_currently_set[ stage ].cend() ,
                        scenario_index );
 if( iter != scenario_currently_set[ stage ].end() )
  return std::distance( scenario_currently_set[ stage ].cbegin() , iter );
 return Inf< Index >();
}

/*--------------------------------------------------------------------------*/

void ParallelSDDPSolver::ParallelSDDPOptimizer::prepare_new_stage
( Index stage , bool backward ) const {

 const auto num_sub_blocks = get_num_sub_blocks_per_stage();

 // Resize some containers if necessary

 last_sub_block_solved.resize( sddp_solver->get_time_horizon() );

 if( scenario_currently_set.empty() )
  scenario_currently_set.resize( sddp_solver->get_time_horizon() ,
                                 std::vector< Index >() );

 if( scenario_currently_set[ stage ].empty() )
  scenario_currently_set[ stage ].resize( num_sub_blocks , Inf< Index >() );

 // Fill the preferred, reserved, and non-reserved vectors

 const auto num_simulations = get_number_simulations( backward );

 preferred_block.resize( num_simulations );

 reserved_blocks.clear();
 reserved_blocks.reserve( num_simulations );

 non_reserved_block.assign( num_sub_blocks , true );

 auto i = decltype( num_simulations ){ 0 };
 for( auto & block_index : preferred_block ) {
  const auto scenario_index = get_scenario_index( i++ , backward );
  block_index = get_sub_block_with_scenario( stage , scenario_index );
  if( block_index < Inf< Index >() ) {
   reserved_blocks.push_back( block_index );
   non_reserved_block[ block_index ] = false;
  }
 }
}

/*--------------------------------------------------------------------------*/

void ParallelSDDPSolver::ParallelSDDPOptimizer::synchronize_cuts() const {
 if( mesh_provided() )
  return; // A mesh was provided. Nothing needs to be synchronized.

 const auto num_sub_block_per_stage = get_num_sub_blocks_per_stage();

 for( Index stage = 0 ; stage < sddp_solver->get_time_horizon() - 1 ; ++stage ) {
  // The sub-Blocks at the last stage do not need to be synchronized, as their
  // cuts are added only once in the beginning of SDDPSolver::compute().

  // This function is guaranteed to have all cuts
  const auto complete_function =
   get_polyhedral_function( stage , 0 , last_sub_block_solved[ stage ] );
  const auto num_complete_rows = complete_function->get_nrows();

  const auto & complete_A = complete_function->get_A();
  const auto & complete_b = complete_function->get_b();

  // Updates all functions that are out of sync
  for( Index i = 0 ; i < num_sub_block_per_stage ; ++i ) {

   auto polyhedral_function = get_polyhedral_function( stage , 0 , i );
   auto num_rows = polyhedral_function->get_nrows();

   if( num_rows < num_complete_rows ) {
    // This PolyhedralFunction is out of sync. The last cuts of the "complete"
    // PolyhedralFunction are added to it.

    const auto num_new_rows = num_complete_rows - num_rows;
    PolyhedralFunction::MultiVector A( num_new_rows );
    PolyhedralFunction::RealVector b( num_new_rows );

    for( Index k = 0 ; k < num_new_rows ; ++k ) {
     b[ k ] = complete_b[ num_rows + k ];
     A[ k ] = complete_A[ num_rows + k ];
    }

    polyhedral_function->add_rows( std::move( A ) , std::move( b ) );
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

Index ParallelSDDPSolver::ParallelSDDPOptimizer::lock
( Index stage , Index simulation_id , bool backward ) const {

 auto sub_block_index = Inf< Index >();

 while( true ) {

#pragma omp critical (ParallelSDDPSolver)
  {

   if( backward ) {
    // This is a backard pass. The sub-Blocks will possibly become out of sync.

    if( cuts_synchronized ) {
     // This means that this is the first subproblem to be solved in the
     // backward pass after a forward pass has been performed.
     assert( stage == sddp_solver->get_time_horizon() - 1 );

     // Thus, the iteration number is updated.
     current_iteration++;

     // Indicate that the sub-Blocks are possibly out of sync.
     cuts_synchronized = false;
    }
   }
   else if( ! cuts_synchronized ) {
    // This means that this is the first subproblem to be solved in the
    // forward pass.
    assert( stage == 0 );

    // Synchronize cuts if necessary.
    synchronize_cuts();
    cuts_synchronized = true;

    // If it is time to output the cuts, do it.
    if( ( get_output_frequency() > 0 ) &&
        ( current_iteration % get_output_frequency() == 0 ) ) {
     sddp_solver->file_output();
    }
   }

   // Prepare for a new stage.
   if( new_stage ) {
    prepare_new_stage( stage , backward );
    new_stage = false;
   }

   // Try to find an unlocked sub-Block
   sub_block_index = find_available_sub_block( simulation_id , stage );

   if( sub_block_index < Inf< Index >() ) {
    // A sub-Block is available. Lock it.
    lock_sub_block( sub_block_index );
   }

  } // end omp critical (ParallelSDDPSolver)

  if( sub_block_index < Inf< Index >() )
   // An unlocked sub-Block has been found. Return its index.
   return sub_block_index;
  else
   // No sub-Block is available. Wait.
   std::this_thread::sleep_for
    ( std::chrono::duration< double >( waiting_time ) );
 }
}

/*--------------------------------------------------------------------------*/

void ParallelSDDPSolver::ParallelSDDPOptimizer::unlock
( Index stage , Index sub_block_index , Index scenario_index ,
  bool scenario_was_set ) const {

#pragma omp critical (ParallelSDDPSolver)
 {
  if( scenario_was_set )
   // If a scenario was set, mark it
   scenario_currently_set[ stage ][ sub_block_index ] = scenario_index;

  // Keep track of the last sub-Block solved
  last_sub_block_solved[ stage ] = sub_block_index;

  // Unlock the sub-Block
  locked[ sub_block_index ] = false;

 } // end omp critical (ParallelSDDPSolver)
}

/*--------------------------------------------------------------------------*/

bool ParallelSDDPSolver::ParallelSDDPOptimizer::is_scenario_set
( Index stage , Index sub_block_index , Index scenario_index ) const {
 if( scenario_currently_set.empty() ||
     scenario_currently_set[ stage ].empty() )
  return false;
 return scenario_currently_set[ stage ][ sub_block_index ] == scenario_index;
}

/*--------------------------------------------------------------------------*/

Index ParallelSDDPSolver::ParallelSDDPOptimizer::find_available_sub_block
( Index simulation_id , Index stage ) const {

 // Check if the preferred sub-Block for the given simulation is available
 const auto preferred_block_index = preferred_block[ simulation_id ];
 if( preferred_block_index < Inf< Index >() &&
     ( ! is_locked( preferred_block_index ) ) ) {
  return preferred_block_index;
 }

 // Search for an unlocked sub-Block among the non-reserved sub-Blocks
 Index non_reserved_index = Inf< Index >();
 for( Index i = 0 ; i < non_reserved_block.size() ; ++i )
  if( non_reserved_block[ i ] && ( ! is_locked( i ) ) ) {
   if( scenario_currently_set[ stage ][ i ] == Inf< Index >() )
    // Give preference to sub-Blocks without a scenario
    return i;
   if( non_reserved_index == Inf< Index >() )
    // Keep the index of the first non-reserved sub-Block
    non_reserved_index = i;
  }

 if( non_reserved_index < Inf< Index >() )
  // No sub-Block without a scenario is available, but a non-reserved
  // sub-Block is. So, that is the one we select.
  return non_reserved_index;

 // Finally, search for an unlocked sub-Block in the list of reserved ones
 for( auto index : reserved_blocks )
  if( ! is_locked( index ) )
   return index;

 // All sub-Blocks are locked
 return Inf< Index >();
}

/*--------------------------------------------------------------------------*/
/*------------------- End File ParallelSDDPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
