/*--------------------------------------------------------------------------*/
/*--------------------------- File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SDDPBlock class.
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

#include "AbstractPath.h"
#include "BendersBlock.h"
#include "SDDPBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( SDDPBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS of SDDPBlock -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------- CONSTRUCTING AND DESTRUCTING SDDPBlock ------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::deserialize_random_cuts( const std::string & filename ) {

 if( filename.empty() )
  return;

 netCDF::NcFile file( filename.c_str() , netCDF::NcFile::read );

 const auto TimeHorizon = file.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the dimension TimeHorizon "
           "was not provided." ) );

 const auto time_horizon = TimeHorizon.getSize();

 if( time_horizon != get_time_horizon() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the expected TimeHorizon "
           "dimension is " + std::to_string( get_time_horizon() ) +
           ", but " + std::to_string( time_horizon ) + " was given." ) );

 const auto NumberScenarios = file.getDim( "NumberScenarios" );

 if( NumberScenarios.isNull() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the dimension "
           "NumberScenarios was not provided." ) );

 const auto number_scenarios = NumberScenarios.getSize();

 if( number_scenarios != scenario_set.size() )
  throw( std::invalid_argument
         ( "SDDPBlock::deserialize_random_cuts: the expected NumberScenarios "
           "dimension is " + std::to_string( scenario_set.size() ) +
           ", but " + std::to_string( number_scenarios ) + " was given." ) );

 // Possibly clear the previous random cuts
 random_cuts.resize( boost::extents[ 0 ][ 0 ] );

 // Create the random cuts
 random_cuts.resize( boost::extents[ time_horizon ][ number_scenarios ] );

 for( Index t = 0 ; t < time_horizon ; ++t ) {

  // Collect the active Variables of the PolyhedralFunction at stage t
  const auto polyhedral_function = get_polyhedral_function( t );
  PolyhedralFunction::VarVector active_variables
   ( polyhedral_function->get_num_active_var() );
  for( Index i = 0 ; i < polyhedral_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( polyhedral_function->get_active_var( i ) );

  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   // Set the active Variables of the PolyhedralFunction
   auto variables = active_variables;
   random_cuts[ t ][ s ].set_variables( std::move( variables ) );

   // Deserialize the PolyhedralFunction (if provided)
   auto group_name = "PolyhedralFunction_" +
    std::to_string( t ) + "_" + std::to_string( s );
   auto group = file.getGroup( group_name );
   if( group.isNull() )
    continue;
   random_cuts[ t ][ s ].deserialize( group );
  }
 }
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_Modification( sp_Mod mod , Observer::ChnlName chnl ) {
 // TODO
 if( anyone_there() )
  Block::add_Modification( std::make_shared< NBModification >( this ) , chnl );
}

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE SDDPBlock --------------*/
/*--------------------------------------------------------------------------*/

int SDDPBlock::get_objective_sense() const {
 try {
  auto sub_Block = get_sub_Block( 0 );
  if( sub_Block )
   return sub_Block->get_objective_sense();
 }
 catch( ... ) {}
 return Objective::eUndef;
}

/*--------------------------------------------------------------------------*/

StochasticBlock * SDDPBlock::get_sub_Block
( Index stage , Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid stage " +
                                std::to_string( stage ) ) );
 if( sub_block_index >= num_sub_blocks_per_stage )
  throw( std::invalid_argument( "SDDPBlock::get_sub_Block: invalid sub-Block in"
                                "dex " + std::to_string( sub_block_index ) ) );
 const auto index = stage * num_sub_blocks_per_stage + sub_block_index;
 return static_cast< StochasticBlock * >( v_Block[ index ] );
}

/*--------------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF AN SDDPBlock -------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::add_cuts( PolyhedralFunction::MultiVector && A ,
                          PolyhedralFunction::RealVector && b , Index stage ,
                          Index sub_block_index ,
                          Index number_cuts_to_keep ) {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::add_cuts: invalid stage index: " +
                                std::to_string( stage ) ) );

 auto polyhedral_function =
  get_polyhedral_function( stage , 0 , sub_block_index );

 const auto num_rows = polyhedral_function->get_nrows();

 // Possibly remove the last (num_rows - number_cuts_to_keep) cuts
 if( number_cuts_to_keep < num_rows )
  polyhedral_function->delete_rows( Range( number_cuts_to_keep , num_rows) );

 // Add the given cuts
 polyhedral_function->add_rows( std::move( A ) , b );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::store_random_cut( std::vector< double > && coefficients ,
                                  double alpha , Index stage ,
                                  Index scenario_index ) {

 assert( stage < get_time_horizon() );
 assert( scenario_index < scenario_set.size() );

 if( ! f_random_cuts_initialized ) {
  // Create the PolyhedralFunctions that will store the random cuts.

  // Ensure the random cuts are initialized by only one thread.
#pragma omp critical (SDDPBlock_random_cut)
  {
   if( ! f_random_cuts_initialized ) {
    initialize_random_cuts();
    f_random_cuts_initialized = true;
   }
  }
 }

 // Store the given random cut.

 random_cuts[ stage ][ scenario_index ].add_row( std::move( coefficients ) ,
                                                 alpha );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::initialize_random_cuts() {

 const auto time_horizon = get_time_horizon();
 const auto number_scenarios = scenario_set.size();
 random_cuts.resize( boost::extents[ 0 ][ 0 ] );
 random_cuts.resize( boost::extents[ time_horizon ][ number_scenarios ] );

 for( Index t = 0 ; t < time_horizon ; ++t ) {

  const auto polyhedral_function = get_polyhedral_function( t );
  PolyhedralFunction::VarVector active_variables
   ( polyhedral_function->get_num_active_var() );
  for( Index i = 0 ; i < polyhedral_function->get_num_active_var() ; ++i )
   active_variables[ i ] = static_cast< ColVariable * >
    ( polyhedral_function->get_active_var( i ) );

  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   auto variables = active_variables;
   random_cuts[ t ][ s ].set_variables( std::move( variables ) );
   random_cuts[ t ][ s ].set_is_convex( polyhedral_function->is_convex() );
  }
 }
}

/*--------------------------------------------------------------------------*/

double SDDPBlock::get_future_cost( Index stage , Index sub_block_index ) const {
 if( stage >= get_time_horizon() )
  throw( std::invalid_argument( "SDDPBlock::get_future_cost: invalid "
                                "stage index: " + std::to_string( stage ) ) );

 auto function = get_polyhedral_function( stage , 0 , sub_block_index );
 if( ! function )
  // If this SDDPBlock has no PolyhedralFunction, the future cost must be 0
  return 0;
 function->compute();
 return function->get_value();
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const Eigen::ArrayXd & values , Index stage ,
                           Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_state( const std::vector< double > & values , Index stage ,
                           Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 benders_block->set_variable_values( values );
}

/*--------------------------------------------------------------------------*/

std::vector< double > SDDPBlock::get_state( Index stage ,
                                            Index sub_block_index ) const {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );
 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );
 return benders_block->get_variable_values();
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::set_admissible_state( Index stage , Index sub_block_index ) {
 assert( stage < get_time_horizon() );
 assert( sub_block_index < get_num_sub_blocks_per_stage() );

 auto benders_block = static_cast< BendersBlock * >
  ( get_sub_Block( stage , sub_block_index )->get_nested_Block( 0 ) );

 auto admissible_state = get_admissible_state( stage );
 benders_block->set_variable_values( admissible_state );
}

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR PRINTING & SAVING THE SDDPBlock ---------------*/
/*--------------------------------------------------------------------------*/

void SDDPBlock::print( std::ostream & output , char vlvl ) const
{
 output << std::endl << "SDDPBlock with ";

 if( v_Block.empty() )
  output << "no inner Block";
 else
  output << v_Block.size() << " sub-Blocks" << std::endl;
 }

/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 // type

 group.putAtt( "type" , "SDDPBlock" );

 // TimeHorizon

 const auto time_horizon = get_time_horizon();
 auto TimeHorizon_dim = group.addDim( "TimeHorizon" , time_horizon );

 // NumSubBlocksPerStage

 group.addDim( "NumSubBlocksPerStage" , num_sub_blocks_per_stage );

 // StochasticBlock_i

 for( Index i = 0 ; i < get_time_horizon() ; ++i ) {
  auto sub_group = group.addGroup( "StochasticBlock_" +
                                   std::to_string( i ) );
  get_sub_Block( i )->serialize( sub_group );
 }

 // AbstractPaths to PolyhedralFunctions

 std::vector< AbstractPath > paths;
 paths.reserve( get_time_horizon() * num_polyhedral_per_sub_block );

 for( Index t = 0 ; t < get_time_horizon() ; ++t ) {
  for( Index i = 0 ; i < num_polyhedral_per_sub_block ; ++i ) {
   auto reference_block = get_sub_Block( t )->get_nested_Block( 0 );
   assert( reference_block );
   paths.emplace_back( get_polyhedral_function( t , i ) , reference_block );
  }
 }

 AbstractPath::serialize( paths , group );

 if( num_polyhedral_per_sub_block != 1 )
  group.addDim( "NumPolyhedralFunctionsPerSubBlock" ,
                num_polyhedral_per_sub_block );

 // Scenarios

 scenario_set.serialize( group );

 // Initial state

 auto InitialStateSize = group.addDim( "InitialStateSize" ,
                                       initial_state.size() );

 ::SMSpp_di_unipi_it::serialize( group , "InitialState" ,
                                 netCDF::NcDouble() , InitialStateSize ,
                                 initial_state , false );

 // StateSize

 std::vector< Index > state_size( time_horizon );
 for( Index t = 0 ; t < time_horizon - 1 ; ++t )
  state_size[ t ] = admissible_state_begin[ t+1 ] - admissible_state_begin[ t ];
 if( time_horizon > 0 )
  state_size.back() = admissible_states.size() - admissible_state_begin.back();

 ::SMSpp_di_unipi_it::serialize( group , "StateSize" , netCDF::NcUint64() ,
                                 TimeHorizon_dim , state_size , false );

 // AdmissibleState

 auto AdmissibleState_dim = group.addDim( "AdmissibleState_dim" ,
                                          admissible_states.size() );

 ::SMSpp_di_unipi_it::serialize( group , "AdmissibleState" ,
                                 netCDF::NcDouble() , AdmissibleState_dim ,
                                 admissible_states , false );
}

/*--------------------------------------------------------------------------*/

void SDDPBlock::serialize_random_cuts( const std::string & filename ) const {

 if( filename.empty() || ( random_cuts.num_elements() == 0 ) )
  return;

 netCDF::NcFile file( filename , netCDF::NcFile::replace );

 const auto time_horizon = get_time_horizon();
 file.addDim( "TimeHorizon" , time_horizon );

 const auto number_scenarios = random_cuts[ 0 ].size();
 file.addDim( "NumberScenarios" , number_scenarios );

 for( Index t = 0 ; t < time_horizon ; ++t ) {
  for( Index s = 0 ; s < number_scenarios ; ++s ) {
   auto group_name = "PolyhedralFunction_" +
    std::to_string( t ) + "_" + std::to_string( s );
   auto group = file.addGroup( group_name );
   random_cuts[ t ][ s ].serialize( group );
  }
 }
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File SDDPBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
