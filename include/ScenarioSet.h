/*--------------------------------------------------------------------------*/
/*------------------------ File ScenarioSet.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of ScenarioSet, a class to represent a set of scenarios.
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

#ifndef __ScenarioSet
#define __ScenarioSet
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"
#include <boost/multi_array.hpp>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ScenarioSet_CLASSES Classes in ScenarioSet.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ScenarioSet ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// ScenarioSet, a class for representing a set of scenarios
/** ScenarioSet is a class that represents a set of scenarios
 */

class ScenarioSet {

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 using Index = unsigned int;

/**@} ----------------------------------------------------------------------*/
/*--------------- CONSTRUCTING AND DESTRUCTING ScenarioSet -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioSet
 *  @{ */

 /// Constructs a ScenarioSet
 /** Constructs a ScenarioSet
  */

 ScenarioSet() { }

/*--------------------------------------------------------------------------*/

 /// Destructor
 virtual ~ScenarioSet() { }

/*--------------------------------------------------------------------------*/

 /**
  * - The "TimeHorizon" dimension, containing the time horizon.
  *
  * - The "NumberScenarios" dimension specifying the number of
  *   scenarios.
  *
  * - The "ScenarioSize" dimension containing the size of a single
  *   scenario, which spans all stages.
  *
  * - The "SubScenarioSize" variable, of type netCDF::NcUint
  *   and indexed over dimension "TimeHorizon". This dimension is
  *   optional. If it is not provided, then all sub-scenarios are assumed to
  *   have the same size, i.e.,
  *
  *     s_t = ScenarioSize / TimeHorizon
  *
  *   for all t in {0, ..., "TimeHorizon - 1"}, and "ScenarioSize" is a multiple
  *   of "TimeHorizon". If this dimension is provided, then SubScenarioSize[t]
  *   is the size of the sub-scenario associated with stage t, i.e., s_t =
  *   SubScenarioSize[t], for each t in {0, ..., TimeHorizon-1}. In the latter
  *   case, the following must hold:
  *
  *   \f[
  *     \text{ScenarioSize} = \sum_{t = 0}^{\text{TimeHorizon} - 1}
  *                           \text{SubScenarioSize}[t].
  *   \f]
  *
  * - The two-dimensional variable "Scenarios" of type netCDF::NcDouble and
  *   indexed over the dimensions "NumberScenarios" and "ScenarioSize",
  *   containing the scenarios. The i-th row of "Scenarios" contains the i-th
  *   scenario, so that Scenarios[i][j] is the j-th component of the i-th
  *   scenario.
  *
  * - The "NumberRandomDataGroups" dimension containing the number of
  *   groups of related random data within each sub-scenario. This dimension
  *   is optional. If it is not provided, then we assume that there is a
  *   single group of related random data. Also, this dimension is meaningful
  *   only if all sub-scenarios have the same size.
  *
  * - The "SizeRandomDataGroups" variable, of type netCDF::Uint and indexed
  *   over the "NumberRandomDataGroups" dimension, containing the size of each
  *   group of related random data in each sub-scenario. For each i in {0,
  *   ..., "NumberRandomDataGroups - 1"}, NumberRandomDataGroups[i] is the size
  *   of the i-th group of a sub-scenario. This variable is optional. It is
  *   required only if "NumberRandomDataGroups" is provided and
  *   "NumberRandomDataGroups" > 1.
  *
  * @param group A netCDF::NcGroup holding the data describing this
  *        ScenarioSet.
  */

 void deserialize( const netCDF::NcGroup & group ) {

  // TimeHorizon

  ::SMSpp_di_unipi_it::deserialize_dim( group , "TimeHorizon" ,
                                        time_horizon , false );

  // NumberScenarios

  ::SMSpp_di_unipi_it::deserialize_dim( group , "NumberScenarios" ,
                                        num_scenarios , false );

  // ScenarioSize

  ::SMSpp_di_unipi_it::deserialize_dim( group , "ScenarioSize" ,
                                        scenario_size , false );

  // SubScenarioSize

  sub_scenario_start_index.resize( time_horizon + 1 );

  if( ::SMSpp_di_unipi_it::deserialize( group , "SubScenarioSize" ,
                                        time_horizon , sub_scenario_size ,
                                        true , false ) ) {
   // SubScenarioSize was provided

   if( scenario_size !=
       std::accumulate( sub_scenario_size.begin() ,
                        sub_scenario_size.end() ,
                        decltype( sub_scenario_size )::value_type( 0 ) ) )
    throw( std::logic_error( "ScenarioSet::deserialize: The sum of the "
                             "elements in 'SubScenarioSize' must be equal to "
                             "'ScenarioSize'" ) );
  }
  else {
   // SubScenarioSize was not provided

   if( scenario_size % time_horizon != 0 )
    throw( std::logic_error( "ScenarioSet::deserialize: 'SubScenarioSize' was "
                             "not provided. Thus, 'ScenarioSize' must be a "
                             "multiple of 'TimeHorizon'." ) );

   sub_scenario_size.resize( time_horizon , scenario_size / time_horizon );
  }

  // sub_scenario_start_index

  Index next_index = 0;
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   sub_scenario_start_index[ i ] = next_index;
   next_index += sub_scenario_size[ i ];
  }
  sub_scenario_start_index[ time_horizon ] = scenario_size;

  // Scenarios

  deserialize_scenarios( group );

  // NumberRandomDataGroups and SizeRandomDataGroups

  if( std::adjacent_find( sub_scenario_size.begin() , sub_scenario_size.end() ,
                          std::not_equal_to<>() ) != sub_scenario_size.end() ) {
   // Not all sub-scenarios have the same size. In this case, we consider a
   // single random data group.
   num_random_data_groups = 1;
  }
  else {

   // All sub-scenarios have the same size. Thus, we check
   // NumberRandomDataGroups and SizeRandomDataGroups.

   if( ! ::SMSpp_di_unipi_it::deserialize_dim
       ( group , "NumberRandomDataGroups" , num_random_data_groups , true ) ) {
    // NumberRandomDataGroups was not provided. Hence, there must be a single
    // random data group.
    num_random_data_groups = 1;
   }
   else {
    // NumberRandomDataGroups was provided. Now, we check SizeRandomDataGroups.

    if( ::SMSpp_di_unipi_it::deserialize
        ( group , "SizeRandomDataGroups" , num_random_data_groups ,
          size_random_data_groups , true , false ) ) {

     // SizeRandomDataGroups was provided.

     if( ( scenario_size / time_horizon ) != std::accumulate
         ( size_random_data_groups.begin() , size_random_data_groups.end() ,
           decltype( size_random_data_groups )::value_type(0) ) )
      throw( std::logic_error( "ScenarioSet::deserialize: The sum of the "
                               "sizes in 'SizeRandomDataGroups' must be "
                               "equal to 'ScenarioSize' / 'TimeHorizon'." ) );
    }
    else {
     // SizeRandomDataGroups was not provided.
     size_random_data_groups = { ( scenario_size / time_horizon ) };
     if( num_random_data_groups > 1 ) {
      throw( std::logic_error( "ScenarioSet::deserialize: 'NumberRandomData"
                               "Groups' must be provided since "
                               "'NumberRandomDataGroups' > 1." ) );
     }
    }
   }
  }
 }

/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR Saving THE DATA OF THE ScenarioSet ------------*/
/*--------------------------------------------------------------------------*/
/** @name Saving the data of the ScenarioSet
 *  @{ */

 /// serialize an ScenarioSet into a netCDF::NcGroup
 /** Serialize an ScenarioSet into a netCDF::NcGroup with the format explained
  * in the comments of the deserialize() function.
  *
  * @param group The NcGroup in which this ScenarioSet will be serialized.
  */

 void serialize( netCDF::NcGroup & group ) const {

  // TimeHorizon

  auto TimeHorizon_dim = group.getDim( "TimeHorizon" );
  if( TimeHorizon_dim.isNull() )
   TimeHorizon_dim = group.addDim( "TimeHorizon" , get_time_horizon() );

  // NumberScenarios

  auto NumberScenarios_dim = group.addDim( "NumberScenarios" , num_scenarios );

  // ScenarioSize

  auto ScenarioSize_dim = group.addDim( "ScenarioSize" , scenario_size );

  // SubScenarioSize

  ::SMSpp_di_unipi_it::serialize( group , "SubScenarioSize" ,
                                  netCDF::NcUint() , TimeHorizon_dim ,
                                  sub_scenario_size , false );

  // Scenarios

  auto scenarios_var = group.addVar
   ( "Scenarios" , netCDF::NcDouble() ,
     { NumberScenarios_dim , ScenarioSize_dim } );

  for( decltype(scenarios)::size_type i = 0 ; i < scenarios.size() ; ++i )
   scenarios_var.putVar( { i , 0 } , { 1 , scenarios[ i ].size() } ,
                         scenarios[ i ].data() );

  // NumberRandomDataGroups and SizeRandomDataGroups

  auto NumberRandomDataGroups_dim = group.addDim( "NumberRandomDataGroups" ,
                                                  num_random_data_groups );

  if( ! size_random_data_groups.empty() )
   ::SMSpp_di_unipi_it::serialize( group , "SizeRandomDataGroups" ,
                                   netCDF::NcUint() ,
                                   NumberRandomDataGroups_dim ,
                                   size_random_data_groups , false );
 }

/**@} ----------------------------------------------------------------------*/
/*---------------- METHODS FOR MODIFYING THE ScenarioSet -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the ScenarioSet
 *  @{ */

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS DESCRIBING THE BEHAVIOR OF A ScenarioSet ------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of a ScenarioSet
 * @{ */

/*--------------------------------------------------------------------------*/

 /// returns the number of scenarios in this ScenarioSet
 /** This function returns the number of scenarios in this ScenarioSet.
  *
  * @return The number of scenarios. */

 Index size() const {
  return num_scenarios;
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of a scenario in this ScenarioSet
 /** This function returns the size of a scenario in this ScenarioSet. Notice
  * that all scenarios have the same size.
  *
  * @return The size of a scenario. */

 Index get_scenario_size() const {
  return scenario_size;
 }

/*--------------------------------------------------------------------------*/

 /// returns the time horizon
 /** This function returns the time horizon
  *
  * @return The time horizon. */

 Index get_time_horizon() const {
  return time_horizon;
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of each random data group
 /** This function returns a reference to the vector containing the size of
  * each group of related random data. The i-th element in this vector is the
  * size of the i-th group of related random data. If this vector is empty,
  * then all random data groups have the same size, which is the same as the
  * size of the corresponding sub-scenario.
  *
  * @return The size of each group of related random data. */

 const std::vector< Index > & get_size_random_data_groups() const {
  return size_random_data_groups;
 }

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the array containing the scenario \p i
 /** This function returns a (const) pointer to the array containing the data
  * of the \p i-th scenario. The size of this array is given by
  * get_scenario_size().
  *
  * @param i The index of a scenario, which must be between 0 and size() - 1.
  *
  * @return A pointer to the array containing the data of scenario \p i. */

 const double * scenario( Index i ) const {
  if( i >= size() )
   throw( std::invalid_argument
          ( "ScenarioSet::scenario: Invalid scenario index " +
            std::to_string( i ) + ". The total number of scenarios is " +
            std::to_string( size() ) + "." ) );

  return scenarios[ i ].data();
 }

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the array containing the sub-scenario (i, t)
 /** This function returns a (const) pointer to the array containing the data
  * of the sub-scenario associated with time \p t of the \p i-th scenario. The
  * size of this array is given by get_sub_scenario_size( t ).
  *
  * @param i The index of a scenario, which must be between 0 and size() - 1.
  *
  * @param t A time instant, which must be between 0 and
  *          get_time_horizon() - 1.
  *
  * @return A pointer to the array containing the sub-scenario of scenario \p
  *         i associated with time instant \p t. */

 const double * sub_scenario( Index i , Index t ) const {

  if( i >= size() )
   throw( std::invalid_argument
          ( "ScenarioSet::sub_scenario: Invalid scenario index " +
            std::to_string( i ) + ". The total number of scenarios is " +
            std::to_string( size() ) + "."  ) );

  assert( t < get_time_horizon() );
  return( scenarios[ i ].data() + sub_scenario_start_index[ t ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns an iterator to the first element of the sub-scenario (i, t)
 /** This function returns a const iterator to the first element of the
  * sub-scenario associated with time \p t of the \p i-th scenario. The size
  * of this vector is given by get_sub_scenario_size( t ).
  *
  * @param i The index of a scenario, which must be between 0 and size() - 1.
  *
  * @param t A time instant, which must be between 0 and
  *          get_time_horizon() - 1.
  *
  * @return An iterator to the first element of the sub-scenario of scenario
  *         \p i associated with time instant \p t. */

 std::vector< double >::const_iterator
 sub_scenario_begin( Index i , Index t ) const {
  if( i >= size() )
   throw( std::invalid_argument
          ( "ScenarioSet::sub_scenario_begin: Invalid scenario index " +
            std::to_string( i ) + ". The total number of scenarios is " +
            std::to_string( size() ) + "." ) );

  assert( t < get_time_horizon() );
  return std::next( scenarios[ i ].cbegin() , sub_scenario_start_index[ t ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns an iterator to the element following the last in sub-scenario (i,t)
 /** This function returns a const iterator to the element following the last
  * element of the vector containing the data of the sub-scenario associated
  * with time \p t of the \p i-th scenario. The size of this vector is given by
  * get_sub_scenario_size( t ).
  *
  * @param i The index of a scenario, which must be between 0 and size() - 1.
  *
  * @param t A time instant, which must be between 0 and
  *        get_time_horizon() - 1.
  *
  * @return An iterator to the element following the last element of the
  *         vector containing the sub-scenario of scenario \p i associated
  *         with time instant \p t. */

 std::vector< double >::const_iterator
 sub_scenario_end( Index i , Index t ) const {
  if( i >= size() )
   throw( std::invalid_argument
          ( "ScenarioSet::sub_scenario_end: Invalid scenario index " +
            std::to_string( i ) + ". The total number of scenarios is " +
            std::to_string( size() ) + "." ) );

  assert( t < get_time_horizon() );
  return std::next( scenarios[ i ].cbegin() ,
                    sub_scenario_start_index[ t + 1 ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of the sub-scenario associated with time \p t
 /** This function returns the size of the sub-scenario associated with time
  * \p t.
  *
  * @param t A time instant, which must be between 0 and
  *        get_time_horizon() - 1.
  *
  * @return The size of the sub-scenario associated with time instant \p t. */

 auto get_sub_scenario_size( Index t ) const {
  assert( t < get_time_horizon() );
  return sub_scenario_start_index[ t + 1 ] - sub_scenario_start_index[ t ];
 }

/*--------------------------------------------------------------------------*/

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods
    @{ */

/**@} ----------------------------------------------------------------------*/
/*-------------------------- PROTECTED FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected fields
    @{ */

 /// Number of scenarios
 Index num_scenarios;

 /// The size of a scenario (spanning the whole time horizon)
 Index scenario_size;

 /// The size of each sub-scenario
 /** A scenario is divided into sub-scenarios, each sub-scenario being
  * associated with a time instant. This vector stores the size of each
  * sub-scenario. For each t in {0, TimeHorizon -1}, sub_scenario_size[ t ] is
  * the size of the sub-scenario associated with time t.
  */
 std::vector< Index > sub_scenario_size;

 /// Matrix storing the scenarios
 /** The number of rows is the number of scenarios and the number of columns
  * is the size of a scenario.
  */
 std::vector< std::vector< double > > scenarios;

 /// The number of groups of related random data
 Index num_random_data_groups;

 /// The size of each group of related random data
 /** If there are more than one group of related random data, then
  * size_random_data_groups[ i ] is the size of the i-th group of related
  * random data. If this vector is empty, then there is a single group of
  * related random data.
  */
 std::vector< Index > size_random_data_groups;

 Index time_horizon;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRTIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private fields
    @{ */

 std::vector< Index > sub_scenario_start_index;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRTIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private methods
    @{ */

 void deserialize_scenarios( const netCDF::NcGroup & group ) {

  scenarios.resize( num_scenarios , std::vector< double >( scenario_size ) );

  auto scenarios_var = group.getVar( "Scenarios" );

  if( scenarios_var.isNull() )
   throw( std::invalid_argument
          ( "ScenarioSet::deserialize_scenarios: 'Scenarios' "
            "variable has not been provided." ) );

  auto dims = scenarios_var.getDims();

  if( ( dims.size() != 2 ) || ( dims[ 0 ].getSize() != num_scenarios ) ||
      ( dims[ 1 ].getSize() != scenario_size ) )

   throw( std::logic_error
          ( "ScenarioSet::deserialize_scenarios: 'Scenarios' must be a two-"
            "dimensional array whose first and second dimensions have sizes "
            "'NumberScenarios' and 'ScenarioSize', respectively." ) );

  for( decltype(scenarios)::size_type i = 0 ; i < scenarios.size() ; ++i )
   scenarios_var.getVar( { i , 0 } , { 1 , scenarios[ i ].size() } ,
                         scenarios[ i ].data() );
 }

/**@} ----------------------------------------------------------------------*/

};   // end( class ScenarioSet )

/** @} end( group( ScenarioSet_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ScenarioSet.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File ScenarioSet.h ---------------------------*/
/*--------------------------------------------------------------------------*/
