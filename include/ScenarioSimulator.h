/*--------------------------------------------------------------------------*/
/*--------------------- File ScenarioSimulator.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file of ScenarioSimulator, a class that derives from
 * StOpt::SimulatorSDDPBase and serves as a simple simulator for a set of
 * scenarios.
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

#ifndef __ScenarioSimulator
#define __ScenarioSimulator
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioSet.h"
#include "StOpt/sddp/SimulatorSDDPBase.h"

#include <algorithm>
#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include <random>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ScenarioSimulator_CLASSES Classes in ScenarioSimulator.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioSimulator ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// ScenarioSimulator, a SimulatorSDDPBase for a set of scenarios
/** ScenarioSimulator is a class that derives from StOpt::SimulatorSDDPBase
 * and serves as a simple simulator for a set of scenarios.
 */

class ScenarioSimulator : public StOpt::SimulatorSDDPBase {

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 using Index = unsigned int;

/**@} ----------------------------------------------------------------------*/
/*------------ CONSTRUCTING AND DESTRUCTING ScenarioSimulator --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioSimulator
 *  @{ */

 /// Constructs an empty ScenarioSimulator
 /** Constructs an empty ScenarioSimulator.
  *
  * @param backward It indicates whether this is a simulator for the SDDP
  *        backward sweep. */

 ScenarioSimulator( const bool backward = true ) :
  SimulatorSDDPBase() , backward_simulator( backward ) {
  random_number_engine.seed( initial_seed );
 }

/*--------------------------------------------------------------------------*/

 /// Constructs a ScenarioSimulator for the given set of scenarios
 /** Constructs a ScenarioSimulator for the given set of scenarios.
  *
  * @param scenario_set A set of scenarios.
  *
  * @param backward It indicates whether this is a simulator for the SDDP
  *        backward sweep. */

 ScenarioSimulator( const ScenarioSet & scenario_set , const bool backward ) :
  SimulatorSDDPBase() , backward_simulator( backward ) {
  const auto number_scenarios = scenario_set.size();

  if( backward_simulator )
   set_number_simulations( number_scenarios );
  else
   set_number_simulations( std::min( 3u , number_scenarios ) );

  set_scenarios( scenario_set );
  random_number_engine.seed( initial_seed );
 }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~ScenarioSimulator() {}

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR MODIFYING THE ScenarioSimulator ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the ScenarioSimulator
 *  @{ */

 /// defines the number of simulations to be produced
 void set_number_simulations( int n ) {
  number_simulations = n;
  indices_selected_particles.resize( number_simulations );
 }

/*--------------------------------------------------------------------------*/

 /// indicates if the sampling of scenarios must be with or without replacement
 /** This function determines if the sampling of scenarios must be done with
  * or without replacement.
  *
  * @param with_replacement If true, the sampling is done with
  *        replacement. Otherwise, the sampling is done without
  *        replacement. */

 void set_sampling_replacement( bool with_replacement ) {
  sampling_with_replacement = with_replacement;
 }

/*--------------------------------------------------------------------------*/

 /// indicates if the sampling must be done at at time instant
 /** This function determines if the sampling of scenarios must be done at
  * each time instant as opposed to only at time 0 (for a forward simulator)
  * or time get_number_dates() - 1 (for a backard simulator).
  *
  * @param resampling If true, the sampling is done at each time step within
  *        updateDateIndex(). If false, the sampling is performed only at
  *        extreme time steps. */

 void set_intermediate_sampling( bool intermediate_sampling ) {
  this->intermediate_sampling = intermediate_sampling;
 }

/*--------------------------------------------------------------------------*/

 /// defines the set of scenarios
 void set_scenarios( const ScenarioSet & scenario_set ) {

  // Construct the particles

  const auto num_scenarios = scenario_set.size();

  distribution = std::uniform_int_distribution< Index >
   ( 0 , num_scenarios - 1 );

  const auto time_horizon = scenario_set.get_time_horizon();

  all_particles.resize( time_horizon );

  const auto & size_random_data_groups =
   scenario_set.get_size_random_data_groups();

  const auto particle_length = size_random_data_groups.empty() ? 1 :
   size_random_data_groups.size();

  for( Index t = 0 ; t < time_horizon ; ++t ) {

   all_particles[ t ].resize( particle_length , num_scenarios );

   for( Index i = 0 ; i < num_scenarios ; ++i ) {

    auto sub_scenario_begin = scenario_set.sub_scenario_begin( i , t );
    auto sub_scenario_end = scenario_set.sub_scenario_end( i , t );

    if( particle_length == 1 ) {
     all_particles[ t ]( 0 , i ) =
      std::accumulate( sub_scenario_begin , sub_scenario_end , double( 0.0 ) ) /
      std::distance( sub_scenario_begin , sub_scenario_end );
    }
    else {
     Index start = 0;
     for( Index k = 0 ; k < size_random_data_groups.size() ; ++k ) {
      all_particles[ t ]( k , i ) =
       std::accumulate( sub_scenario_begin + start ,
                        sub_scenario_begin + start +
                        size_random_data_groups[ k ] ,
                        double( 0.0 ) ) / size_random_data_groups[ k ];
      start += size_random_data_groups[ k ];
     }
    }
   }
  }
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS DESCRIBING THE BEHAVIOR OF A ScenarioSimulator ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of a ScenarioSimulator
 * @{ */

 /// returns the number of particles (used in regression part)
 int getNbSimul() const override {
  return number_simulations;
 }

/*--------------------------------------------------------------------------*/

 /// returns the number of samples
 /** This function returns the number of samples. For the ScenarioSimulator,
  * the number of samples is fixed to 1. */

 int getNbSample() const override {
  return 1;
 }

/*--------------------------------------------------------------------------*/

 /// update the simulator for the given date
 /** Update the simulator for the date associated with the given index.
  *
  * @param date_index The index of the date to be considered. */

 void updateDateIndex( const int & date_index ) override {
  current_date_index = date_index;

  if( backward_simulator ) {
   if( ( current_date_index == int( get_number_dates() ) - 1 ) ||
       intermediate_sampling )
    sample();
  }
  else {
   if( ( current_date_index == 0 ) || intermediate_sampling )
    sample();
  }
 }

/*--------------------------------------------------------------------------*/

 /// returns the particle associated with the given index
 /** This function returns the particle associated with the given \p index.
  *
  * @param index The index of the particle to be returned.
  *
  * @return The particle associated with the given \p index. */

 Eigen::VectorXd getOneParticle( const int & index ) const override {
  return all_particles[ current_date_index ].col
   ( indices_selected_particles[ index ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns all particles associated with the current date
 /** This function returns the matrix containing all particles associated with
  * the current date. The number of particles is equal to the number of
  * columns in this matrix and each column of this matrix is a particle.
  *
  * @return The matrix containing all particles. */

 Eigen::MatrixXd getParticles() const override {

  assert( decltype( indices_selected_particles )::size_type( number_simulations )
          == indices_selected_particles.size() );

  Eigen::MatrixXd particles( all_particles[ current_date_index ].rows() ,
                             number_simulations );

  for( int i = 0 ; i < number_simulations ; ++i ) {
   particles.col( i ) = all_particles[ current_date_index ].col
    ( indices_selected_particles[ i ] );
  }

  return particles;
 }

/*--------------------------------------------------------------------------*/

 /// reset the simulator to use it again for another SDDP sweep
 /** This function resets the simulator to use it again for another SDDP
  * sweep. For a backward simulator, it resets the random seed and sets the
  * current date index to get_number_dates() - 1. For a forward simulator, it
  * sets the current date index to 0. */

 virtual void resetTime() override {
  if( backward_simulator ) {
   random_number_engine.seed( initial_seed );
   updateDateIndex( get_number_dates() - 1 );
  }
  else {
   updateDateIndex( 0 );
  }
 }

/*--------------------------------------------------------------------------*/

 /** This function simply updates the the number of simulations and calls
  * resetTime(). This function should only be called by a forward simulator.
  *
  * @param number_simulations The new number of simulations. */

 virtual void updateSimulationNumberAndResetTime
 ( const int & number_simulations ) override {
  assert( ! backward_simulator );
  set_number_simulations( number_simulations );
  resetTime();
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE ScenarioSimulator ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ScenarioSimulator
    @{ */

 /// returns the index of the scenario associated with the given simulation id
 Index get_scenario_index( Index simulation_id ) const {
  assert( simulation_id < indices_selected_particles.size() );
  return indices_selected_particles[ simulation_id ];
 }

/*--------------------------------------------------------------------------*/

 /// returns the size of a particle
 Index get_particle_length() const {
  if( all_particles.empty() )
   return 0;
  return all_particles.front().rows();
 }

/*--------------------------------------------------------------------------*/

 /// returns the total number of scenarios available
 Index get_number_scenarios() const {
  if( all_particles.empty() )
   return 0;
  return all_particles.front().cols();
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods
    @{ */

 /// returns the number of dates, which is equal to the time horizon
 Index get_number_dates() const {
  return all_particles.size();
 }

/*--------------------------------------------------------------------------*/

 /// sample the particles
 void sample() {
  indices_selected_particles.resize( number_simulations );

  if( sampling_with_replacement ) {
   std::generate( indices_selected_particles.begin() ,
                  indices_selected_particles.end() ,
                  [ this ]() { return this->distribution
                    ( this->random_number_engine ); } );
   std::sort( indices_selected_particles.begin() ,
              indices_selected_particles.end() );
  }
  else {
   const auto number_scenarios = get_number_scenarios();

   assert( decltype( number_scenarios )( number_simulations ) <=
           number_scenarios );

   if( decltype( number_scenarios )( number_simulations ) ==
       number_scenarios ) {
    std::iota( indices_selected_particles.begin() ,
               indices_selected_particles.end() , 0 );
   }
   else {

    if( indices_all_particles.size() != number_scenarios ) {
     indices_all_particles.resize( number_scenarios );
     std::iota( indices_all_particles.begin() ,
                indices_all_particles.end() , 0 );
    }

    std::sample( indices_all_particles.begin() , indices_all_particles.end() ,
                 indices_selected_particles.begin() ,
                 number_simulations , random_number_engine );
   }
  }
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- PROTECTED FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// The vector containing the matrices representing the particles
 /** This is a vector containing the particles for each time step. The t-th
  * element in this vector is a matrix containing the particles associated
  * with the time step t. For each of these matrices, each column is a
  * particle. So, the number of particles is equal to the number of columns in
  * a matrix and the dimension of each particle is equal to the number of rows
  * in a matrix. */

 std::vector< Eigen::MatrixXd > all_particles;

/*--------------------------------------------------------------------------*/

 /// The vector containing the indices of the particles
 /** This is a vector containing the indices of all particles, i.e., contains
  * the set {0, 1, ..., number_scenarios - 1}. */

 std::vector< Index > indices_all_particles;

/*--------------------------------------------------------------------------*/

 /// The vector containing the indices of the selected particles

 std::vector< Index > indices_selected_particles;

/*--------------------------------------------------------------------------*/

 /// The index associated with the current date
 int current_date_index = 0;

 /// The number of simulations to be produced
 int number_simulations = 0;

 /// Indicates whether this is a simulator for the SDDP backward sweep
 bool backward_simulator = true;

 /// Random number engine to select the simulations
 std::mt19937 random_number_engine;

 /// Initial seed for the random number engine
 unsigned int initial_seed = 93645u;

 /// Indicates whether the sampling must be performed with replacement
 bool sampling_with_replacement = false;

 /// Indicates whether a new sample must be drawn at each time instant
 bool intermediate_sampling = false;

 /// Distribution for selecting the particles
 std::uniform_int_distribution< Index > distribution;

};   // end( class ScenarioSimulator )

/** @} end( group( ScenarioSimulator_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ScenarioSimulator.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ScenarioSimulator.h ------------------------*/
/*--------------------------------------------------------------------------*/
