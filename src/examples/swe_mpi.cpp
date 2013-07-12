/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */

#include <mpi.h>

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "blocks/SWE_WavePropagationBlock.hh"

#include "writer/VtkWriter.hh"

#include "scenarios/SWE_simple_scenarios.hh"

#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

// Exchanges the left and right ghost layers.
void exchangeGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                          const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                          const int i_boundarySize);

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  int l_myRank, l_numProc;
  
  // Initialise MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &l_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &l_numProc);
  
  // Tell the logger our rank
  tools::Logger::logger.setProcessRank(l_myRank);
  
  /**
   * Initialization
   */
  // Parse command line parameters
  tools::Args args;
  args.addOption("grid-size-x", 'x', "Number of cell in x direction");
  args.addOption("grid-size-y", 'y', "Number of cell in y direction");
  args.addOption("output-basepath", 'o', "Output base file name");

  tools::Args::Result ret = args.parse(argc, argv);

  switch (ret)
  {
  case tools::Args::Error:
	  return 1;
  case tools::Args::Help:
	  return 0;
  }

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  l_baseName = args.getArgument<std::string>("output-basepath");

  // create a simple artificial scenario
  SWE_BathymetryDamBreakScenario l_scenario;

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints = 20;

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;
  
  //! local number of grid cells in x- and y-direction
  int l_localnX, l_localnY;
  
  // Compute local number of grid cells and origin
  l_localnX = l_nX / l_numProc;
  l_localnY = l_nY;   // We only have one row at the moment

  // create a single wave propagation block
  SWE_WavePropagationBlock l_wavePropgationBlock(l_localnX,l_localnY,l_dX,l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);
  
  // Set the local origin
  l_originX += (l_dX * l_localnX) * l_myRank;

  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario);


  //! time when the simulation ends.
  float l_endSimulation = l_scenario.endSimulation();

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,l_myRank,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};
  
  // consturct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_localnX, l_localnY,
		  l_dX, l_dY,
		  l_myRank * l_localnX, 0 * l_localnY);
  
  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          (float) 0.);


  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.0;
  progressBar.update(l_t);

  unsigned int l_iterations = 0;
  
  // Compute the neighbor ranks
  int l_leftNeighbor = l_myRank - 1;
  if (l_leftNeighbor < 0)
    l_leftNeighbor = MPI_PROC_NULL;
  int l_rightNeighbor = l_myRank + 1;
  if (l_rightNeighbor >= l_numProc)
    l_rightNeighbor = MPI_PROC_NULL;
  
  // Get the proxy objects for communication
  SWE_Block1D* l_leftInflow  = l_wavePropgationBlock.grabGhostLayer(BND_LEFT);
  SWE_Block1D* l_leftOutflow = l_wavePropgationBlock.registerCopyLayer(BND_LEFT);
  if (l_myRank == 0)
    // grabGhostLayer is automatically changing the boundary type
    // change this back to OUTFLOW if we are the left most neighbor
    l_wavePropgationBlock.setBoundaryType(BND_LEFT, OUTFLOW);

  SWE_Block1D* l_rightInflow  = l_wavePropgationBlock.grabGhostLayer(BND_RIGHT);
  SWE_Block1D* l_rightOutflow = l_wavePropgationBlock.registerCopyLayer(BND_RIGHT);
  if (l_myRank == l_numProc-1)
    l_wavePropgationBlock.setBoundaryType(BND_RIGHT, OUTFLOW);

  // loop over checkpoints
  for(int c=1; c<=l_numberOfCheckPoints; c++) {

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
      // Exchange ghost layers
      exchangeGhostLayers(l_leftNeighbor,  l_leftInflow,  l_leftOutflow,
		          l_rightNeighbor, l_rightInflow, l_rightOutflow,
			  l_localnY+2);
      
      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();
      
      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();
      
      // Find the global maxTimeStepWidth
      MPI_Allreduce(MPI_IN_PLACE, &l_maxTimeStepWidth, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

      // update the cpu time in the logger
      tools::Logger::logger.updateTime("Cpu");

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      l_iterations++;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);
    }

    // print current simulation time of the output
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);

    // write output
    l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                            l_wavePropgationBlock.getDischarge_hu(),
                            l_wavePropgationBlock.getDischarge_hv(),
                            l_t);
  }

  /**
   * Finalize.
   */
  // write the statistics message
  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printTime("Cpu", "CPU time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);
  
  // Finalize MPI
  MPI_Finalize();

  return 0;
}

/**
 * Exchanges the left and right ghost layers with MPI's SendReceive.
 *
 * @param i_leftNeighborRank MPI rank of the  left neighbor.
 * @param o_leftInflow ghost layer, where the left neighbor writes into.
 * @param i_leftOutflow layer where the left neighbor reads from.
 * @param i_rightNeighborRank MPI rank of the right neighbor.
 * @param o_rightInflow ghost layer, where the right neighbor writes into.
 * @param i_rightOutflow layer, where the right neighbor reads form.
 * @param i_boundarySize Number of cells on the boundary.
 */
void exchangeGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                                   const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                                   const int i_boundarySize) {

  MPI_Status l_status;

  // send to left, receive from the right:
  MPI_Sendrecv( i_leftOutflow->h.elemVector(), i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  1,
                o_rightInflow->h.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 1,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hu.elemVector(), i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  2,
                o_rightInflow->hu.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 2,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hv.elemVector(), i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  3,
                o_rightInflow->hv.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 3,
                MPI_COMM_WORLD, &l_status );

  // send to right, receive from the left:
  MPI_Sendrecv( i_rightOutflow->h.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 4,
                o_leftInflow->h.elemVector(),   i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  4,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_rightOutflow->hu.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 5,
                o_leftInflow->hu.elemVector(),   i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  5,
                MPI_COMM_WORLD, &l_status);

  MPI_Sendrecv( i_rightOutflow->hv.elemVector(), i_boundarySize, MPI_FLOAT, i_rightNeighborRank, 6,
                o_leftInflow->hv.elemVector(),   i_boundarySize, MPI_FLOAT, i_leftNeighborRank,  6,
                MPI_COMM_WORLD, &l_status );

}