/*
 * This file is part of SWE.
 *
 *         Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#ifndef CUDA
#include "blocks/SWE_WavePropagationBlock.hh"
#else
#include "blocks/cuda/SWE_WavePropagationBlockCuda.hh"
#endif

#ifdef WRITENETCDF
#include "writer/NetCdfWriter.hh"
#else
#include "writer/VtkWriter.hh"
#endif

#ifdef ASAGI
#include "scenarios/SWE_AsagiScenario.hh"
#else
#include "scenarios/SWE_simple_scenarios.hh"
#endif

#ifdef READXML
#include "tools/CXMLConfig.hpp"
#endif

#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"

/*
 * A simple progress bar using stdout
 */

#include <cassert>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <limits>

#ifndef WIN32
 #include <unistd.h>
 #include <sys/ioctl.h>
#endif

namespace tools
{

class ProgressBar
{
private:
	/* Local rank (we only do work on rank 0) */
	int m_rank;

	/* Total amount of work */
	float m_totalWork;

	/* Progress bar initialization time */
	time_t m_startTime;

	/* Terminal size */
	unsigned int m_terminalSize;

	/* Rotating bar char */
	unsigned char m_rotatingBar;

public:
	ProgressBar(float totalWork = 1., int rank = 0)
		: m_rank(rank),
		  m_totalWork(totalWork),
		  m_startTime(time(0)),
		  m_rotatingBar(0)
	{
		if (rank != 0)
			return;

#ifdef TIOCGSIZE
		struct ttysize ts;
		ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
		m_terminalSize = ts.ts_cols;
#elif defined(TIOCGWINSZ)
		struct winsize ts;
		ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
		m_terminalSize = ts.ws_col;
#else
		m_terminalSize = 0;
#endif
		if (m_terminalSize > 300)
			// Probably an error due to MPI
			m_terminalSize = MIN_TERM_SIZE;
	}

	/*
	 * Update the progress bar with done as the amount of work already done
	 */
	void update(float done)
	{
		if (m_rank != 0 || m_terminalSize < MIN_TERM_SIZE)
			return;

		unsigned int printed = 2;
		std::cout << '\r';
		printed += printTimeLeft(done);
		std::cout << ' ';
		printed += printPercentage(done);
		std::cout << ' ';
		printProgressBar(done, m_terminalSize-printed-2);
		std::cout << ' ';
		printRotatingBar();
		std::cout << std::flush;
	}

	void clear()
	{
		if (m_rank != 0 || m_terminalSize < MIN_TERM_SIZE)
			return;

		std::cout << '\r';
		for (unsigned int i = 0; i < m_terminalSize; i++)
			std::cout << ' ';
		std::cout << '\r';
	}

private:
	/*
	 * Return the number of characters printed
	 */
	unsigned int printTimeLeft(float done)
	{
		float timeLeft;
		if (done <= 0)
			timeLeft = std::numeric_limits<float>::max();
		else
			timeLeft = (time(0) - m_startTime) * (m_totalWork - done) / done;

		std::cout << "Time left: ";

		if (timeLeft < 1) {
			for (int i = 3; i < TIME_SIZE; i++)
				std::cout << ' ';
			std::cout << "< 1";
		} else {
			int digits = ceil(log(timeLeft)/log(10));
			if (digits > TIME_SIZE) {
				// Maximum number we can show
				for (int i = 0; i < TIME_SIZE; i++)
					std::cout << '9';
			} else {
				streamsize oldPrec = std::cout.precision();
				std::ios::fmtflags oldFlags = std::cout.flags();
				streamsize oldWidth = std::cout.width();

				std::cout.precision(std::max(0, TIME_SIZE-digits-2));
				std::cout.setf(std::ios::fixed);
				std::cout.width(TIME_SIZE);

				std::cout << timeLeft;

				std::cout.precision(oldPrec);
				std::cout.flags(oldFlags);
				std::cout.width(oldWidth);
			}
		}

		std::cout << " sec";

		return 11+TIME_SIZE+4;
	}

	/*
	 * Return the number of characters printed
	 */
	unsigned int printPercentage(float done)
	{
		int per = floor(done/m_totalWork*100);

		std::cout << '(';

		streamsize oldWidth = std::cout.width();

		std::cout.width(3);
		std::cout << per;

		std::cout.width(oldWidth);

		std::cout << "% done)";

		return 1+3+7;
	}

	void printProgressBar(float done, unsigned int size)
	{
		if (size < 3)
			return;

		size -= 2; // leave space for []
		unsigned int per = floor(done/m_totalWork * size);

		std::cout << '[';

		for (unsigned int i = 0; i < per; i++)
			std::cout << '=';

		if (per < size) {
			std::cout << '>';
			per++;
		}

		for (unsigned int i = per; i < size; i++)
			std::cout << ' ';

		std::cout << ']';
	}

	void printRotatingBar()
	{
		static const char* CHARS = "|/-\\";

		std::cout << CHARS[m_rotatingBar];

		m_rotatingBar = (m_rotatingBar + 1) % 4;
	}

	static const unsigned int MIN_TERM_SIZE = 80;
	static const int TIME_SIZE = 8;
};

}

/*
 * Main program (obviously)
 */
int main( int argc, char** argv ) {
  #ifndef READXML
  // Parse command line parameters
  tools::Args args;
  args.addOption("grid-size-x", 'x', "Number of cell in x direction");
  args.addOption("grid-size-y", 'y', "Number of cell in y direction");
  args.addOption("output-basepath", 'o', "Output base file name");
  #endif

  tools::Args::Result ret = args.parse(argc, argv);

  switch (ret)
  {
  case tools::Args::Error:
	  return 1;
  case tools::Args::Help:
	  return 0;
  }

  int l_nX, l_nY;
  std::string l_baseName;

  #ifndef READXML
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  l_baseName = args.getArgument<std::string>("output-basepath");
  #endif

  #ifdef READXML
  assert(false); //TODO: not implemented.
  if(argc != 2) {
    s_sweLogger.printString("Aborting. Please provide a proper input file.");
    s_sweLogger.printString("Example: ./SWE_gnu_debug_none_augrie config.xml");
    return 1;
  }
  s_sweLogger.printString("Reading xml-file.");

  std::string l_xmlFile = std::string(argv[1]);
  s_sweLogger.printString(l_xmlFile);

  CXMLConfig l_xmlConfig;
  l_xmlConfig.loadConfig(l_xmlFile.c_str());
  #endif

  #ifdef ASAGI
  /* Information about the example bathymetry grid (tohoku_gebco_ucsb3_500m_hawaii_bath.nc):
   *
   * Pixel node registration used [Cartesian grid]
   * Grid file format: nf = GMT netCDF format (float)  (COARDS-compliant)
   * x_min: -500000 x_max: 6500000 x_inc: 500 name: x nx: 14000
   * y_min: -2500000 y_max: 1500000 y_inc: 500 name: y ny: 8000
   * z_min: -6.48760175705 z_max: 16.1780223846 name: z
   * scale_factor: 1 add_offset: 0
   * mean: 0.00217145586762 stdev: 0.245563641735 rms: 0.245573241263
   */

  float simulationArea[4];
  simulationArea[0] = -450000;
  simulationArea[1] = 6450000;
  simulationArea[2] = -2450000;
  simulationArea[3] = 1450000;

  SWE_AsagiScenario l_scenario( ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
                                ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
                                (float) 28800., simulationArea);
  #else
  // create a simple artificial scenario
  SWE_BathymetryDamBreakScenario l_scenario;
  #endif

  int ncpts = 20;

  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  #ifndef CUDA
  SWE_WavePropagationBlock wpb(l_nX,l_nY,l_dX,l_dY);
  #else
  SWE_WavePropagationBlockCuda wpb(l_nX,l_nY,l_dX,l_dY);
  #endif

  float l_originX, l_originY;

  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

  wpb.initScenario(l_originX, l_originY, l_scenario);

  float l_endSimulation = l_scenario.endSimulation();
  float* ckpts = new float[ncpts+1];

  for(int cp = 0; cp <= ncpts; cp++) {
  ckpts[cp] = cp*(l_endSimulation/ncpts);
  }

  tools::ProgressBar progressBar(l_endSimulation);

  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string fn = generateBaseFileName(l_baseName,0,0);
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};
#ifdef WRITENETCDF
  //construct a NetCdfWriter
  io::NetCdfWriter l_writer( fn,
		  wpb.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY,
		  l_originX, l_originY);
#else
  // consturct a VtkWriter
  io::VtkWriter l_writer( fn,
		  wpb.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY );
#endif
  // Write zero time step
  l_writer.writeTimeStep( wpb.getWaterHeight(),
                          wpb.getDischarge_hu(),
                          wpb.getDischarge_hv(),
                          (float) 0.);


  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  float l_t = 0.0;
  progressBar.update(l_t);

  unsigned int iter = 0;

  for(int c=1; c<=ncpts; c++) {

    while( l_t < ckpts[c] ) {
      wpb.setGhostLayer();
      
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // approximate the maximum time step
      // TODO: This calculation should be replaced by the usage of the wave speeds occuring during the flux computation
      // Remark: The code is executed on the CPU, therefore a "valid result" depends on the CPU-GPU-synchronization.
//      wpb.computeMaxTimestep();

      // compute numerical flux on each edge
      wpb.computeNumericalFluxes();

      float l_maxTimeStepWidth = wpb.getMaxTimestep();

      wpb.updateUnknowns(l_maxTimeStepWidth);

      tools::Logger::logger.updateTime("Cpu");

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      iter++;

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
    l_writer.writeTimeStep( wpb.getWaterHeight(),
                            wpb.getDischarge_hu(),
                            wpb.getDischarge_hv(),
                            l_t);
  }

  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();

  tools::Logger::logger.printTime("Cpu", "CPU time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(iter);

  return 0;
}
