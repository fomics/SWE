/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * Benchmarks scenarios for the SWE code.
 */

#ifndef __BENCHMARKS_H
#define __BENCHMARKS_H

#include <cmath>

#include "SWE_Scenario.hh"

/**
 * Dam-break with wet dry interface, aligned to x-dimension
 */
class SWE_DamBreakWetDryXaligned : public SWE_Scenario {

  public:
    float getBathymetry(float x, float y) {
       return 0;
    };

    float getWaterHeight(float x, float y) { 
      if( x < 0 ) {
        return 10;
      }
      else {
        return 0;
      }
    };

	  virtual float endSimulation() { return (float) 4; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)-10;
       else if ( i_edge == BND_RIGHT)
         return (float) 10;
       else if ( i_edge == BND_BOTTOM )
         return (float)-1;
       else
         return (float) 1;
    };
};

/**
 * Dam-break with wet dry interface, not aligned to a spatial dimension.
 */
class SWE_DamBreakWetDryRotated: public SWE_Scenario {

  public:
    float getBathymetry(float x, float y) {
       return 0;
    };

    float getWaterHeight(float x, float y) { 
      if( y < 1.5 - 0.5*x ) {
        return 10;
      }
      else {
        return 0;
      }
    };

	  virtual float endSimulation() { return (float) 4; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)-10;
       else if ( i_edge == BND_RIGHT)
         return (float) 10;
       else if ( i_edge == BND_BOTTOM )
         return (float)-10;
       else
         return (float) 10;
    };
};
#endif
