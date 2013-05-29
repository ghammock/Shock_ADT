/******************************************************************************
||  shock.cpp               (Definition file)                                ||
||===========================================================================||
||                                                                           ||
||    Author: Gary Hammock, PE                                               ||
||    Creation Date: 2013-05-20                                              ||
||    Last Edit Date: 2013-05-28                                             ||
||                                                                           ||
||===========================================================================||
||  FILE DESCRIPTION                                                         ||
||===========================================================================||
||    This file contains the definition details for calculating and storing  ||
||    the thermodynamic properties both upstream and downstream of a normal  ||
||    shock in air.                                                          ||
||                                                                           ||
||===========================================================================||
||  REQUIREMENTS                                                             ||
||===========================================================================||
||    The Air ADT from Hammock.                                              ||
||    shock.cpp                                                              ||
||    shock.h                                                                ||
||                                                                           ||
||===========================================================================||
||  REFERENCES                                                               ||
||===========================================================================||
||    Anderson, John D., "Modern Compressible Flow with Historical           ||
||        Perspective".  Third Edition.  McGraw-Hill Companies, Inc.  2003.  ||
||                                                                           ||
||    Emanuel, George.  "Analytical Fluid Dynamics".  2nd Edition.           ||
||        CRC Press, Boca Raton, FL 2001.  pp 751-753.                       ||
||                                                                           ||
||===========================================================================||
||  LICENSE    (MIT/X11 License)                                             ||
||===========================================================================||
||    Copyright (C) 2013 Gary Hammock                                        ||
||                                                                           ||
||    Permission is hereby granted, free of charge, to any person obtaining  ||
||    a copy of this software and associated documentation files (the        ||
||    "Software"), to deal in the Software without restriction, including    ||
||    without limitation the rights to use, copy, modify, merge, publish,    ||
||    distribute, sublicense, and/or sell copies of the Software, and to     ||
||    permit persons to whom the Software is furnished to do so, subject to  ||
||    the following conditions:                                              ||
||                                                                           ||
||    The above copyright notice and this permission notice shall be         ||
||    included in all copies or substantial portions of the Software.        ||
||                                                                           ||
||    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        ||
||    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     ||
||    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. ||
||    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   ||
||    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   ||
||    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ||
||    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 ||
||                                                                           ||
******************************************************************************/

/**
 *  @file shock.h
 *  @author Gary Hammock, PE
 *  @date 2013-05-28
*/

#ifndef _GH_SHOCK_DEF_H
#define _GH_SHOCK_DEF_H

// We need the math library for sqrt() and pow().
#ifndef _CMATH_
#include <cmath>
#endif

// This include pre-processor gives us access to
// the air properties library.
#include "../air/air.h"

/**
 *  @class Shock An ADT to calculate and store the downstream
 *         parameters of a standing shock based on input upstream parameters.
 *
*/
class Shock
{
  public:
    /******************************************************
    **           Constructors / Destructors              **
    ******************************************************/

    /** Default constructor.  */
    Shock ();

    /** Copy constructor.
     *
     *  @pre none.
     *  @post A new Shock object is instantiated with the copied values.
     *  @param copyFrom The shock object whose values are to be copied.
    */
    Shock (const Shock &copyFrom);

    /** Initialization constructor for a normal shock.
     *
     *  @pre none.
     *  @post A new Shock object is instantiated based on the given
     *        upstream inputs.
     *  @param upstreamStaticPressure The static pressure upstream
     *         of the shock in MPa.
     *  @param upstreamStaticEnthlapy The static enthalpy upstream
     *         of the shock in kJ/kg.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
    */
    Shock (double upstreamStaticPressure, double upstreamStaticEnthalpy,
           double upstreamMach);

    /** Initialization constructor for a normal shock.
     *
     *  @pre none.
     *  @post A new normal Shock object is instantiated based on the
     *        given upstream inputs.
     *  @param upstreamAirProperties An instantiated Air ADT that
     *         contains the thermodynamic properties upstream of the shock.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
    */
    Shock (const Air &upstreamAirProperties, double upstreamMach);

    /** Initialization constructor for an oblique shock.
     *
     *  @pre none.
     *  @post A new Shock object is instantiated based on the given
     *        upstream inputs.
     *  @param upstreamStaticPressure The static pressure upstream
     *         of the shock in MPa.
     *  @param upstreamStaticEnthlapy The static enthalpy upstream
     *         of the shock in kJ/kg.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param shockAngle The shock angle in radians.
    */
    Shock (double upstreamStaticPressure, double upstreamStaticEnthalpy,
           double upstreamMach, double shockAngle);

    /** Initialization constructor for an oblique shock.
     *
     *  @pre none.
     *  @post A new oblique Shock object is instantiated based on the
     *        given upstream inputs.
     *  @param upstreamAirProperties An instantiated Air ADT that
     *         contains the thermodynamic properties upstream of the shock.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param shockAngle The shock angle in radians.
    */
    Shock (const Air &upstreamAirProperties, double upstreamMach,
           double shockAngle);

    /** Default destructor.  */
    ~Shock ();

    /******************************************************
    **               Accessors / Mutators                **
    ******************************************************/

    ////////////////////
    //    Getters
    ////////////////////

    /** Retrieve the thermodynamic conditions of the air upstream of the shock.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return An Air object with values copied from _upstream.
    */
    Air getUpstreamConditions (void) const;

    /** Retrieve the thermodynamic conditions downstream of the shock.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return An Air object with values copied from _downstream.
    */
    Air getDownstreamConditions (void) const;

    /** Retrieve the Mach number of the flow upstream of the shock.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value stored in _upstreamMach.
    */
    double getUpstreamMach (void) const;

    /** Retrieve the Mach number of the flow downstream of the shock.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value stored in _downstreamMach.
    */
    double getDownstreamMach (void) const;

    /** Retrieve the upstream total pressure [units: MPa].
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _upstreamTotalPressure.
    */
    double getUpstreamTotalPressure (void) const;

    /** Retrieve the downstream total pressure [units: MPa].
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _downstreamTotalPressure.
    */
    double getDownstreamTotalPressure (void) const;

    /** Retrieve the total pressure ratio of the shock
     *  (downstream total pressure to upstream total pressure).
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _totalPressureRatio.
    */
    double getTotalPressureRatio (void) const;

    /** Retrieve the flow deflection angle (the angular deflection
     *  due to the flow passing through a compression wave).
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _downstreamFlowAngle [units: radians].
    */
    double getDownstreamFlowAngle (void) const;

    /** Retrieve the shock angle.
     *
     *  @pre The object is instantiated.
     *  @post none.
     *  @return The value of _shockAngle [units: radians].
    */
    double getShockAngle (void) const;

    ////////////////////
    //    Setters
    ////////////////////

    /** Calculate the downstream properties of the normal shock.
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamStaticPressure The static pressure upstream
     *         of the shock in MPa.
     *  @param upstreamStaticEnthlapy The static enthalpy upstream
     *         of the shock in kJ/kg.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @return none.
    */
    void calculateDownstream (double upstreamStaticPressure,
                              double upstreamStaticEnthalpy,
                              double upstreamMach);

    /** Calculate the downstream properties of the normal shock.
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamAirProperties An instantiated and filled Air ADT that
     *         contains the thermodynamic properties upstream of the shock
     *         (specifically temperature and pressure).
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @return none.
    */
    void calculateDownstream (const Air &upstreamAirProperties,
                              double upstreamMach);

    /** Calculate the downstream properties of the oblique shock.
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamStaticPressure The static pressure upstream
     *         of the shock in MPa.
     *  @param upstreamStaticEnthlapy The static enthalpy upstream
     *         of the shock in kJ/kg.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param shockAngle The shock angle in radians.
     *  @return none.
    */
    void calculateDownstream (double upstreamStaticPressure,
                              double upstreamStaticEnthalpy,
                              double upstreamMach,
                              double shockAngle);

    /** Calculate the downstream properties of the oblique shock.
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamAirProperties An instantiated and filled Air ADT that
     *         contains the thermodynamic properties upstream of the shock
     *         (specifically temperature and pressure).
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param shockAngle The shock angle in radians.
     *  @return none.
    */
    void calculateDownstream (const Air &upstreamAirProperties,
                              double upstreamMach, double shockAngle);

    /** Calculate the downstream properties of a weak oblique shock given
     *  a downstream flow deflection angle (the weak shock is favored and
     *  usually occurs in nature).
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamStaticPressure The static pressure upstream
     *         of the shock in MPa.
     *  @param upstreamStaticEnthlapy The static enthalpy upstream
     *         of the shock in kJ/kg.
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param flowDeflectionAngle The flow deflection angle in radians.
     *  @return none.
    */
    void calculateDownstream_weakShock (double upstreamStaticPressure,
                                        double upstreamStaticEnthalpy,
                                        double upstreamMach,
                                        double flowDeflectionAngle);

    /** Calculate the downstream properties of a weak oblique shock given
     *  a downstream flow deflection angle (the weak shock is favored and
     *  usually occurs in nature).
     *
     *  @pre The object is instantiated.
     *  @post The downstream properties are calculated.
     *  @param upstreamAirProperties An instantiated and filled Air ADT that
     *         contains the thermodynamic properties upstream of the shock
     *         (specifically temperature and pressure).
     *  @param upstreamMach The Mach number of the flow upstream
     *         of the shock.
     *  @param flowDeflectionAngle The flow deflection angle in radians.
     *  @return none.
    */
    void calculateDownstream_weakShock (const Air &upstreamAirProperties,
                                        double upstreamMach,
                                        double flowDeflectionAngle);

    /** Clear/Reinitialize the members of the shock to zero.
     *
     *  @pre The object is instantiated.
     *  @post All the members have their values cleared/reset.
     *  @return none.
    */
    void clear (void);

  private:
    /******************************************************
    **                      Members                      **
    ******************************************************/

    Air _upstream,    // Upstream air properties.
        _downstream;  // Downstream air properties.

    double _upstreamMach,     // The upstream Mach number (pre-shock).
           _downstreamMach,   // The downstream Mach number (post-shock).
           _upstreamTotalPressure,   // Upstream total P [units: MPa].
           _downstreamTotalPressure, // Downstream total P [units: MPa].
           _totalPressureRatio,  // Downstream-to-Upstream total P ratio.
           _downstreamFlowAngle, // Downstream flow deflection angle in rad.
           _shockAngle;       // The shock angle in radians.
                              //   pi/2 = normal shock

    static const double pi;  // 3.1415926...

    /******************************************************
    **                   Helper Methods                  **
    ******************************************************/

    /** Calculate an oblique shock angle given the gas ratio of specific
     *  heats, upstream mach number, and flow deflection angle.
     *
     *  @pre none.
     *  @post none.
     *  @param gamma The ratio of specific heats of the gas.
     *  @param upstreamMach The mach number of the flow upstream of the shock.
     *  @param flowDeflectionAngle The angle (in radians) of the deflected
     *         flow downstream of the shock.
     *  @param weakShock true to calculate the weak shock (most common)
     *         solution, false to calculate the strong shock (less common).
     *  @return A double precision value of the calculated shock angle
     *          [units: radians].
    */
    double _calculateShockAngle (double gamma, double upstreamMach,
                                 double flowDeflectionAngle,
                                 bool weakShock) const;

};  // End class Shock

#endif