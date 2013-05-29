/******************************************************************************
||  shock.cpp               (Implementation file)                            ||
||===========================================================================||
||                                                                           ||
||    Author: Gary Hammock, PE                                               ||
||    Creation Date: 2013-05-20                                              ||
||    Last Edit Date: 2013-05-28                                             ||
||                                                                           ||
||===========================================================================||
||  FILE DESCRIPTION                                                         ||
||===========================================================================||
||    This file contains the implementation details for calculating and      ||
||    storing the thermodynamic properties both upstream and downstream of   ||
||    a normal shock in air.                                                 ||
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
 *  @file shock.cpp
 *  @author Gary Hammock, PE
 *  @date 2013-05-28
*/

#include "shock.h"

const double Shock::pi = 3.141592653589793246;

/******************************************************
**           Constructors / Destructors              **
******************************************************/

/** Default constructor.  */
Shock::Shock ()
    : _upstreamMach(0.0), _downstreamMach(0.0),
      _upstreamTotalPressure(0.0), _downstreamTotalPressure(0.0),
      _totalPressureRatio(0.0), _downstreamFlowAngle(0.0),
      _shockAngle(pi / 2.0)
{}

/** Copy constructor.
 *
 *  @pre none.
 *  @post A new Shock object is instantiated with the copied values.
 *  @param copyFrom The shock object whose values are to be copied.
*/
Shock::Shock (const Shock &copyFrom)
    : _upstream(copyFrom._upstream), _downstream(copyFrom._downstream),
      _upstreamMach(copyFrom._upstreamMach),
      _downstreamMach(copyFrom._downstreamMach),
      _upstreamTotalPressure(copyFrom._upstreamTotalPressure),
      _downstreamTotalPressure(copyFrom._downstreamTotalPressure),
      _totalPressureRatio(copyFrom._totalPressureRatio),
      _downstreamFlowAngle(copyFrom._downstreamFlowAngle),
      _shockAngle(copyFrom._shockAngle)
{}

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
Shock::Shock (double upstreamStaticPressure, double upstreamStaticEnthalpy,
              double upstreamMach)
{
    calculateDownstream(upstreamStaticPressure,
                        upstreamStaticEnthalpy,
                        upstreamMach);
}

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
Shock::Shock (const Air &upstreamAirProperties, double upstreamMach)
{
    calculateDownstream(upstreamAirProperties, upstreamMach);
}

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
Shock::Shock (double upstreamStaticPressure, double upstreamStaticEnthalpy,
              double upstreamMach, double shockAngle)
{
}

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
Shock::Shock (const Air &upstreamAirProperties, double upstreamMach,
              double shockAngle)
{
}

/** Default destructor.  */
Shock::~Shock () {}

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
Air Shock::getUpstreamConditions (void) const
{  return _upstream;  }

/** Retrieve the thermodynamic conditions downstream of the shock.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return An Air object with values copied from _downstream.
*/
Air Shock::getDownstreamConditions (void) const
{  return _downstream;  }

/** Retrieve the Mach number of the flow upstream of the shock.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value stored in _upstreamMach.
*/
double Shock::getUpstreamMach (void) const
{  return _upstreamMach;  }

/** Retrieve the Mach number of the flow downstream of the shock.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value stored in _downstreamMach.
*/
double Shock::getDownstreamMach (void) const
{  return _downstreamMach;  }

/** Retrieve the upstream total pressure [units: MPa].
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _upstreamTotalPressure.
*/
double Shock::getUpstreamTotalPressure (void) const
{  return _upstreamTotalPressure;  }

/** Retrieve the downstream total pressure [units: MPa].
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _downstreamTotalPressure.
*/
double Shock::getDownstreamTotalPressure (void) const
{  return _downstreamTotalPressure;  }

/** Retrieve the total pressure ratio of the shock
 *  (downstream total pressure to upstream total pressure).
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _totalPressureRatio.
*/
double Shock::getTotalPressureRatio (void) const
{  return _totalPressureRatio;  }

/** Retrieve the flow deflection angle (the angular deflection
 *  due to the flow passing through a compression wave).
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _downstreamFlowAngle [units: radians].
*/
double Shock::getDownstreamFlowAngle (void) const
{  return _downstreamFlowAngle;  }

/** Retrieve the shock angle.
 *
 *  @pre The object is instantiated.
 *  @post none.
 *  @return The value of _shockAngle [units: radians].
*/
double Shock::getShockAngle (void) const
{  return _shockAngle;  }

////////////////////
//    Setters
////////////////////

/** Calculate the downstream properties of the normal shock.
 *
 *  @pre The object is instantiated.
 *  @post The downstream properties are calculated.
 *  @param upstreamTotalPressure The total pressure upstream
 *         of the shock in MPa.
 *  @param upstreamTotalEnthlapy The total enthalpy upstream
 *         of the shock in kJ/kg.
 *  @param upstreamMach The Mach number of the flow upstream
 *         of the shock.
 *  @return none.
*/
void Shock::calculateDownstream (double upstreamStaticPressure,
                                 double upstreamStaticEnthalpy,
                                 double upstreamMach)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    // Calculate the upstream properties
    Air upstream;
    upstream.calculateProps_PH(upstreamStaticPressure, upstreamStaticEnthalpy);

    // Now use the upstream properties to calculate the downstream properties
    // using the Air ADT and the polymorphic calculateDownstream() method.
    calculateDownstream(upstream, upstreamMach);

    return;
}

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
void Shock::calculateDownstream (const Air &upstreamAirProperties,
                                 double upstreamMach)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    double p2_p1,   // Static pressure ratio.
           pt2_pt1, // Total pressure ratio.
           p_pt,    // Static-to-total pressure ratio.
           T2_T1,   // Static temperature ratio.
           T_Tt;    // Static-to-total temperature ratio.

    // Store the input upstream Mach number.
    _upstreamMach = upstreamMach;

    // Since we're dealing with a normal shock, ensure that the shock
    // angle is 90 degrees (pi / 2 radians).  The downstream flow is
    // not deflected by the shock wave so the flow angle is 0 degrees.
    _shockAngle = pi / 2.0;
    _downstreamFlowAngle = 0.0;

    // Store the upstream state.
    _upstream = upstreamAirProperties;

    // Store the upstream gamma value
    double gamma = _upstream.getGamma();

    // Store a few temporary values to aid readability.
    double gm1 = gamma - 1.0,  // (used often)
           gp1 = gamma + 1.0,  // (used often)
           gamMach[2];  // 0 = upstream, 1 = downstream

    // We need to calculate the downstream Mach number.
    gamMach[0] = 0.5 * gm1 * pow(_upstreamMach, 2.0);

    _downstreamMach = sqrt((1.0 + gamMach[0])
                       / ((gamma * pow(_upstreamMach, 2.0)) - (0.5 * gm1)));

    gamMach[1] = 0.5 * gm1 * pow(_downstreamMach, 2.0);

    // Calculate the upstream static-to-total temperature and pressure ratios.
    T_Tt = 1.0 / (1.0 + (0.5 * gm1 * pow(_upstreamMach, 2.0)));
    p_pt = pow(T_Tt, (gamma / gm1));

    // Calculate the upstream total pressure [Units: MPa].
    double upstreamStaticPressure = _upstream.getPressure();
    _upstreamTotalPressure = upstreamStaticPressure / p_pt;

    // Calculate the ratio of static pressures across a shock.
    p2_p1 = 1 + (2.0 * gamma * (pow(_upstreamMach, 2.0) - 1.0) / gp1);

    // Calculate the ratio of total pressure across a shock to determine
    // the Pitot pressure.
    pt2_pt1 = pow((1.0 + gamMach[1]) / (1.0 + gamMach[0]), (gamma / gm1));
    pt2_pt1 *= p2_p1;

    _totalPressureRatio = pt2_pt1;

    // Calculate the ratio of static temperatures across a shock.
    T2_T1 = (1.0 + gamMach[0]) / (1.0 + gamMach[1]);

    double downstreamStaticPressure,    // The calculated downstream static P.
           downstreamStaticTemperature; // The calculated downstream static T.

    // Calculate the total pressure using the calculated
    // total pressure ratio.  [Units: MPa]
    _downstreamTotalPressure  = pt2_pt1 * _upstreamTotalPressure;
    downstreamStaticPressure = p2_p1 * upstreamStaticPressure;

    // Calculate the static temperature using the
    // calculated temperature ratios.  [Units: K]
    downstreamStaticTemperature = T2_T1 * _upstream.getTemperature();

    // Calculate the thermodynamic properties downstream of the shock.
    _downstream.calculateProperties(downstreamStaticPressure,
                                    downstreamStaticTemperature);

    return;
}

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
void Shock::calculateDownstream (double upstreamStaticPressure,
                                 double upstreamStaticEnthalpy,
                                 double upstreamMach,
                                 double shockAngle)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    // Calculate the upstream properties
    Air upstream;
    upstream.calculateProps_PH(upstreamStaticPressure, upstreamStaticEnthalpy);

    // Now use the upstream properties to calculate the downstream properties
    // using the Air ADT and the polymorphic calculateDownstream() method.
    calculateDownstream(upstream, upstreamMach, shockAngle);

    return;
}

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
void Shock::calculateDownstream (const Air &upstreamAirProperties,
                                 double upstreamMach, double shockAngle)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    double p2_p1,   // Static pressure ratio.
           pt2_pt1, // Total pressure ratio.
           p_pt,    // Static-to-total pressure ratio.
           T2_T1,   // Static temperature ratio.
           T_Tt;    // Static-to-total temperature ratio.

    ///////////////////////////////////////////////////////
    // Check that the shock angle is in the valid range.
    // Correcting if necessary.
    ///////////////////////////////////////////////////////

    // Corrects for "nearly normal".
    if (fabs(fabs(shockAngle * 180.0 / pi) - 90.0) < 1E-4)
        shockAngle = pi / 2.0;

    // The arguments for atan(x) are -pi/2 <= x <= pi/2.
    // Also, the maximum wave angle is pi/2.
    if (shockAngle > (pi / 2.0))
        shockAngle -= pi;
    else if (shockAngle < (-pi / 2.0))
        shockAngle += pi;

    // The minimum wave angle is the mach angle.  Correct if a
    // smaller angle is given.
    if (shockAngle < asin(1.0 / upstreamMach))
        shockAngle = asin(1.0 / upstreamMach);

    ///////////////////////////////////////////////////////
    //  End checking/correcting
    ///////////////////////////////////////////////////////
    
    // Store the upstream state.
    _upstream = upstreamAirProperties;

    // Store the upstream gamma value
    double gamma = _upstream.getGamma();

    // Store the input upstream Mach number.
    _upstreamMach = upstreamMach;

    // Since we're dealing with an oblique shock, we need to know the
    // shock angle.  This will determine the downstream flow deflection angle.
    _shockAngle = shockAngle;

    // The oblique shock relations are governed by the velocity
    // of the fluid that is *normal* to the flow.  This allows us
    // to use the same isentropic equations as for the normal shock
    // substituting the normal component of the velocity in place of
    // the upstream Mach number.
    double us_machNormal = _upstreamMach * sin(shockAngle);

    double A, B, C;  // Temporary variables used in deflection angle calc.

    A = 2.0 / tan(shockAngle);
    B = (pow(upstreamMach, 2.0) * pow(sin(shockAngle), 2.0)) - 1.0;
    C = (pow(upstreamMach, 2.0) * (gamma + cos(2.0 * shockAngle))) + 2.0;

    _downstreamFlowAngle = atan(A * B / C);  // Units: radians

    // Store a few temporary values to aid readability.
    double gm1 = gamma - 1.0,  // (used often)
           gp1 = gamma + 1.0;  // (used often)

    // We need to calculate the downstream Mach number.
    double ds_machNormal;  // The downstream mach number normal to the shock.

    ds_machNormal  = pow(us_machNormal, 2.0) + (2.0 / gm1);
    ds_machNormal /= ((2.0 * gamma * pow(us_machNormal, 2.0)) / gm1) - 1.0;
    ds_machNormal  = sqrt(ds_machNormal);

    _downstreamMach = ds_machNormal / sin(shockAngle - _downstreamFlowAngle);

    // Used to aid readability for the ratio calculation.
    double gamMach[2];  // 0 = upstream, 1 = downstream
    gamMach[0] = 0.5 * gm1 * pow(_upstreamMach, 2.0);
    gamMach[1] = 0.5 * gm1 * pow(_downstreamMach, 2.0);

    // Calculate the upstream static-to-total temperature and pressure ratios.
    T_Tt = 1.0 / (1.0 + (0.5 * gm1 * pow(_upstreamMach, 2.0)));
    p_pt = pow(T_Tt, (gamma / gm1));

    // Calculate the upstream total pressure [Units: MPa].
    double upstreamStaticPressure = _upstream.getPressure();
    _upstreamTotalPressure = upstreamStaticPressure / p_pt;

    // Calculate the ratio of static pressures across a shock.
    p2_p1 = 1 + (2.0 * gamma * (pow(us_machNormal, 2.0) - 1.0) / gp1);

    // Calculate the ratio of total pressure across a shock to determine
    // the Pitot pressure.
    pt2_pt1 = pow((1.0 + gamMach[1]) / (1.0 + gamMach[0]), (gamma / gm1));
    pt2_pt1 *= p2_p1;

    _totalPressureRatio = pt2_pt1;

    // Calculate the ratio of static temperatures across a shock.
    T2_T1 = (1.0 + gamMach[0]) / (1.0 + gamMach[1]);

    double downstreamStaticPressure,    // The calculated downstream static P.
           downstreamStaticTemperature; // The calculated downstream static T.

    // Calculate the total pressure using the calculated
    // total pressure ratio.  [Units: MPa]
    _downstreamTotalPressure  = pt2_pt1 * _upstreamTotalPressure;
    downstreamStaticPressure = p2_p1 * upstreamStaticPressure;

    // Calculate the static temperature using the
    // calculated temperature ratios.  [Units: K]
    downstreamStaticTemperature = T2_T1 * _upstream.getTemperature();

    // Calculate the thermodynamic properties downstream of the shock.
    _downstream.calculateProperties(downstreamStaticPressure,
                                    downstreamStaticTemperature);

    return;
}

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
void Shock::calculateDownstream_weakShock (double upstreamStaticPressure,
                                           double upstreamStaticEnthalpy,
                                           double upstreamMach,
                                           double flowDeflectionAngle)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    // Calculate the upstream properties
    Air upstream;
    upstream.calculateProps_PH(upstreamStaticPressure, upstreamStaticEnthalpy);

    // Now use the upstream properties to calculate the downstream properties
    // using the Air ADT and the polymorphic calculateDownstream_weakShock()
    // method.
    calculateDownstream_weakShock(upstream, upstreamMach, flowDeflectionAngle);

    return;
}

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
void Shock::calculateDownstream_weakShock (const Air &upstreamAirProperties,
                                           double upstreamMach,
                                           double flowDeflectionAngle)
{
    ///////////////////////////////////////////////////////
    //    NOTE: Always evaluate properties using the     //
    //    static values!                                 //
    ///////////////////////////////////////////////////////

    // We need to calculate the shock angle...
    double shockAngle = _calculateShockAngle(upstreamAirProperties.getGamma(),
                                             upstreamMach,
                                             flowDeflectionAngle,
                                             true);

    // ...then we can calculate the downstream properties as we normally would.
    calculateDownstream(upstreamAirProperties, upstreamMach, shockAngle);

    return;
}

/** Clear/Reinitialize the members of the shock to zero.
 *
 *  @pre The object is instantiated.
 *  @post All the members have their values cleared/reset.
 *  @return none.
*/
void Shock::clear (void)
{
    _upstream.reset();
    _downstream.reset();

    _upstreamMach            = 0.0;
    _downstreamMach          = 0.0;
    _upstreamTotalPressure   = 0.0;
    _downstreamTotalPressure = 0.0;
    _totalPressureRatio      = 0.0;
    _downstreamFlowAngle     = 0.0;
    _shockAngle              = pi / 2.0;

    return;
}

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
double Shock::_calculateShockAngle (double gamma, double upstreamMach,
                                    double flowDeflectionAngle,
                                    bool weakShock) const
{
    // Calculate the shock angle using the equation developed by Emanuel.
    //
    // Reference:
    // ------------------------------------------------------------------------
    // Emanuel, George.  "Analytical Fluid Dynamics".  2nd Edition.
    //    CRC Press, Boca Raton, FL 2001.  pp 751-753.
    //
    // Anderson, John D., "Modern Compressible Flow with Historical
    //    Perspective".  Third Edition.  McGraw-Hill Companies, Inc.  2003.

    double gp1 = gamma + 1.0,
           gm1 = gamma - 1.0,
           lambda,      // Argument for shock angle solution.
           chi,         // Argument for shock angle solution.
           delta,       // Flag for which solution (weak/strong) to solve.
           shockAngle;  // The calculated shock angle [units: radians].

    if (weakShock)
        delta = 1.0;  // Calculates the weak shock solution.
    else
        delta = 0.0;  // Calculates the strong shock solution.

    // Temporary/readability variables.
    double M2 = pow(upstreamMach, 2.0),
           A,
           B,
           C;

    // Calculate lambda using the temporary variables.
    A  = 1.0 + (0.5 * gp1 * M2);
    A  = A * pow(tan(flowDeflectionAngle), 2.0);
    B  = 1.0 + (0.5 * gm1 * M2);
    B  = B * 3.0;
    C  = pow((M2 - 1.0), 2.0) - (A * B);
    lambda  = sqrt(C);

    // Calculate chi using the temporary variables.
    A = 1.0 + (0.5  * gm1 * M2) + (0.25 * gp1 * pow(M2, 2.0));
    A = A * pow(tan(flowDeflectionAngle), 2.0);

    B = 9.0 * (1.0 + (0.5 * gm1 * M2));

    C = pow((M2 - 1.0), 3.0);

    chi = (C - (B * A)) / pow(lambda, 3.0);

    // Now we can calculate the shock angle using lambda and chi.
    A = ((4.0 * pi * delta) + acos(chi)) / 3.0;
    B = M2 - 1.0 + (2.0 * lambda * cos(A));
    C = 3.0 * (1.0 + (0.5 * gm1 * M2)) * tan(flowDeflectionAngle);

    shockAngle = atan(B / C);

    return shockAngle;
}