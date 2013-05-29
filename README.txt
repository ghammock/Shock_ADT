Shock ADT
Author: Gary Hammock, PE

================================================================================
                                DESCRIPTION
================================================================================

This C++ abstract data type is used to calculate the thermodynamic and
transport properties of equilibrium air upstream and downstream of a shock wave.
The shock can either be a normal shock or an oblique shock.

The current implementation requires the Air ADT
(https://github.com/ghammock/Air_ADT) for the calculation of properties.

================================================================================
                               IMPLEMENTATION
================================================================================

To create a new Shock object:

                   Shock sw;

List of accessor methods used by the ADT:

    
================================================================================
                              DESIRED UPDATES
================================================================================
Rewrite the air-specific parts of the ADT to be fluid independent.  This may
require using a generic fluid ADT as a base class from which the Air ADT could
be a derivative.

================================================================================
                                 REFERENCES
================================================================================
1.) Anderson, John D., "Modern Compressible Flow with Historical Perspective".
      Third Edition.  McGraw-Hill Companies, Inc.  2003.

2.) Emanuel, George.  "Analytical Fluid Dynamics".  2nd Edition. CRC Press,
      Boca Raton, FL 2001.  pp 751-753.

3.) Liepmann, H. W., A. Roshko.  "Elements of Gasdynamics".  Dover Publications.
      2001.  (original copyright: New York.  John Wiley & Sons.  1957.)
      ISBN 978-0-486-41963-3.