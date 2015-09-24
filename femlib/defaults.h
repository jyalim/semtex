/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.
 * All parameters are now in a single table, everything is double precision.
 *
 * $Id: defaults.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

static struct {
  char*  name;
  double cval;
} consts[] = {

  /* -- Mathematical constants. */

  "E"           ,   2.71828182845904523536 ,
  "DEG"         ,  57.29577951308232087721 ,
  "PI"          ,   3.14159265358979323844 ,
  "TWOPI"       ,   6.28318530717958647688 ,
  "EULER"       ,   0.57721566490153286061 ,
  "GOLDEN"      ,   1.61803398874989484820 ,
  
  /* -- Default named parameters. */

  "t"           ,   0.0    ,	/* -- Time.                               */
  "D_T"         ,   0.01   ,	/* -- Time step.                          */

  "TOL_REL"     ,   1.0e-8 ,	/* -- Relative tolerance (PCG)            */
  "TOL_ABS"     ,   1.0e-8 ,	/* -- Absolute tolerance.                 */
  "TOL_POS"     ,   1.0e-5 ,    /* -- Positional tolerance.               */

  "z"           ,   0.0    ,	/* -- z-plane location.                   */
  "BETA"        ,   1.0    ,	/* -- TWOPI / Lz (Fourier constant).      */
  "LAMBDA2"     ,   0.0    ,	/* -- Helmholtz constant.                 */

  "KINVIS"      ,   1.0    ,	/* -- Kinematic viscosity.                */
  "REFVIS"      ,   1.0    ,	/* -- Reference kinematic viscosity.      */
  "RHO"         ,   1.0    ,	/* -- Density.                            */
  "GRAVITY"     ,   9.81   ,	/* -- Gravitational acceleration.         */
  "T_REF"       ,   288.15 ,	/* -- Reference temperature (15C).        */
  "PRANDTL"     ,   0.72   ,	/* -- Prandtl number for air at STP.      */

  "C_SMAG"      ,   0.1    ,	/* -- Smagorinsky's constant (RNG value). */
  "LAMBDA_M"    ,   2.0    ,    /* -- Assumed difference in mesh lengths. */

  "FFX"         ,   0.0    ,	/* -- Body force per unit mass (x).       */
  "FFY"         ,   0.0    ,	/* -- y component.                        */
  "FFZ"         ,   0.0    ,	/* -- z component.                        */
  "FFC"         ,   0.0    ,    /* -- (Analogous) scalar production rate. */

  "X_SHIFT"     ,   0.0    ,    /* -- Optional shift to mesh in x.        */
  "Y_SHIFT"     ,   0.0    ,    /* -- Optional shift to mesh in y.        */

  "X_SCALE"     ,   1.0    ,    /* -- Optional factor to scale mesh in x. */
  "Y_SCALE"     ,   1.0    ,    /* -- Optional factor to scale mesh in y. */

  "SVV_MN"      ,  -1      ,    /* -- SVV SEM Mode num (-1 is off) < N_P. */
  "SVV_EPSN"    ,   0.0    ,    /* -- SVV SEM Eps, usually propto KINVIS. */

  "SVV_MZ"      ,  -1      ,    /* -- SVV Fourier mode start, < (N_Z/2).  */
  "SVV_EPSZ"    ,   0.0    ,    /* -- SVV Fourier Eps, as for SEM.        */

  "UODelta"     ,   0.05   ,    /* -- Outflow boundary velocity scale.    */

  /* -- Option switches. */

  "ITERATIVE"   ,   0   ,	/* -- Select PCG solver for velocities.   */
  "TBCS"        ,   0   ,	/* -- Select time-varying BCs.            */
  "CYLINDRICAL" ,   0   ,	/* -- Select cylindrical coordinates.     */
  "VERBOSE"     ,   0   ,	/* -- Set verbose output.                 */
  "CHKPOINT"    ,   1   ,	/* -- Set checkpointing of field dumps.   */
  "AVERAGE"     ,   0   ,	/* -- Select averaging of fields.         */
  "SPAWN"       ,   0   ,	/* -- Set respawning for particle tracks. */
  "N_PHASE"     ,   0   ,       /* -- Toggle phase averaging.             */
  "RANSEED"     ,   0   ,       /* -- Set wall-clock random seeding.      */
  
  /* -- Default integer values. */

  "IO_FLD"      ,   500 ,	/* -- Step interval for field dumps.     */
  "IO_HIS"      ,   10  ,	/* -- Step interval for history points.  */
  "IO_CFL"      ,   50  ,	/* -- Step interval for CFL + divergence.*/
  "IO_MDL"      ,   50  ,	/* -- Step interval for modal energy.    */
  "IO_WSS"      ,   0   ,       /* -- Step interval + toggle of WSS out. */

  "N_P"         ,   5   ,	/* -- No. of points along element edge.  */
  "N_TIME"      ,   2   ,	/* -- Order of timestepping scheme.      */
  "N_STEP"      ,   1   ,	/* -- Number of timesteps to integrate.  */
  "N_Z"         ,   1   ,	/* -- Number of planes of data.          */
  "I_PROC"      ,   0   ,	/* -- Process index for parallel soln.   */
  "N_PROC"      ,   1   ,	/* -- Number of processes for parallel.  */
  "STEP_MAX"    ,   500 ,	/* -- Max number of iterations for PCG.  */
  "NR_MAX"      ,   20  ,       /* -- Max iterations for Newton-Raphson. */
  
  0             ,   0.0
};




