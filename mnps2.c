/*******************************************************************************
*     PROGRAM WildlifeDensity
*
*     (File mnps2.c, Version 2.4b2)
*
*     This program is designed to return population density estimates
*     from 'distance' data collected using either line transect or fixed
*     point methods, together with estimates of other parameters of the
*     observing situation.  It is intended primarily for use in
*     ground surveys of mammal and bird populations in terrestrial
*     habitats, and with aerial surveys using helicopters.  It is also
*     capable, with modification, of being used with aerial survey
*     data collected from fixed-wing aircraft.
*
*     WildlifeDensity works by mathematical modelling of the frequency
*     distributions of animal numbers detected under uniform observing
*     conditions over a range of distances from either an observer
*     or a transect line.  Curve-fitting is achieved by using the
*     simplex method to seek parameter values that minimize an ordinary
*     least squares fit between the observed and calculated numbers across
*     the range of detection distances.
*
*     Based on the 'best-fit' solution so obtained, the program returns
*     population density estimates, estimates of other parameters of the
*     model and related information on the system in which the original
*     observations were made. Standard errors of the parameters are
*     obtained by 'bootstrap' resampling of data from the original data
*     set, unless the number of sets of iterations is set at 1, when the
*     program computes standard errors by a quadratic surface fitting
*     method (qsf).
*
*
*     The WildlifeDensity model was designed by David Morgan, Department
*     of Zoology, The University of Melbourne, Vic. 3010, Australia, and
*     engineered by James Clough.
*
*
*     ******************************************************************
*
*
*     File mnps2.c has been written in C, and consists of a main
*     subroutine, calculate_density(), and five subsidiary subroutines:
*     givef(), resample(), srandnag(), srandgen() and qsf().  Only the
*     calculate_density() routine should be called by external programs.
*
*
*     The main subroutine determines parameter values for a least
*     squares fit between calculated and observed frequencies,
*     using the simplex method, a numerical method for minimizing an
*     objective function in a many-dimensional space.  Its design is
*     derived from Nelder & Mead (1965) Computer Journal 7:308-313, as
*     used in the 'MINIM' program by D.E.Shaw, Divn. of Mathematical
*     Statistics, CSIRO, Australia, though extensively modified since.
*     The quadratic surface fitting method within it is based on a Hessian
*     matrix, as set out in the appendix to the Nelder and Mead paper.
*
*     Subroutine GIVEF calculates the sum of squares by using the
*     mathematical model, and returns this to the main program.
*     Subroutines SRANDNAG and SRANDGEN supply a pseudo-random number used
*     in the data resampling procedure in response to a seed value supplied
*     by the main program.
*
*
*     The program computes density estimates based on the numbers of
*     animals detected either at various horizontal radial distances (r)
*     from an observer or at various horizontal distances (y) perpendicular
*     to a transect line.  It accepts several alternative types of distance
*     data, as determined by appropriate values of the control parameters
*     IFX and IRY supplied with the data, viz:
*
*       . radial distance data supplied (IFX=0 and IRY=0);
*       . perpendicular distances calculated from radial distance
*           and horizontal detection angle data supplied
*           (IFX=0 and IRY=1);
*       . perpendicular distance data supplied (IFX=0 and IRY=2);
*       . fixed-observer distance data supplied (IFX=1 and IRY=0).
*
*
*     Data are entered into the program by means of an appropriately
*     formatted tab- or comma-separated dataset in a file named
*     '<filename>.dat'.  It outputs two files:  the first, named
*     '<filename>.results', prints out final values of the parameters
*     of the dataset supplied, together with some of the main attributes
*     of the observing situation;  the second, named '<filename>.graphData',
*     sets out the original frequency distribution and the calculated
*     frequency distribution of best fit in a column format suitable
*     as input to a graphing application such as MS 'Excel'.
*
*
*     Key to Nomenclature:
*
*     The code names of the parameters supplied as input
*     via the data file are, in alphabetical order:
*
*
*     ANGLE[] - the horizontal angle between the direction of a transect
*          and the bearing to a cluster of animals at the moment of
*          detection, and measured in degrees.  These data are
*          required if computations are to be based on perpendicular
*          distance modelling; otherwise they are not needed.
*          The program accepts negative angles in the input (e.g.
*          for observation left of a transect line): these are
*          pooled with positive angles during computations.
*
*     CLINT - the class interval to be used within the program,
*          chosen so that 80 x CLINT is at least equal to the
*          maximum radial detection distance.
*
*     DIST - the overall transect length (L), expressed either in
*          metres or kilometres, or the distance corrected for the
*          effect of animal movement (LJ).  (If F[NUM_SHAPE_PARAMS-1] is set at
*          zero, only the overall transect length is needed and
*          the program calculates its own LJ value.)
*
*     DURN - the duration of a fixed-point census (in min.).
*
*     F  - the model's parameter values used in the search for the
*          minimum point.  Four initial F values are required
*          as input to the program:
*
*          F[0]:  an estimate of the conspicuousness coefficient
*                 of the species under the conditions of the
*                 census;
*
*          F[1]:  an estimate of either the overall lateral vegetation
*                 cover (c) between observer and animals, or of the
*                 attenuation coefficient (b) of animal sounds in the
*                 observing situation under the census conditions;
*
*          F[2]:  an estimate of the population density (D),
*                 expressed in number of individuals per hectare
*
*          F[3]:  an estimate of the maximum distance (in m) from the
*                 observer at which species recognition is
*                 possible under the conditions of the census.
*
*          Those supplied with the data are the starting values
*          for the search; those printed on exit are the parameter
*          values which specify the minimum point.  Initial F[]
*          values supplied are best based on existing knowledge.
*
*     FUNC - the sum of squares of the differences between
*          observed and estimated D2L values;
*
*     IFX - a parameter to modify the computations appropriately
*          where data are collected by an observer sampling from
*          a fixed point, viz:
*
*            = 0  handles data from line transects;
*            = 1  handles data from fixed point sampling.
*
*     IPRINT - a parameter which directs the program to print out
*          progress sum of squares and parameter values to allow
*          perusal of the progress of the minimization process.
*
*            = 0  no progress evaluations;
*            = 1  If JPRINT=0, reports initial convergence and the
*                 function and parameter values of the initial
*                 simplex from each bootstrap iteration series;
*                 If JPRINT=1, directs printing of these values at
*                 individual steps through the function minimization
*                 process, and informs the  user if the process
*                 failed to converge on a minimum within 750
*                 steps (the limit set).
*
*     ISHOW - a parameter which prints out the values of the
*          following functions, in sequence, for each class interval
*          across the distance range for the final result: PRC,
*          PRR, QR, E, TOTE (total E), ED, EXPD, EXPDV and the
*          detectability coefficient S (see comments within
*          Subroutine GIVEF.
*
*            = 0  no output;
*            = 1  produces output.
*
*     IRY - a parameter to control computations to model different
*           types of line transect data, viz:
*
*            = 0  models radial distance transect data N(r);
*            = 1  models perpendicular distance data N(y), using
*                   radial distance and angle data supplied;
*            = 2  models perpendicular distance data N(y), using
*                    perpendicular distance data supplied.
*
*     JPRINT - a parameter to enable perusal of the original and
*           bootstrapped frequency distributions, viz:
*
*            = 0  no frequency distributions or internal output;
*            = 1  prints the frequency distributions of the original
*                 and bootstrapped data and makes possible perusal
*                 of some of the function minimization process.
*
*     KDT - a control parameter to indicate attributes of the
*          observational data supplied to the program, viz:
*
*            = 0  for visual and flushing data;
*            = 1  for auditory data or long distance (no vegn.) data;
*            > 1  the maximum class interval boundary distance (in m)
*                 when visual data are from a limited distance range.
*                 The unit is still m even if KM is set.
*
*     KM - a parameter to modify the program if transect lengths and
*          detection distances are measured in kilometres, viz:
*
*            = 0  all distances are measured in metres;
*            = 1  transect lengths are measured in kilometres
*                   and detection distances in metres;
*            = 2  both transect lengths and detection distances
*                 are measured in kilometres.
*
*     LTMIN - the approximate minimum distance at which an animal
*          may be obscured from the observer by topography. Set
*          at 999 if the topography is approximately level.
*
*     LTMAX - the p=0.001 upper-limit distance at which an animal
*          may be detected by the observer in uneven topography
*          in the absence of vegetation cover.  Set at 0 if this
*          distance is to be calculated from the data.
*
*     MAXJB  -  the number of sets of iterations used to derive
*          the backstrapped series of model parameters and
*          calculate standard errors.  This is best set at between
*          500 and 2000.  Fewer sets saves computation time but
*          increases the size of confidence intervals of parameters.
*          If MAXJB is set at 1, no bootstrapping occurs; instead
*          the program calls Subroutine qsf and calculates variances
*          and standard errors using surface fitting algorithms.
*
*     NS (& SNS) - the number of sides of the transect line used for
*          observations during a line transect survey, viz:
*
*            = 0  not transect data (fixed point), so not relevant;
*            = 1  handles data from only one side of the line;
*            = 2  handles data from both sides of the transect line.
*
*     NSIZE[] - the size of each cluster (group) of animals at the
*          moment of detection, listed in the data input in the same
*          order as the corresponding R[] values.
*
*     NUMA (& NUMAIN) - the total number of individual animals detected
*          ahead of a moving observer on a transect.
*
*     NUMO (& NUMOIN) - the number of individuals overtaking the observer
*          from behind during a transect, obtained by totalling
*          the number of individuals in R=0 observations.
*
*     NVALS - the total number of individual animal detections
*          in the sample submitted for modelling.
*
*     PD - the proportion of the population observable at the
*          time of a census, usually taken =1.
*
*     PS - in fixed-point censuses, the proportion of the circle
*          scanned by the observer during the census (usually
*          either 0.5 or 1).
*
*     R[] - the radial detection distances (r) originally measured
*          in the field, submitted in the order in which the data
*          were collected.  Or the perpendicular distances from the
*          transect line if these are precalculated and entered directly.
*
*     RATE - the overall mean rate of animal movement (in m/min)
*
*     STEP - the step sizes used initially in the function
*          minimization process to modify the F values. Typical
*          STEP values are:
*
*          STEP[0]:  just under 50% of the initial F[0] value;
*
*          STEP[1]:  equal to the initial F[1] value;
*
*          STEP[2]:  half the initial F[2] value;
*
*          STEP[3]:  usually set at zero.
*
*          Setting the STEP value at 0.0 for a parameter fixes an F[]
*          value at that supplied to the program.  Where a data set
*          (NVALS) is not great (say, < 250 detections), STEP[1] is
*          set at 0.0 by the program.  Where the set is very small
*          (say, < 30 detections), STEP[0] should also be set at 0.0 .
*
*     STT - the detection distance value (r or y) from which
*          class intervals begin (usually zero). Entered in m even if
*          KM is set.
*
*     THH - the vertical distance between observer eye level
*          and the median horizontal plane occupied by
*          the population, best represented by the root mean
*          square of the measured vertical height differences.
*
*
*
*    The names of the temporary variables used within the main subroutine
*    are as follows, essentially in alphabetical order:
*
*
*     A, B, C - Three coefficients used in the search for a minimum
*          value, viz:
*
*          A - a reflection coefficient;
*          B - a contraction coefficient; and
*          C - an expansion coefficient.
*
*     BOOTSTRAP - an index variable used in DO loops in the program.
*
*     COEFF[], DEN[] - intermediate values of the parameter estimates.
*
*     COEFFNT1 - estimated value of the conspicuousness coefficient.
*
*     COEFFNT2 - estimated value of the attenuation coeffiocient or lateral
*          cover.
*
*     COEFFNT3 - estimated value of the detectability coefficient.
*
*     CF1DIF, CF2DIF, CF3DIF, CF1SUM, CF2SUM, CF3SUM - temporary variables
*          used to calculate standard errors of the best-fit parameters.
*
*     D2L - the product of the population density and twice the
*          total transect length - a convenient variable in
*          calculations.
*
*     DCOEFF - estimated value of the intermediate coefficient D2L.
*
*     DENDIF, DSUM - temporary variables used in the calculation of
*          standard errors of the best-fit parameters.
*
*     DMAX - the maximum direct-line detection distance ('dmax'),
*          either submitted to the program as F[NUM_SHAPE_PARAMS-1] or calculated.
*
*     ERMAX - a computational upper-limit maximum detection distance
*          used as a range limit in certain cases.
*
*     ESTDEN  - the estimated population density ('D').
*
*     ESTJ - an overall movement correction factor (J) for line
*          transects, calculated when transect duration and animal
*          movement rate data are entered, using an approximation for
*          the J[k] range 0 < k < 5.
*
*     ESTDMAX - the value of DMAX (= F[3]) calculated from the radial
*          distance data in the program input if F[3] is set at zero.
*
*     FNK - the ratio k=u/w, used to estimate a movement correction factor.
*
*     G[] - initial values of the estimated parameters used to derive
*           simplexes.
*
*     H, H[], HMIN, HMAX, HSTAR, HSTST - sum of squares totals used in
*           the search for a function minimum.
*
*     HSTD - a temporary variable used in estimating the
*          standard errors of function values in the current simplex.
*
*     I, IA, IC, IE, IFLAG, IH, IN, IR, J, JS, JV, K - index
*          variables used in DO loops within the program.
*
*     LOOP - a parameter used to control the search process.
*
*     MAX - the maximum number of function evaluations to be
*          allowed (arbitrarily set at 750);
*
*     MFAIL, MSFAIL - two variables used to count the numbers of
*          times when convergence fails (MFAIL) or when the
*          detectability coefficient (S) cannot be calculated (MSFAIL).
*          They are used to correct the numbers of counts used in
*          computing means and standard deviations of key parameters.
*
*     MTEST - a variable used to flag a successful run through Loop 1410,
*          and assist in the computation of MSFAIL at the end of
*          Subroutine GIVEF.
*
*     NAP - the number of parameters to be varied with STEP not equal
*          to zero.
*
*     NBSZ[] - numbers of animals in a cluster produced by the
*          bootstrapping process.
*
*     NCLASS - the number of distances classes in the range used for the
*          analysis.  It is either equal to DMAX-STT, or KDT-STT,
*          depending on whether or not KDT is greater than 1.
*
*     NEVAL - the number of evaluations used to date in the search for
*          a minimum difference.
*
*     NGROUPS - the total number of observations within the selected
*          range of distances (= NVALS - NOTIN).
*
*     NLOOP - convergence is tested for every NLOOP times the
*          process changes the simplex.  After initial
*          convergence, NLOOP further changes are allowed
*          before testing for final convergence.
*
*     NOP - the number of parameters potentially varied during
*          the minimization process.
*
*     NOTIN - a count variable that records the number of occasions
*          in which an observation falls outside a selected data
*          range, i.e. where distance is less than STT or more than KDT.
*
*     NOVTKS - the number of groups overtaking the observer during a
*          line transect.
*
*     NP1 - NAP incremented by 1.
*
*     NUMGRA - the number of individuals or groups detected ahead of
*          observer during a transect.
*
*     OBSW - the overall observer rate of travel along the line
*          transects, calculated from the total distance travelled and
*          the time taken, and expressed in m/min.
*
*     P  - the conspicuousness coefficient ('a').
*
*     PBAR[], PSTAR[], PSTST[] - temporary function values used in simplexes
*          in the search for a function minimum.
*
*     Q  - for visual data (where KDT=0 or >1), the mean vegetation
*          cover proportion (c) in the habitat between animal and
*          observer.
*        - for auditory data (where KDT=1), the conspicuousness
*          coefficient ('b').
*
*     RESAMP_DIST[] - a value of either R[] or Y[] reselected at
*          random from the original data as the key part of the
*          bootstrapping process.
*
*     RLMEAN - the mean of the natural logarithms of radial distances.
*
*     RLSD - the standard deviation of the natural logarithms of the
*          radial distances in the data set.
*
*     RMAX - the maximum horizontal detection distance ('rmax'),
*          calculated from DMAX and THH.
*
*     S  - a coefficient of detectability, usable in density estimation
*          by direct calculation.
*
*     SAVEMN - a temporary variable used in estimating the
*          standard errors of function values in the current simplex.
*
*     SCF1, SCF2, SCF3 - temporary variables used in the calculation of
*          standard errors of the best-fit parameters.
*
*     SDEN - the estimated standard error of the 'best-fit' density
*          estimate.
*
*     STOPC - a stopping criterion used to initiate a test for
*          convergence on a minimum point.
*
*     T, TDSUM, TTRDENTIN, TDENDIF, TCL1, TCL2, CL1, CL2, STRDEN -
*          intermediate variables used during transformations in
*          estimating confidence limits.
*
*     T001, RLTST, RLSUM, RDIFSQ, RLF4, RDIFSQ - intermediate variables
*          used in estimating a maximum detection distance from data
*          submitted to the program.
*
*     TCOV - a topographical cover value, i.e. the estimated proportion
*          of topographical cover in the line of sight between observer
*          and animal.
*
*     TDEN, TCOEFF1, TCOEFF2, TCOEFF3, TRDEN - temporary variables used
*          in accumulating parameter estimates.
*
*     TEST - an intermediate variable that tests the closeness of the
*          ratio SAVEMN/HMEAN to 1.
*
*     TOPDD - the proportion of the population at a given radial
*          distance that is unobscured by topographical cover (such as
*          hills and ridges).
*
*     useMedianValues (IMV) - a parameter to select the method
*          of calculating probability used in Subroutine GIVEF.
*          It either:
*
*            = false uses the mean value of P[r] in an interval;
*            = true  uses the median value of P[r] in the interval.
*
*          useMedianValues is enabled within the program if the range of
*          observations is restricted (KDT > 1) or if the coefficient
*          becomes negative. (The two approaches produce almost identical
*          results.)
*
*     VALT[] - the initial class totals.
*
*
*     Temporary variables used exclusively in the other subroutines are
*     listed in the subroutine statements, and the relevant explanatory
*     notes are included within the program itself.
*
*     This program uses double precision for most real numbers.
*
*	Translated from Fortran by Jeremy Begg, VSM Software Services Pty Ltd
*	The '//' comment syntax is used to highlight the original Fortran code
*	where it was deemed appropriate.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

#include "mnps2.h"

FILE *output_results = NULL;


/*
*  The following output construct is used at several places in the program
*  Assumes NUM_SHAPE_PARAMS is 4 (same as the 'nop' variable).
*/

#define DUMP_SS(ss,p) fprintf(output_results,"\n   %4i    %13.6e %13.6e %13.6e %13.6e %13.6e\n",neval,ss,p[0],p[1],p[2],p[3]);
#if (NUM_SHAPE_PARAMS != 4)
#error "Check DUMP_SS definition!"
#endif


/*******************************************************************************
*
*	SRANDNAG	Initialise the random number generator
*	RANDNAG		Generate a random number
*
*******************************************************************************/

/*
*   Static local storage -- these variables are not visible to calling programs
*
** JMB: In FORTRAN these were declared as follows in each routine:
**
**	DOUBLE PRECISION jseed,ifrst
**	COMMON /seed/ jseed,ifrst
**
**     plus a BLOCK DATA definition as follows:
**
**	BLOCK DATA randbd
**	DOUBLE PRECISION jseed,ifrst
**	COMMON /seed/ jseed,ifrst
**
**	DATA jseed,ifrst/123456789,0/
**
**	END BLOCK DATA randbd
*
*/

static double jseed=123456789.0;
static int ifrst=0;


/*
*	srandnag (iseed)
*
*  This subroutine sets the integer seed to be used with the
*  companion RANDNAG function to the value of ISEED.  A flag is
*  set to indicate that the sequence of pseudo-random numbers
*  for the specified seed should start from the beginning.
*
*/

void srandnag (int iseed)
{
      /* printf("Srandnag(%d)\n", iseed); */
      jseed = iseed;
      ifrst = 0;

      return;
}
//      END SUBROUTINE srandnag


/*
*	randgen (result)
*
*  This routine returns a pseudo-random number for each invocation.
*  It is an adaptation of the "Integer Version 2" minimal standard
*  number generator whose Pascal code appears in the article:
*     Park, Steven K. and Miller, Keith W., "Random Number
*     Generators: Good Ones are Hard to Find", Communications of
*     the ACM, October, 1988.
*  It is used in the form recommended by the Numerical
*  Algorithms Group (NAG) for UNIX systems.
*
*/

void randgen (double *result)
{
      const int mplier=16807;
      const int modlus=2147483647;
      const int mobymp=127773;
      const int momdmp=2836;

      int hvlue, lvlue, testv;
      static int nextn;


      if (ifrst == 0) {
        nextn = jseed;
        ifrst = 1;
      }

      hvlue = nextn / mobymp;
      lvlue = nextn % mobymp;
      testv = mplier*lvlue - momdmp*hvlue;
      nextn = (testv > 0) ? testv : testv + modlus;

      *result = (nextn * 1.0) / (modlus * 1.0);
}
//      END SUBROUTINE randgen


/*******************************************************************************
 *
 *	FREQ_DISTRIB     Calculate frequencies of each class
 *
 *       Inputs: NCLASS, STT, CLINT, NVALS, KDT - as described above
 *               Y  - distance to observed groups, whether radial or perpendicular,
 *                    original or resampled
 *               ABSOLUTE_DISTANCES  - if true, consider negative distances as positive
 *                                     and consider zero distances to be in the first class
 *               NSIZE  - size of the observed groups, original or resampled
 *
 *       Outputs: NUMA, VAL, NGROUPS - as described above
 *                NUMO  - counted if ABSOLUTE_DISTANCES is false
 *
 *******************************************************************************/

void freq_distrib(int nclass, double stt, double clint, int nvals, int kdt, float y[],
                  bool absolute_distances, int nsize[],
                  int *numa, int *numo, double val[], int *ngroups)
{
    *numo = 0;
    *numa = 0;
    *ngroups = 0;
    for (int ic=0; ic < nclass; ic++) {			//  DO ic=1,nclass
        val[ic]=0;
    }

    double max_distance = (nclass * clint) + stt;
    
    /*
     *      Numbers ahead (NUMA), overtaking (NUMO) and groups are totalled.
     *      Numbers in each class, VAL[IC], are also accumulated for the groups
     *      included in the data set.
     */
    for (int irb=0; irb < nvals; irb++) {	//  DO irb=1,nvals

        float group_distance = y[irb];
        if (absolute_distances) {
            /*
             * When data from the two sides of a transect line are pooled,
             * absolute values of Y() are used.
             */
            group_distance = fabs(y[irb]);
        }
        
        if ((kdt > 1) && (group_distance > kdt)) {
            // group distance is beyond maximum detection distance
        }
        else if (group_distance < stt) {
            // group distance is before the first class (only possible without absolute_distances)
        }
        else if (group_distance >= max_distance) {
            // group distance is beyond the final class
        }
        else if (!absolute_distances && (group_distance == 0.)) {
            *numo += nsize[irb];
            (*ngroups)++;
        }
        else {
            // class IC includes groups at clint x IC, excludes groups at clint x (IC+1)
            int ic = (group_distance - stt) / clint;
            *numa += nsize[irb];
            (*ngroups)++;
            val[ic] += nsize[irb];
        }
    }					//  END DO irb
}


/*******************************************************************************
 *
 *       SELECT_SEARCH_PARAMETERS     Calculate parameters to use for the density calculation
 *
 *       Inputs: most values from the GUI are considered
 *
 *       Outputs: NUMA, VAL, NCLASS, NGROUPS, STT, CLINT, ESTJ - as described above
 *                NUMO     - counted if ABSOLUTE_DISTANCES is false
 *                DISTANCE - distance to observed groups, corrected to metres if entered in km
 *                F, STEP  - start point and step size for the search
 *
 *******************************************************************************/

void select_search_parameters(calc_params *params, search_params *result)
{
    result->f[0] = params->enteredValue[0];
    result->f[1] = params->enteredValue[1];
    result->f[2] = params->enteredValue[2];
    /*
     *     If no value of the maximum detection distance F[3] has been
     *     entered (i.e. f[3]=0), and radial detection distances are
     *     provided in the data, an estimated maximum distance is
     *     calculated based on the assumption that the logarithm of
     *     the radial detection distance r is distributed according
     *     to a normal distribution, with a maximum value at the
     *     mean + (t.001 x s.d.), which matches observed values from
     *     large data sets.  This is calculated, then a back-transformation
     *     undertaken to give an F[3] value.
     *
     *     If f[3]=0 and only perpendicular distance data are supplied,
     *     the logarithms of the calculated perpendicular distances are
     *     assumed to follow a half-normal distribution (a special case of
     *     the folded normal distribution with a mean of zero). Its maximum
     *     value is then estimated as (t.001 x s.d.), where s.d. is the
     *     standard deviation of the half-normal distribution.  If only
     *     perpendicular distances are supplied AND a maximum class boundary
     *     distance has been set which is within the observed distribution,
     *     f[3] should NOT be set at 0 but given an approximate (maximum)
     *     value instead.
     *
     */

    if (params->enteredValue[3] == 0) {
        int numgra = 0;
        float rltot = 0, rlsum = 0;
        double t001 = 3.33256 + 33.0731/pow(params->nvals, 1.39337);

        if (params->iry < 2) {

            /*
             *     This option handles all situations with radial data supplied.
             *     The first step is to calculate a logarithmic detection distance
             *     total RLTOT, then a logarithmic mean value RLMEAN, using data
             *     from ahead of the observer only.
             */

            for (int ih=0; ih < params->nvals; ih++) {
                if (params->r[ih] > 0) {
                    rltot += log(params->r[ih]+1);
                    numgra++;
                }
            }

            float rlmean = rltot/numgra;

            /*
             *     A standard error of the logarithmic r (RLSD) is now calculated,
             *     followed by the estimated logarithmic maximum distance RLF4,
             *     which is then backtransformed to give an F[3] value.
             */

            for (int ih=0; ih < params->nvals; ih++) {
                if (params->r[ih] > 0) {
                    float tmp = log(params->r[ih] + 1) - rlmean;
                    float rdifsq = (tmp * tmp) / (numgra - 1);
                    rlsum += rdifsq;
                }
            }

            float rlsd = sqrt(rlsum);
            float rlf4 = rlmean + t001*rlsd;

            result->f[3] = exp(rlf4) - 1;
        }

        else {

            /*
             *   The following option is used only if recalculated perp. data have been
             *   supplied to the program (iry == 2).  It involves normalising their
             *   distribution by square root transformation, estumating the standard
             *   deviation of the half-normal distribution, with Bessel's correction,
             *   multiplying it by the 0.001 value of Student's t, then back-transforming
             *   by squaring.  The first step is to sum the absolute perpendicular distance
             *   values in Loop_43, then calculate the Bessel-corrected mean of these
             *   distances, and finally estimate the maximum detection distance by squaring
             *   the result.
             */

            for (int ih=0; ih < params->nvals; ih++) {
                rlsum += fabsf((float) params->r[ih]);
            }

            float rlsd = sqrt(rlsum/(params->nvals-1));
            float rlf4 = rlsd*t001;
            result->f[3] = rlf4 * rlf4;
        }
        
    }  /* end if (f[3] == 0) */

    else {
        result->f[3] = params->enteredValue[3];
    }

    for (int ih = 0; ih < params->nvals; ih++) {
        result->distance[ih] = params->r[ih];
    }
    if (params->ifx == 0 && params->km > 0) {
        /* Transect lengths have been expressed in kilometres. Convert to metres */

        // dist *= 1000;
        if (params->km == 2) {
            result->f[3] *= 1000;
            /*
             *     If detection distances were entered in kilometres (km=2),
             *     then these distances are also converted to metres.
             */

            for (int ih = 0; ih < params->nvals; ih++) {
                result->distance[ih] *= 1000;
            }
        }
    }
    
    double min_distance = params->kdt > 1 ? params->stt : 0;
    double max_distance = params->kdt > 1 ? fmin(params->kdt, result->f[3]) : result->f[3];
    double distance_range = max_distance - min_distance;
    result->stt = min_distance;

   /*
    *      Compute the class interval.
    */

    if (params->clint > 0) {
        // To do: manual CLINT values
    }

    int ngroups = 0;
    double nearest_observation = max_distance;
    double furthest_observation = min_distance;
    for (int ih = 0; ih < params->nvals; ih++) {
        if (result->distance[ih] != 0.) {
            ngroups++;
            nearest_observation = fmin(nearest_observation, fabs(result->distance[ih]));
            furthest_observation = fmax(furthest_observation, fabs(result->distance[ih]));
        }
    }
    printf("select_clint: [%.1f - %.1f) => [%.1f - %.1f)\n", min_distance, max_distance, nearest_observation, furthest_observation);
    if (furthest_observation < nearest_observation + 0.06) {
        // nonsensical data
        result->clint = 1;
    } else {
        double num_intervals = fmax(5, fmin(20, ngroups / 5.33));
        result->clint = ceil((furthest_observation - nearest_observation + 0.1) / num_intervals);
    }

    /*
     *     If a value of either the maximum detection distance F[3]
     *     or the maximum of the selected interval KDT has been
     *     entered as more than 80 times the class interval, CLINT is
     *     reset at (F[3]-STT or KDT-STT)/80 to avoid computation problems.
     */
    
    result->clint = fmax(result->clint, distance_range / MAX_INTERVALS);

    /*
     *     Case_150 treats the situation where calculations are to be based
     *     on perpendicular distances (y) from the transect line, and the
     *     data supplied are either radial distances and angles, indicated
     *     by IRY=1, or pre-calculated perpendicular distances, indicated
     *     by IRY=2.  If IRY=2, distances entered as r(in) are reassigned
     *     as y(in) values.  If IRY=1, perpendicular distances are
     *     calculated from distances and lateral angles using trigonometry.
     *     With both situations, any angle data supplied as negative numbers
     *     (e.g. those to the left of a transect line) are converted to
     *     positive and pooled with the remainder.
     *     If calculations are to be based on radial detection distances,
     *     and radial distances only are supplied (IFX=0), no changes are
     *     made (either by calculation or reassignment).
     *     If a radial distance value of precisely zero (r==0) was supplied
     *     to the program in the Observations, then each such value will be
     *     recognized as an ‘overtake’ at a later point in the program.
     */
    
Case_150:
    switch(params->iry) {
            
        case 2:
            /*
             *     If perp. distance data were entered as r values (IRY=2), they
             *     are renamed as y values at this stage, unless r=0
             *     when they are 'overtakes' and omitted.  Negative y values
             *     submitted to the program (as negative r value) are
             *     converted to positive and pooled with the rest.
             */
Loop_190:
            for (int in = 0; in < params->nvals; in++) {
                result->distance[in] = fabsf(result->distance[in]);
            }
            break;
            
            
        case 1:
            /*
             *     For IRY=1, perpendicular distances are calculated from radial
             *     distances and lateral observing angles.
             */
Loop_150:
            for (int in = 0; in < params->nvals; in++) {
                result->distance[in] = fabs(result->distance[in] * sin((params->angle[in] * 3.14159265) / 180.));
            }
            break;


        default:
            /*
             *     If IFX=0, no reassignment is needed.
             */
            break;
            
    } /* switch */
    
    /*
     *     NCLASS is the number of distance classes in the selected range.
     *     0.49 is added to avoid counting errors due to 'chopping'.
     */

    result->nclass = ((result->f[3] - result->stt) / result->clint) + 0.49;
    if (result->nclass > MAX_INTERVALS)
        result->nclass = MAX_INTERVALS;

    /*
     *     The program now computes the first distribution of the numbers
     *     detected within each class interval, based on the detection
     *     distances (R[IN]), the numbers in each group (NSIZE[IN]) and
     *     the class interval (CLINT) preset in the input.
     *
     *     In subsequent runs through Loop 1410, bootstrapping applies
     *     (JBSTP=1) and a bootstrapped distribution is used instead.
     *
     *     The number detected within each class (VAL[IC]), is the sum of the
     *     numbers in each class in which the R[IN] or Y[IN] values fall.
     *     Calculating the various VAL[IC] values first requires
     *     finding which R[IN] or Y[IN] values fall within the interval
     *     concerned, then adding all the NSIZE[IN] values which fall
     *     within that class.  This will be done for each class interval
     *     in turn, beginning with the calculation of VAL[0] for the
     *     nearest class to r=0 or STT or y=0 or STT.  If a minimum value in the
     *     range (STT) has been specified, or a maximum value for r or y of KDT
     *     (>1), the program also computes the number of data clusters (NOTIN)
     *     below STT and above KDT and subtracts it from NVALS to give the
     *     correct magnitude of NVALS for use in later calculations.
     */
    
    freq_distrib(result->nclass, result->stt, result->clint, params->nvals, params->kdt, result->distance, false,
                 params->nsize, &result->numa, &result->numo, result->val, &result->ngroups);

    result->step[0] = params->enteredStep[0];
    result->step[1] = params->enteredStep[1];
    result->step[2] = params->enteredStep[2];

    if (params->ifx == 0) {
        /*
         *     The program now calculates the mean overall observer movement
         *     rate (w) as OBSW=DIST/DURN (DIST in m and DURN in min), provided
         *     that an DURN value has been entered and the data are from line
         *     transects (IFX=0).  Also, k is defined as k=u/w (FNK=RATE/OBSW).
         *     If not, this entire step is bypassed.
         *
         * Comment by J.Begg: in FORTRAN, dividing by 0 gives Infinity but this
         * doesn't seem to be the case in C.  So the logic is extended to include
         * a test for DURN=0 and set FNK accordingly.  (Note that OBSW is not used
         * for any purpose other than to calculate FNK.)
         */
        double fnk;
        if (params->durn <= 0) {
            fnk = 0;
        }
        else {
            double obsw = params->dist/params->durn;
            fnk = params->rate/obsw;
        }
        
        /*
         *     Assuming that the actual distance (and not LJ) has been
         *     entered, the movement correction factor (J) is calculated
         *     using an approximation.  Different approximations are used
         *     if k=u/w is less than or greater than 5.
         */
        if (fnk == 0) {
            result->estj = 1;
        } else if ((fnk > 0) && (fnk <= 1)) {
            result->estj = 1.000 + (0.00495562*fnk) + (0.0995941*fnk*fnk) + (0.0324447*fnk*fnk*fnk);
        } else if ((fnk > 1) && (fnk <= 5)) {
            result->estj = 1.0051 - (0.0212978*fnk) + (0.0002868*fnk*fnk)  + (0.279106*fnk*fnk*fnk)
            - (0.12982*fnk*fnk*fnk*fnk) + (0.0234994*fnk*fnk*fnk*fnk*fnk) - (0.00153274*fnk*fnk*fnk*fnk*fnk*fnk) ;
        } else {
            result->estj = 0.8183*fnk;
        }
    }

    /*
     *     Unless the number of iterations has been set at 1 and bottom option 3 has not
     *     been selected, the Line_100 sequence computes revised initial estimates of the
     *     parameters ‘a’‘c’ and ‘D’, together with initial step sizes for them,
     *     using a detectability coefficient estimate (ests), whenever the initial step
     *     size for either or both of f[0] or f[1] is set at zero or all step sizes
     *     are set at zero.
     */

/* Line_100: */
    if  ((params->maxjb > 1) && (params->ishow == 0))      {
        double estdmax = result->f[3];

        /* True if either or both of the initial step sizes are zero, so also
         * true if all three initial step sizes are zero */
        bool stepSizeIsZero = ((params->enteredStep[0] == 0) || (params->enteredStep[1] == 0));
        if ((params->nvals < 250) || stepSizeIsZero) {
            result->f[0] = 1.39103 * pow(estdmax, 0.334858);
            result->f[1] = 5.45005 * pow(estdmax, -0.958841);
        }

        /*
         *     In these modes the conspicuousness coefficient (‘a’, f[0]) is always searched
         *     with a preset step size based on estdmax.
         */
        result->step[0] = (0.3 * result->f[0]) ;

        /*
         *     The lateral cover proportion (‘c’, f[1]) is searched only if the sample size
         *     is large enough.
         */
        if (params->nvals >= 250) {
            result->step[1] = result->f[1] ;
        }
        else {
            result->step[1] = 0.0 ;
        }

        /*
         *     Computation of an initial values for f[2] and step[2] depends on whether
         *     line transect data (ifx=0) or fixed point (ifx=1) data are provided.
         *
         */
        double ests = 3.83527 + (0.0982708 + 5.76169e-6 * estdmax) * estdmax;
        if (params->ifx == 0) {
            result->f[2] = (1.e4 * (result->numa + result->numo)) / (params->ns * params->dist * result->estj * params->pd * ests);
        }
        else {
            result->f[2] = (1.e4 * (result->numa + result->numo)) / (2 * params->rate * params->durn * params->ps * params->pd * ests);
        }
        result->step[2] = 0.5 * result->f[2];
    }
}
//      END SUBROUTINE select_search_parameters


/*******************************************************************************
*
*	RESAMPLE		Calculate frequencies across classes
*
*******************************************************************************/

void resample (float orig_dist[], int orig_size[], int nvals,
	       float resamp_dist[], int resamp_size[])
{
      /* printf("\n");
         printf("Resampling ...\n"); */
      for (int js=0; js < nvals; js++) {	//  DO js=1,nvals
        resamp_size[js] = 0;
      }
      for (int ix=0; ix < nvals; ix++) {	//    DO ix = 1, nvals
        double x;
	randgen(&x);
	int n = x * nvals;
	resamp_dist[ix] = orig_dist[n];
	resamp_size[ix] = orig_size[n];
	/* printf("resamp(%4d) = orig(%4d)\n", ix+1, n+1); */
      }					//  END DO
}
//      END SUBROUTINE resample


/*******************************************************************************
*
*	DETECTION_PROBABILITY		Compute the estimated probability curve
*                                       for detections
*
*       Inputs: Q - The mean vegetation cover proportion (c)
*               D - Distance to detection under consideration
*               PA - The square of the conspicuousness coefficient (a)
*
*       Outputs: YY - The area under the detectability curve from the current
*                     D value to infinity
*
*******************************************************************************/
/*
*     This method of computation uses one of two approximations to
*     the exponential integral (APREXI).  If bd < 1, Approximation
*     2 is used; if bd is equal to or greater than 1, Approximation
*     1 is used.
*/

double detection_probability(double q, double d, double pa) {
/*
*     So the logarithm in APREXI doesn't get "too negative", safe_d is
*     restricited to be at least 0.03.
*/
    double safe_d = fmax(d, 0.03);
    double z = q * safe_d;
    z = fmin(68, z);
    assert(q >= 0 && z >= 0);
    double ss = exp(z);
    double y = pa / (safe_d * ss);
    double yy;

    if (z <= 0.0) {
        yy = y; // yy = y * (1 - vow) with vow=0;
    }
/*
*     The choice of approximation depends on the value of bd.
*/
    else if (z < 1.0) {
/*
*     Approximation 2 is used if bd<1.
*/
        double aprexi = -0.577216 + 0.999992*z - 0.249910*z*z + 0.055200*z*z*z
                        -0.009760*z*z*z*z + 0.0010792*z*z*z*z*z - log(z);
        yy = y - pa * q * aprexi;
    }
    else {
/*
*     Approximation 1 is used if bd=1 or bd>1.
*/
        double v = z*z + 2.334733*z + 0.250621;
        double w = z*z + 3.330457*z + 1.681534;
        double vow = v / w;
        yy = y * (1.0 - vow);
    }

    return yy;
}
//      END SUBROUTINE detection_probability


/*******************************************************************************
*
*	GIVEF		Calculate frequencies across classes
*
*******************************************************************************/

/*
*     This version of Subroutine GIVEF will handle ground survey
*     data and also aerial survey data for which there is complete
*     visibility ahead of the observer.  The subroutine
*     handles data divided into up to 80 classes, beginning its
*     computations at radial distances which exceed the maximum
*     recognition distance, and working inwards.
*
*     The subroutine calculates an 'expected' frequency in each
*     class, using parameter values supplied from the main
*     program.  It then compares each with the relevant
*     observed value supplied to the program, computes a sum
*     of squares of the differences, and returns this to the main
*     program as FUNC.
*
*/

void givef (double f[NUM_SHAPE_PARAMS],
	    double *func,
	    double *s,
	    double val[80],
	    double clint,
	    /* double pd, */
	    double stt,
	    double tcov,
	    double thh,
	    /* double vgh, */
	    int ifx,
	    int iry,
	    int ishow,
	    int kdt,
	    int kprint,
	    /* double dmax, */
	    double ltmax,
	    double ltmin,
	    int nclass,
	    int numa,
	    int numo,
	    int *msfail,
	    int maxjb,
	    int mtest,
	    bool iqsf,
            calc_results *results)
{
    if ((kprint > 0) && (ishow > 0)) {
        if (iry > 0) {
            fprintf(output_results, "\n Columns in order: P(r), Q(r), pdf, E{N(r)}, E{N(y)}\n\n");
        } else {
            fprintf(output_results, "\n\n Columns, in order: P(r), Q(r), pdf, E{N(r)}, E{N(y)}\n\n");
        }
    }

/*
*     Double precision for all real numbers is set, together with
*     common values, dimensions and symbols for key variables.
*
*/
      int iermax, jl;
      int l10,l20;

      double cint,d2l,dd,dds,ddsm,dh,dif;
      double difsq,dint,dmax,dnr,dnrl,dnrh,e,ed;
      double corrn,ermax,expd,expdv;
      double hcint,htot,obsd,p,pa;
      double pr,prc,prr,prmax,q,qr;
      double rlow,rmax;
      double texpd, tot,tote,tr, wl;

      bool useMedianValues, rlow_found = false;

/*
*     TOT is the progressive value of the sum of squares of the
*     differences between observed and calculated values.  It is
*     initially set at zero.
*/
      tot = 0.0;

/*
*     The program now sets upper and lower limits to the
*     parameter values supplied, usually by effectively giving
*     FUNC (through HTOT) a very high value if a parameter falls
*     outside predetermined limits.
*/
      p = f[0];

/*
*     HTOT is a variable used to deflect the search away from
*     highly improbable values of the parameters; HTOT=1 000 000.
*/
      htot = 0.0;
      dmax = f[3];


      q = f[1];

/*
*     Set mean values to be used for the probability distribution unless
*     the data range is not truncated or if logarithms of negative values
*     would appear in Approximation 1 above due to a negative value of Q.
*
*     The program now sets IMV=1 as a default value, unless the
*     expected frequency distribution is likely to have a sharp peak,
*     which may happen if the data are radial, NCLASS has a low
*     value (say, <20) and the data range is not truncated (KDT is not
*     greater than 1).
*/

      useMedianValues = (q < 0) || (kdt > 1);//|| (iry != 0) || (nclass >= 20);

/*
*     If Q is negative, the program sets useMedianValues and so
*     uses the median value of d in an interval as the basis of
*     computations, in order to avoid logarithms of negative
*     values appearing in Approximation 1 above.
*/
    
    if (q < 0) useMedianValues = true;

/*
*
*     F[2] is renamed D2L for its run through Subroutine GIVEF.
*/

      d2l = f[2];

/*
*     HTOT is set to a higher value if D2L<0.001 and MAXJB is not 1
*/

      if ((d2l < 0.001) && !iqsf) htot += 1.0e+6;


/*
*     HTOT is set high if a>400 or <0.01
*/
      if (((p > 400) || (p < 0.01)) && !iqsf) htot += 1.0e+6;


/*
*     RMAX - the maximum horizontal recognition distance - is
*     computed from the direct line distance and height difference,
*     unless DMAX is less than or equal to THH.
*/

      rmax = (dmax <= thh) ? 0.0 : sqrt(dmax*dmax - thh*thh);


/*
*     A range of 0 - 0.3 is set for the attenuation coefficient
*/

      if (((q >= 0.3) || (q < 0.)) && !iqsf) htot += 1.0e+6;


/*
*
*     The theoretical probability of detecting an individual at
*     the maximum recognition distance [PRMAX=Pr(rmax)] is now
*     calculated.
*/

/*
*     PA is the square of the conspicuousness coefficient.
*
*/
      pa = p*p;

    if (kdt != 1) {

/*
*     Visibility and audibility will be affected by topographical
*     features in habitats where the ground is not level. The total
*     probability of visibility at DMAX (VISMAX) will be the product of
*     the probability VEGMAX that the animal is unobscured by vegetation
*     at distance DMAX, and the probability TOPMAX that it is unobscured
*     by topography.  VEGMAX will be a function of dMAX (VEGMAX=(1-Q)**DMAX),
*     while TOPMAX is approximated by the function
*     TOPMAX=(1-T)**(LTMAX-DMAX)=(EXP(ln(0.001)/(LTMAX-LTMIN))**(LTMAX-DMAX),
*     and VISMAX=VEGMAX*TOPMAX.  A first step is to calculate TOPMAX,
*     provided that either LTMIN is not very large or RMAX is currently
*     less than LTMIN.
*/
      double topmax;
      if ((ltmin >= 999.) && (rmax <= ltmin)) {
        topmax = 1.0;
      } else {
        topmax = pow(exp(-6.9078/(ltmax-ltmin)),(dmax-ltmin));
      }

/*
*     If the calculated cover proportion is less than zero,
*     htot is set to a high value to deflect the search away
*     from such unreal values.
*/
/* Line_170: */
      if ((maxjb >= 1) && (q < 0.)) htot += 1.0e+6;

/*
*     If visual data are supplied, and there is vegetation close
*     to the ground only, and the animals are within that
*     vegetation, the cover proportion correction applies only
*     to the proportion of the direct-line distance potentially
*     obscured by that ground vegetation.
*/
/* Line_190: */

        double vegmax = pow(1.0 - q, dmax);
        double vismax = vegmax*topmax;
      ddsm = dmax*dmax;
      prmax = pa*vismax/ddsm;
    }

    else { // The case kdt = 1 follows:
/* Line_230: */

/*
*     To prevent overflow during computations, an upper limit of
*     70 and a lower limit of -68 are set to QDMAX.
*  JAJC 19/8/2017: Was [-68, 70] but now [-70, 68] to align with other code
*/
        double qdmax = fmin(68, fmax(-70, q * dmax));
        double ssmax = exp(qdmax);
      ddsm = dmax*dmax;
      double pam = pa/ddsm;
      prmax = pam/ssmax;
    }

/*
*     L10, the number of class intervals used for comparisons, is
*     defined.
*/
/* Line_260: */
      l10 = nclass;
      tote = 0.0;

/*
*     TR is the highest r value in the current frequency class;
*     because the first class computed is that furthest
*     from the observer or from the transect line - TR is initially
*     set at CLINT*NCLASS + STT, i.e. equal to or just below the
*     maximum horizontal recognition distance.
*/
      tr = clint*nclass + stt;

/*
*     To compute 'expected' values within each class, each class
*     is subdivided into subclasses, each of width 800/L10, so that
*     CINT is small and L10 x L20=800.  Each subclass
*     corresponds to an observing arc of width CINT ('delta-r')
*     that sweeps forwards ahead of the observer.  L20 is the
*     number of classes in the inner loop.
*     0.5 is added to remove errors due to 'chopping'.
   JAJC 21/8/17: but quotient and lvalue are int, so 0.5 is ignored
*/
      l20 = (800/l10) + 0.5;

/*
*     CINT ('delta-r') is the width of each subclass; it is set
*     at the class width divided by the number of subclasses used.
*/

      cint = clint/(l20);
      hcint = cint/2.0;

/*
*     A function ERMAX is defined to be an exact multiple of CLINT
*/

      iermax = f[3]/clint;
      ermax = iermax * clint;

/*
*     In the event that PRMAX > 1., HTOT is set to a high value.
*/
      if (prmax >= 1.0) htot += prmax*1.0e6;

/*
*     QR is the proportion of the number of animals originally in
*     a subclass still present after animals detected in arcs
*     further out have been disregarded.  It is initially unity.
*/
/* Line_290: */

      qr = 1.0;

/*
*     EXPD is the cumulative expected number within a class; it is
*     set at zero before computations begin; so is EXPDV.
*/
      expd = 0.0;
      expdv = 0.0;
      *s = 0.0;

/*
*     The program computes a correction, CORRN=(NUMA+NUMO)/NUMA, to
*     cater for data sets that contain overtakes.
*/
	  corrn = ((numa+numo)*1.0)/(numa*1.0);

/*
*     Loop 1180, which calculates expected values and compares them
*     with observed values, now begins .....
*/
/* Loop_1180: */
      for (int md = l10 - 1; md >= 0; --md) {	// md = (nclass-1), (nclass-2),..., 0

/*
*     MD is the hth class in the series, calculated from L10
*     downwards, because computations begin around RMAX and continue
*     at progressively decreasing r values.
*/

/*
*     OBSD is the observed value for the hth arc, as supplied
*     in the input to the program.
*     The observed frequency value in the hth class is corrected to
*     to allow for animals that overtake the observer from behind.
*/
	obsd = corrn*val[md];

/*
*     Where the data are radial distance or fixed point data, EXPD
*     accumulates the expected frequency values within a class; it
*     is set at zero before computations for the class begin.
*
*     Where the data are perpendicular distance data (IRY>0), TEXPD
*     accumulates the expected frequency values within each class;
*     it is set at zero before computations for the class begin,
*     while EXPD is allowed to accumulate.
*/

	if (iry <= 0.0) {
	  expd = 0.0;
	} else {
	  texpd = 0.0;
	}

	if (tote < 0.0) tote = 0.0;

/*
*     WL is the lower boundary of the current class.
*/
	wl = tr-clint;

/*
*     If the entire class lies above either the highest calculated
*     or highest preset recognition distances, the expected
*     value (EXPDV) is set at zero and the calculations in Loop
*     870 are bypassed.
*/

	if ((wl > rmax) || (wl > ermax)) {
	  expdv = 0.0;
        }
        else {

/*
*     Two alternative ways of calculating the probability of
*     detection are provided.  If the method is based on the area
*     under the probability curve (!useMedianValues), the area
*     under the curve from y=DH to y=(infinity) is calculated and
*     expressed as YYH.  If IMV=1, this calculation is bypassed.
*/
/* Line_360: */
          double yyh = 0;

          if (!useMedianValues) {
              dh = sqrt(thh * thh + tr * tr);
              yyh = detection_probability(q, dh, pa);
          }
/*
*
*     Loop 870, which calculates the expected number in each
*     arc within the class, and adds them together to produce
*     a progressive total (EXPD and TEXPD), now begins .....
*/

/* Loop_870: */
	dnr = tr-hcint;		// DNR is the central r value in an observing arc.
	dnrh = tr;		// DNRH is the highest r value in an observing arc.

	for (jl=1; jl <= l20; jl++) {    //  DO jl=1,l20

/*
*     DH is the direct-line distance to the outer edge of the
*     observing arc.
*/
	    dh = sqrt(thh*thh+dnrh*dnrh);

/*
*     DNRL is the lowest r value in each arc; it is calculated
*     simply by subtracting the arc width from the highest r
*     value in the arc.
*/
	    dnrl = dnrh-cint;

/*
*     If, in the last class evaluated, DNRL comes to a negative
*     value, it is then set at zero.
*/
	    if (dnrl < 0) dnrl=0.0;

/*
*     DL is the direct-line distance from the observer to the
*     inner edge of the current observing arc.
*/
/* Line_470: */

	    double dl = sqrt(thh*thh+dnrl*dnrl);
            assert(dl <= dh);

/*
*     DINT is the difference between DL and DH, and thus the width
*     of the area under the probability density curve between DH
*     and DL.
*/
	    dint = dh-dl;

/*
*     If useMedianValues is clear, calculation of the probability
*     of detection now follows, using the same approximation to the
*     exponential integral as previously.  Otherwise computation
*     continues from Line_570.
*/

            if (!useMedianValues) {
                  double yyl = detection_probability(q, dl, pa);
                  assert(yyl >= yyh);

/*
*     AUC is the area under the detectability curve between
*     d=DL and d=DH.
*/
                  double auc = yyl - yyh;

/*
*     The mean corrected probability [PR=P(r)] of detecting an
*     individual animal present between d=DL and d=DH is given by
*     dividing the area under the curve by the arc width, then
*     subtracting Pr(rmax).  It is the mean height of the curve.
*/
			pr = auc/dint;

/*
*     YYH is advanced to YYH for the next iteration.
*/
                yyh = yyl;
                
            }       /* End case !useMedianValues */

            else {  /* Start case useMedianValues */

/*
*     Where IMV was set at 1, computation of the probability
*     is based on its median value in the arc rather than
*     the mean height of its density curve.
*
*     DD is the direct-line distance to the centre of the arc.
*/
/* Line_570: */

                dds = thh*thh + dnr*dnr;
                dd = sqrt(dds);

                if (kdt == 1) {
                    double qdd = fmin(68, fmax(-70, q * dd));
                    double pad = pa/dds;
                    pr = pad/exp(qdd);
                }
                
                else {
/*
*      Visibility will be affected not only by lateral vegetation cover but
*      by topographical features too in habitats where the ground is not
*      level. The overall probability of visibility at any distance DD
*      greater than DTMIN (LTMIN) will then be the product (VISDD) of the
*      probability VEGDD that it is unobscured by vegetation and the
*      probability TOPDD that it is also unobscured by topography.  Thus
*      VISDD = VEGDD*TOPDD.  VEGDD is a function of d (VEGDD=(1-Q)**DD)
*      and is calculated at Line_630.
*      A first step in calculating TOPDD is to compute the topographical
*      cover value TCOV from the minimum obscuring distance LTMIN
*      supplied to the program.  This is done after Loop_210 in
*      Subroutine calculate_density.
*      The second step is to calculate TOPDD, provided that DD lies
*      between LTMIN and RMAX.  (If it lies outside those limits,
*      topography will have no effect on visibility, so TCOV is then set
*      equal to 1.) Otherwise TOPDD = EXP(TCOV*(DD-LTMIN)).
*/
                    double topdd;
                    if ((ltmin == 999) || (dd < ltmin) || (ltmax < ltmin)) {
                        topdd = 1.0;
                    } else {
                        topdd = exp(tcov*(dd-ltmin));
                    }

/* Line_630: */
			    double vegdd = pow((1.0-q), dd);
/* Line_640: */
			    double visdd = vegdd*topdd;
			    pr = pa*visdd/dds;
			}

            }  /* End case useMedianValues */

/*
*     PRC [g(r)=P(r)-P(rmax)] is the height of the detectability curve at DD,
*     corrected by subtracting PRMAX.
*/
Line_680:
	    prc = pr-prmax;

/*
*     Because 0<PRC<1, an upper limit of 1 and a lower limit of
*     0 are set to the probability function.
*/
	    if (prc > 1.0) {
		prc = 1.0;
            } else if (prc < 0.0) {
		prc = 0.0;
	    }

/*
*     PRR is the product of PRC=g(r) and the arc of width CINT=∆r.
*/
	    prr = prc*cint;

/*
*     If the median r value in an arc is above both the calculated
*     and preset highest r values, PRR is set at zero and QR at 1.
*     **  JMB: Note that the code does NOT set QR=1 until later **
*
*     If the median r value is below either the calculated
*     or observed r value, PRR is allowed to have a value
*     other than zero.
*/
		
Line_710:
	    if ((dnr >= ermax) && (dnr >= rmax)) prr = 0.0;

/*
*     E [= ∆r.g(r).Q(r)] is calculated by multiplying
*     together the probability of detection in the arc [PRC=g(r)]
*     and the probability [Q(r)] that an individual is still
*     available for detection in that arc. TOTE [=∑(∆r.g(r).Q(r))] is
*     the accumulating probability density function.
*/
Line_720:
	    e = prr*qr;
	    tote += e;

/*
*     The subroutine now calculates the number of detections (ED)
*     expected for radial, fixed-point and perpendicular distance data.
*     For radial data, computation is based on the assumption
*     that the arc of radius r and width ∆r sweeps over
*     a series of plots of unit area.  The total number of detections
*     expected in such a plot will be [plot area] x [apparent density]
*     x [probability of detection].
*
*     With perpendicular distance data, ED is the expected number
*     detected at a distance r from the observer in the strips
*     CINT units wide at a perpendicular distance y=r from the
*     transect line.  Because there are one or two such strips, each
*     L times the area of the individual plot at distance r,
*     the expected no. = ∆r x [DxNsxLJ=d2l] x ∆r x g(r) x Q(r).
*     NS has already been converted to a floating-point number (SNS).
*/
	    if ((ifx == 0) && (iry >= 1)) {
		ed = d2l*cint*e;

/*
*     With radial distance data, each observing arc sweeps out an
*     area L units long and 2r units wide, within which there are
*     2r/CINT x L individual plots at radial distance r.  Thus the
*     expected total number detected  =  [number expected in one
*     plot]x[number of plots]=[DxNsxLJ=d2l] x r x ∆r x g(r) x Q(r).
*/
	    } else if ((ifx == 0) && (iry < 1)) {
		ed = d2l*e*dnrh;

/*
*     With fixed-point data, the expected total number detected
*     = [2ut(=D2LxPd)] x Ps x r x ∆r x g(r) x Q(r).
*/
	    } else {
		ed = d2l*e*dnrh;
	    }

/*
*     The total number expected in a class (EXPD)
*     is obtained by progressively adding together
*     the ED values from each arc with radii between DNRH and DNRL.
*/
/* Line_790: */
		
	    expd += ed;

/*
*     In the case of perpendicular distance data, all the EXPD
*     values within a class are added together to give a TEXPD
*     value for the class (while, for perpendicular distance data, EXPD
*     continues to accumulate from class to class).
*/
	    if (iry >= 1) texpd += expd;

/*
*     In preparation for the next arc in the range, DNRL is
*     now renamed DNRH (because the lowest value in one arc
*     becomes the highest in the next, going inward), while DL
*     similarly becomes DH.
*/
/* Line_810: */

	    dnrh = dnrl;
	    dh = dl;

/*
*     The DNR value for the next arc is DNR minus the
*     arc width (CINT).
*/
/* Line_830: */

	    dnr -= cint;

/*
*     Where PRR is negative, QR is kept at 1, because no animals
*     are in fact removed from the sample at such distances from
*     the observer.
*/
	    if (prr < 0) qr = 1.0;

/*
*     The proportion still present in the next arc will be
*     the proportion present in the present arc less the
*     proportion expected to have been detected already.
   JAJC 23/8/17: this should only be if prr≥0
*/
	    qr = (1.0-prr)*qr;

/*
*     If QR has fallen to a value of less than 0.001 - one
*     animal in a thousand - the program is set to print out the
*     approximate value of rmin where this happens.
*/
	    if ((qr <= 0.001) && (!rlow_found)) {

/*
*     RLOW is the inner boundary of the arc by which 99.9% of the
*     detections are expected to have been made.
*/
		rlow = dnrl-cint;
		rlow_found = true;
	    }

/*
*     Loop 870 now ends.
*/
	}                                          //  END DO Loop_870

/*
*     Loop 1180 is completed by determining the expected value
*     (EXPDV) in the class; this is arrived at differently for
*     radial or fixed-point and perpendicular distance data. In
*     both cases, negative EXPD values (EXPDN) are then added in.
*/
	expdv = (iry > 0) ? texpd : expd;

        } /* End case wl within rmax & ermax */

/*
*     Calculated values are printed out once the program has
*     converged on a minimum, and KPRINT has been set at 1.
*/
		  
/* Line_900: */
        if ((wl <= rmax) && ((kdt <= 1) || (wl <= kdt))) {

/*
*     If D2L has a trial value of zero, calculation of S at this stage
*     is bypassed and a record of this non-computation retained
*     as the temporary variable MSFAIL.  This is done only when the
*     program has converged on a minimum (MTEST=1) and in the final
*     pass through Loop 1180.
*/
	if (d2l > 0) {
	    *s += (expdv/d2l);
	} else if ((mtest == 1) && (md == 0)) {
	    (*msfail)++;
	}

/*
*     Final values of key internal functions are now output if
*     KPRINT=1; otherwise this step is bypassed.
*/
/* Line_930: */

      if (kprint > 0) {

/*
*     For output purposes only, EXPDV is redefined as CALCV.  Negative
*     values of EXPDV are printed as '0.0' in the output because the
*     'observations' are unfloat.
*/
                double midp = tr - (clint/2.0);
                double calcv = fmax(expdv, 0);
                results->midpoints[md] = midp;
                results->calcn[md] = calcv;
                results->obsdn[md] = obsd;
                if (ishow > 0) {
                    char variable = (iry > 0) ? 'y' : 'r';
                    char* method = useMedianValues ? "median" : "mean  ";
                    fprintf(output_results, "  %c=%8.1f     Calc.N(%c)=%9.2f [by %s]    Obsd.N(%c)=%9.1f\n",
                            variable, midp, variable, calcv, method, variable, obsd);
                }

/*
*      If the control variable ISHOW has been set at 1, a variety
*      of intermediate variables is output.
*/
		  
/* Line_1050: */
	if (ishow == 1) {
	    fprintf(output_results, " %12.6f   %12.6f   %12.6f   %12.6f   %12.6f\n\n", prc,qr,tote,expd,expdv);
	}
    }
    }

/*
*     The difference between each observed (OBSD) and expected
*     value (EXPDV) for the class is now calculated.
*
*/
	dif = obsd-expdv;

/*
*      If an upper limit of KDT (>1) has been set for the input
*      data range, the computed difference between observed and
*      computed values is reset at zero.
*/
	if ((kdt > 1) && (wl > kdt)) dif = 0.0;

/*
*     Each difference between the expected (EXPDV) and observed
*     (OBSD) value is squared and a running sum of squares total
*     (TOT) computed.  The final overall sum (FUNC) - which
*     includes HTOT values where parameter values are inappropriate
*     - is then calculated at Step 1300 and returned to the
*     main program.
*/
/* Line_1090: */

	difsq = dif*dif;
	tot += difsq;

/*
*     TR is reduced by CLINT before the subroutine goes to the
*     next class inward.
*/
	tr -= clint;

/*
*     Loop 1180 now ends.
*/
      }							//  END DO Loop_1180

/*
*     The next two outputs occur only if the range of distance
*     values includes the 99.9% r value and/or an r=0 value.
*
*     If the search for a minimum is complete, and KPRINT=1,
*     the program now prints out a value for RLOW and a file
*     which tabulates observed and calculated frequencies.
*/

      if (kprint > 0) {

	if ((rlow_found) && (ishow > 0)) {
	  fprintf(output_results, "\n 99.9%% r value (rmin) =%7.2f m\n", rlow);
	}

/*
*     If the range of distance values includes a g(0) value, the
*     program estimates detectability at g(y)=0.
*/
		  
	if (tr <= 0.1) {
	  fprintf(output_results, " Est.Detectability at g(y=0): %4.2f\n", tote);
	}

      }  /* end (kprint > 0) */

/*
*     The totalled sum of squares of the differences between
*     observed and calculated values is now defined as FUNC.
*/
	
Line_1300:
      *func = tot + htot;

      if (kprint == 1) {
	fprintf(output_results, " Final Difference at Minimum = %15.6f\n\n", *func);
      }
}
//      END SUBROUTINE givef

/*******************************************************************************
*
*	QSF		Quadratic Surface Fit
*
*
*******************************************************************************/


/*******************************************************************
*
*      QSF
*
*      This subroutine enables a surface fitting procedure to be used
*      when the number of model iterations (maxjb) is set at 1.  It
*      estimates the variances of the parameters by fitting a quadratic
*      surface in the region of the minimum (for estimation of the
*      variance-covariance matrix).
*
*/

void qsf (double f[NUM_SHAPE_PARAMS],
	  double func,
	  double approx,
	  double dist,
	  double durn,
	  double rate,
	  double step[NUM_SHAPE_PARAMS],
	  double stopc,
	  double *s,
	  double val[80],
	  double clint,
	  double stt,
	  double tcov,
	  double thh,
	  int ishow,
	  int *kprint,
	  double ltmax,
	  double ltmin,
	  int nclass,
	  int numa,
	  int numo,
	  int *msfail,
	  int mtest,
	  double *estden,
	  double *sden,
	  double *hstst,
	  double pd,
	  double ps,
	  int ifx,
	  int iprint,
	  int iry,
	  int kdt,
	  int km,
	  int nap,
	  /* int neval, */
	  int nop,
	  int ns,
	  /* int nvals, */
	  int np1,
	  double g[21][20],
	  double h[21],
	  int maxjb,
	  /* int jprint, */
	  calc_results *results)
{

/*
*      Double precision is set for most real numbers, together with
*      common values, dimensions and symbols for key variables.
*
*/
       int i, i1, i2, i3, i4, i5, j, j0, j1, k, krun;
       int nless1, in, jk, jless1, l;
       int iplus1, klessi, nu, nl, ij, ndf, jplus1;
       int mnpd = 0, iless1, neval;
       double pstar[20], aval[20], bmat[210], ao/*, dmax*/;
       double temp, ymin, vc[210], var[NUM_SHAPE_PARAMS], vra, sa, sb, sd2l;
       double sd, vrb, vrd2l, den, simp;
       double pmin[20], /* pbar[20], */ pstst[20], t;
       double df, vrs;
       double test, cl1, cl2;


      fprintf(output_results, "\n\n Fitting of quadratic surface in region of minimum\n\n\n");

/*
*      The fitting of the quadratic surface follows the procedure
*      outlined by Nelder and Mead exactly and, where possible,
*      the notation in the comments corresponds to theirs also.
*/

      neval = 0;
      krun = 0;

/*
*      Further function evaluations are counted in NEVAL.
*
*      The final simplex is expanded to overcome rounding errors.
*
*      KRUN is a control variable which allows an automatic rerun
*      through the program if the criterion TEST is negative during
*      quadratic surface fitting.  The program begins again using
*      the computed F[I] values and values of STEP set at one
*      tenth the F[I] values.
*
*      SIMP is set to a higher value than the stopping criterion STOPC.
*/

      simp = 2*stopc;

/* JMB: The loop at Outer0 runs only once, because the routine either returns or it breaks out of the loop. */

Outer0:
      for (i=0; i < np1; i++) {				//  DO i=1,np1
Line_20:
	test = fabs(h[i]) - func;

	if (test <= approx) {	/* Rerun required */
	  if (krun <= 0) {
	    krun = 1;

	    fprintf(output_results, " Rerun needed with new initial & step values \n");
	    for (in=0; in < nop; in++) {	//  DO in=1,nop
	      step[in] = f[in]/10.0;
	    }					//  END DO
	  }
	  else {

	    fprintf(output_results, " Rerun data with class intervals altered \n");
	  }
	  return;
	}

	if (test < simp) {

	  for (j=0; j < nop; j++) {		//  DO j=1,nop
	    pstst[j] = 1000*(g[i][j]-f[j]) + g[i][j];
	  }					//  END DO

	  givef (pstst, &h[i], s, val, clint, /* pd, */ stt, tcov,
		thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, *kprint, /* dmax, */
		ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
		mtest, true, results);

	  neval++;
	  goto Line_20;
	}

        break;
      }							//  END DO outer0

Line_60:
      ao = h[0];

/*
*
*      The function values Y0(I) are calculated and stored in AVAL.
*/
      for (i=0; i < nap; i++) {				//  DO i=1,nap
	i1 = i+1;
	for (j=0; j < nop; j++) {	//  DO j=1,nop
             pstar[j] = (g[0][j]+g[i1][j])/2.0;
	}				//  END DO

	givef (pstar, &aval[i], s, val, clint, /* pd, */ stt, tcov,
	       thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, *kprint, /* dmax, */
	       ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
	       mtest, true, results);

	neval++;
      }							//  END DO

/*
*      The matrix B(I,J) is calculated, and the lower diagonal
*      section stored in the vector BMAT.
*
* Comment by J.Begg: the C code below is based on the FORTRAN code
* but in order to get the array bounds to start at 0 some of the
* arithmetic for calculating the array indices has been subtly
* altered, especially for the BMAT array.
*/

Outer:
      for (i=0; i < nap; i++) {				//  DO i=1,nap

          i1 = i-1;
          i2 = i+1;

		// JMB: 'inner1' block is skipped the first time around 'outer' loop
		//       due to 'j' loop termination test
Inner1:
	  for (j=0; j <= i1; j++) {		//  DO j=1,i1
	    j1 = j+1;

	    for (k=0; k < nop; k++) {	//  DO k=1,nop
	      pstst[k] = (g[i2][k]+g[j1][k])/2.0;
	    }				//  END DO inner2

	    givef(pstst, hstst, s, val, clint, /* pd, */ stt, tcov,
		  thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, *kprint, /* dmax, */
		  ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
		  mtest, true, results);

	    neval++;
	    l = i*i2/2 + j;
	    bmat[l] = 2.0 * (*hstst+ao-aval[i]-aval[j]);
	  }					//  END DO inner1

      }						 //  END DO outer

      l = 0;
      for (i=0; i < nap; i++) {		//  DO i=1,nap
	i1 = i+1;
	l = l+i1;
	bmat[l-1] = 2.0 * (h[i1]+ao-2.0*aval[i]);
      }					//  END DO


/*
*
*      The vector A(I) is calculated and stored in AVAL.
*/
      for (i=0; i < nap; i++) {		//  DO i=1,nap
         i1 = i+1;
         aval[i] = 2.0*aval[i]-(h[i1]+3.0*ao)/2.0;
      }					//  END DO

/*
*      The matrix Q is calculated and stored in the matrix G.
*      Considering the usual orientation of rows and columns,
*      TRANS(Q) is stored in G.
*/
	
      for (i=0; i < nop; i++) {		//  DO i=1,nop
	pmin[i] = g[0][i];	// JMB: g(1,i) changed to g[0][i]
      }					//  END DO

Outer2:
      for (i=0; i < nap; i++) {			//  DO i=1,nap
	i1 = i+1;
	for (j=0; j < nop; j++) {	//  DO j=1,nop
	  g[i1][j] -= g[0][j];	// JMB: g(1,j) changed to g[0][j]
	}				//  END DO
      }						//  END DO Outer2

Outer3:
      for (i=0; i < nap; i++) {			//  DO i=1,nap
	i1 = i+1;
	for (j=0; j < nop; j++) {	//  DO j=1,nop
           g[i][j] = g[i1][j];
	}				//  END DO
      }						//  END DO Outer3

/*
*
*      The matrix B is inverted, using the modified square root
*      method (see SAZONOV, Geodeziya i Aerofotosyemka, No.6, 1962).
*
* Comment by J.Begg: this is another case where the converted FORTRAN
* code required tweaking to produce the array indices in the correct
* order.
*/
      np1 = nap;
      nless1 = nap-1;
      i1 = 1;
Outer4:
      for (i=1; i<np1; i++) {				// DO i=2,np1
	iless1 = i-1;
	i3 = 0;
Inner41:
	for (j=1; j<=iless1; j++) {		// DO j=2,iless1
	  i2 = i1+j;
	  jless1 = j-1;
Inner42:
	  for (j0=0; j0<=jless1; j0++) { // DO j0=1,jless1
	    i4 = i3+j0+1;
	    i5 = i1+j0;
	    bmat[i2] -= bmat[i4]*bmat[i5];
	  }				 // END DO Inner42
	  i3 = i3+j+1;
	}					// END DO  Innter41

	i3 = -1;
	i5 = i1+i;
Inner43:
	for (k=0; k<=iless1; k++) {		// DO k=1,iless1
	  i2 = i1+k;
	  i3 = i3+k+1;
	  temp = bmat[i2]/bmat[i3];
	  bmat[i5] -= bmat[i2]*temp;
	  bmat[i2] = temp;
	}					// END DO Inner43

/*
*        If the matrix B is not positive definite, MNPD is set =-1,
*        the estimation of variances ends, and the parameter
*        estimates are printed out.
*/
		  
	if (bmat[i5] <= 0) {

	  fprintf(output_results, " Matrix to be inverted not positive definite \n");
	  mnpd = -1;
	  goto Line_205;
	}

	i1 = i1+i+1;
      }							// END DO Outer4

      i1=1;
Outer5:
      for (i=1; i < np1; i++) {				//  DO i=2,np1
	iless1 = i-1;
Inner51:
	for (j=0; j <= iless1; j++) {		//  DO j=1,iless1
	  i2 = i1+j;
	  jplus1 = j+1;
Inner52:
	  for (k=jplus1; k <= iless1; k++) {	//  DO k=jplus1,iless1
	    i3 = j+k*(k+1)/2;
	    i4 = i1+k;
	    bmat[i2] += bmat[i3]*bmat[i4];
	  }				//  END DO inner52
	  bmat[i2] = -bmat[i2];
	}					//  END DO inner51
	i1++;
      }							//  END DO outer5

/*
*  Comment by J.Begg: This time we use the FORTRAN loop bounds and just
*  get the array index by subtracting one when referencing the array.
*/
      i1=0;
      for (i=1; i <= np1; i++) {	//  DO i=1,np1
	i1 += i;
	bmat[i1-1] = 1.0/bmat[i1-1];
      }					//  END DO

Outer6:
      for (i=0; i < nless1; i++) {			//  DO i=1,nless1
	iplus1 = i+1;
Inner61:
	for (k=iplus1; k < np1; k++) {		//  DO k=iplus1,np1
	  i1 = k*(k+1)/2;
	  i2 = i1+k;
	  i3 = i1+i;
	  temp = bmat[i2]*bmat[i3];
	  klessi = k-i;
Inner62:
	  for (j=1; j <= klessi; j++) {	//  DO j=1,klessi
	    j0 = i+j-1;
	    i4 = j0*(j0+1)/2+i;
	    i5 = i1+j0;
	    bmat[i4] += bmat[i5]*temp;
	  }				//  END DO inner62
	  bmat[i3] = temp;
	}					//  END DO inner61
      }							//  END DO outer6

/*
*      (B**-1)*A is calculated , and stored in H
*/

Outer7:
      for (i=0; i < nap; i++) {			//  DO i=1,nap
	h[i] = 0.0;
Inner71:
	for (j=0; j < nap; j++) {	//  DO j=1,nap
	  if (j <= i) {
	    ij = i*(i+1)/2+j;
	  } else {
	    ij = j*(j+1)/2+i;
	  }
	  h[i] += bmat[ij]*aval[j];
	}				//  END DO inner71
      }						//  END DO outer7

/*
*      The estimated minimum value (YMIN) and its position (PMIN)
*      are calculated and printed out if IPRINT=1.
*/
      ymin = 0;

      for (i=0; i < nap; i++) {	//  DO i=1,nap
	ymin += h[i]*aval[i];
      }				//  END DO

      ymin = ao-ymin;

Outer8:
      for (i=0; i < nop; i++) {			//  DO i=1,nop
	pstst[i] = 0.0;
	for (j=0; j < nap; j++) {	//  DO j=1,nap
	  pstst[i] += h[j]*g[j][i];
	}				//  END DO
      }						//  END DO outer8

      for (i=0; i < nop; i++) {	//  DO i=1,nop
	pmin[i] -= pstst[i];
      }				//  END DO


#if (NUM_SHAPE_PARAMS != 4)
#error "Check output statements in QSF cover all pmin[] and f[] values!"
#endif
      if (iprint == 1) {

	fprintf(output_results, " Minimum of fitted quadratic surface is %13.6e at: \n", ymin);
	fprintf(output_results, " %13.6e %13.6e %13.6e %13.6e\n", pmin[0], pmin[1], pmin[2], pmin[3]);

	fprintf(output_results, "\nCompare with minimum found by iteration %13.6e at: \n", func);
	fprintf(output_results, " %13.6e %13.6e %13.6e %13.6e\n", f[0], f[1], f[2], f[3]);

	fprintf(output_results, "\n\n If difference is large, information matrix and errors are inaccurate\n\n\n");
      }

/*
*      Q*(B**-1)*TRANSQ is calculated, and its lower diagonal
*      section stored in the vector VC.
*/
Outer9:
      for (i=0; i < nop; i++) {					//  DO i=1,nop
Inner91:
	for (j=0; j < nap; j++) {			//  DO j=1,nap
	  h[j] = 0.0;
Inner92:
	  for (k=0; k < nap; k++) {		//  DO k=1,nap
	    if (k <= j) {
	      jk = j*(j+1)/2+k;
	    } else {
	      jk = k*(k+1)/2+j;
	    }
	    h[j] += bmat[jk]*g[k][i];
	  }					//  END DO inner92
	}						//  END DO inner91

Inner94:
	for (j=i; j < nop; j++) {			//  DO j=i,nop
	  ij = i*(i+1)/2+j;
	  vc[ij] = 0.0;
	  for (k=0; k < nap; k++) {	//  DO k=1,nap
	    vc[ij] += h[k]*g[k][j];
	  }				//  END DO
	}						//  END DO inner94
      }								//  END DO outer9

/*
*      The diagonal elements of VC are stored in VAR.
*/
      for (i=0; i < nop; i++) {		//  DO i=1,nop
	j = i*(i+1)/2+i;
	var[i] = vc[j];
	if (step[i] < approx) {
	  var[i]=0.0;
	}
      }					//  END DO

/*
*      The inverse of the information matrix is printed out.
*/
      if (iprint == 1) {
	 fprintf(output_results, " Inverse of information matrix \n\n");

         for (i=1; i <= nop; i++) {    /* JMB: Was:  DO i=1,nop
				       *       Note these loop bounds must be maintained otherwise the
				       *       calculation of nu & nl goes haywire.
				       */
           nu = i*(i+1)/2;
           nl = i*(i-1)/2+1;

	   for (j=nl-1; j < nu; j++) fprintf(output_results, " %13.6e", vc[j]);	// JMB: now we decrement the loop bounds
	   fprintf(output_results,"\n");
         }  //  END DO

       }

       fprintf(output_results, "\n%2i additional evaluations have been used\n\n", neval);

       if (iprint == 1) {

	 fprintf(output_results, " End  of  quadratic  surface  fitting \n");
       }

/*
*      The estimated 'best fit' values of the parameters, together
*      with their estimated standard errors, are now computed,
*      tabulated and printed out.  Also, KPRINT is set at 1.
*/
	
Line_205:
      *kprint = 1;

/*
*      A density estimate (DEN) is calculated from D2L=F[2] by
*      correcting units to no./ha, and dividing by 2L in the case of
*      line transect data, and by 2Vt in the case of fixed-point data.
*/
      if ((km > 0) && (ifx == 0)) {
	den = (1.0e6*f[2]) / (ns*dist*pd);
      } else if ((km == 0) && (ifx == 0)) {
	den = (1.0e4*f[2]) / (ns*dist*pd);
      } else {
	den = (1.0e4*f[2]) / (2.0*rate*durn*ps*pd);
      }

/*
*      If the variance/covariance matrix is not positive definite
*      (MNPD=-1), the number of degrees of freedom (NDF) is set =0
*/
      if (mnpd < 0) ndf = 0;

/*
*      The next series of steps determines the number of degrees of
*      freedom used to compute the variance estimate.  NDF is set
*      at one less than the number of values (NVALS) supplied to
*      the program, less the number of parameters estimated by
*      the program itself.
*/
	
Line_220:
      ndf = nclass-nap-1;
      df = ndf;

/*
*      If NDF is zero, VRS is set to an arbitrarily high value.
*/
      if (df <= approx) df = approx;
      vrs = *hstst/df;

/*
*
*      The standard errors are computed and the results printed out.
*/
      vra = fabs(var[0]*2.0*vrs);
      sa = sqrt(vra);

      if (var[1] > 1.0e+12) var[0] = 1.0e+12;

      if (var[3] > 1.0e+12) var[3] = 1.0e+12;

      vrb = fabs(var[1]*2.0*vrs);
      sb = sqrt(vrb);

      vrd2l = fabs(var[2]*2.0*vrs);
      sd2l = sqrt(vrd2l);

      if ((km == 1) && (ifx == 0)) {
	sd = (1.0e6*sd2l) / (ns*dist*pd);
      } else if ((km == 0) && (ifx == 0)) {
	sd = (1.0e4*sd2l) / (ns*dist*pd);
      } else {
	sd = (1.0e4*sd2l) / (2.0*rate*durn*ps*pd);
      }

/*
*
*     The program now determines confidence limits for DEN.
*     The p=0.05(2) distribution of t with sample number is first
*     approximated by a function T.
*/
	
      t = (1/(0.093914*pow(nclass,2.09372)))+2.03046;

/*
*     The 95% lower (CL1) and upper (CL2) confidence limits are now
*     calculated from the mean and standard deviation.
*/
      cl1 = den - t*sd;
      cl2 = den + t*sd;


	fprintf(output_results, "\n\n Estimated parameter values: \n\n\n");

      fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results, " x                  x                  x                  x                  x\n");
      fprintf(output_results, " x    Parameter     x      Value       x  Standard Error  x       Unit       x\n");
      fprintf(output_results, " x                  x                  x                  x                  x\n");
      fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results, " x                  x                  x                  x                  x\n");
      fprintf(output_results, " x ESTIMATED        x                  x                  x                  x\n");

      if (km == 0) {
        if (ndf == 0) {

	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x (indeterminate)  x  indivs/hectare  x\n", den);
        } else {

	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x    %10.3f    x  indivs/hectare  x\n", den, sd);
        }

      } else if (km == 1) {
        if (ndf > 0) {

	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x    %10.3f    x   indivs./sq.km. x\n", den, sd);
        } else {

	  fprintf(output_results, " x DENSITY  (D)     x    %10.2f    x    %10.3f    x (indeterminate)  x\n", den, sd);
        }
      }

      fprintf(output_results, " x                  x                  x                  x                  x\n");
      fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results, " x                  x                  x                  x                  x\n");
      fprintf(output_results, " x Conspicuousness  x                  x                  x                  x\n");


      if (ndf == 0) {

	fprintf(output_results, " x Coefficient  (a) x    %10.4f    x (indeterminate)  x      metres      x\n", f[0]);
      } else {

	fprintf(output_results, " x Coefficient  (a) x    %10.3f    x    %10.3f    x      metres      x\n", f[0], sa);
      }

      if (kdt != 1) {

	fprintf(output_results, " x                  x                  x                  x                  x\n");
	fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
	fprintf(output_results, " x                  x                  x                  x                  x\n");
	fprintf(output_results, " x Cover            x                  x                  x                  x\n");

        if (ndf <= 0) {

	  fprintf(output_results, " x Proportion   (c)  x    %10.4f    x (indeterminate)  x                  x\n", f[1]);
        } else {

	  fprintf(output_results, " x Proportion  (c)  x    %10.4f    x    %10.4f    x                  x\n", f[1], sb);
        }

      } else { /* kdt == 1 */

	fprintf(output_results, " x                  x                  x                  x                  x\n");
	fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
	fprintf(output_results, " x                  x                  x                  x                  x\n");
	fprintf(output_results, " x Attenuation      x                  x                  x                  x\n");

        if (ndf <= 0) {

	  fprintf(output_results, " x Coefficient  (b) x    %10.4f    x (indeterminate)  x    per  metre    x\n", f[1]);

        } else {

	  fprintf(output_results, " x Coefficient  (b) x    %10.4f    x    %10.4f    x    per metre     x\n", f[1], sb);
         }

      }

	fprintf(output_results, " x                  x                  x                  x                  x\n");
	fprintf(output_results, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n");

/*
*
*     The next two parameters are output only if at least one parameter
*     is being estimated.
*/
      if (nap >= 1) {

	fprintf(output_results, " 95%% confidence limits for density estimate:%10.3f %10.3f\n", cl1, cl2);
      }

      *estden = den;
      *sden = sd;
}

//      END SUBROUTINE

/******************************************************************************/


//      SUBROUTINE calculate_density (params, header, outfile, graph_file)

void calculate_density (calc_params *params
			,const char *header, const char *outfile, const char *graphfile
			)
{

/*
*     The variables declared below are in rough alphabetical order
*/

      int nsize[MAX_OBSERVATIONS], nbsz[MAX_OBSERVATIONS];
      int bootstrap;
      int i, ia, ie, iflag, ifx, ig, ih;
      int in, iprint, irow, iry, iseed, ishow;
      int j, jprint, jv, kdt, km, kprint, max;
      int maxjb, mfail, msfail, mtest, nap, neval, ngroups;
      int nloop, nop, np1, ns, numa, numest, numo, nclass;
      int numoin, numain, nvals;
      bool iqsf = false;
      double a, approx, b, c, clint, dcoeff, dist, s;
      double estden, func, hmean;
      double hstar, hstd, hstst, durn, ltmin, ltmax;
      double pd, ps, rate, savemn;
      double sns, stopc, stt, test, tcoeff1, tcoeff2;
      double tcoeff3, tcov, tden, thh;
      float estdmax, ttrden, tcl1,tcl2,cl1,cl2;
      float resamp_dist[MAX_OBSERVATIONS], trden[5000];
      double val[80], valt[80];
      double g[21][20], step[NUM_SHAPE_PARAMS], f[NUM_SHAPE_PARAMS], ft[NUM_SHAPE_PARAMS];
      double h[21], pbar[20], pstar[20], pstst[20];
      double coeff1[5000], coeff2[5000], coeff3[5000];
      double den[5000];
      bool dmaxIsPreset = true;

/*
*     The program accepts up to 10000 data values, each being the total
*     number of observations (N[r] or N[y]) within the class
*     intervals, beginning with that nearest r=0 or y=0 (R[0],
*     NSIZE[0] and ANGLE[0], if supplied).  If ANGLE[] is not
*     supplied and IRY==2, computation still proceeds using N[y]
*     values entered in place of N[r] values.
*
*     This program is designed to receive a list-directed data file with
*     a name of the form <filename>.dat and either spaces or commas
*     separating the values.
*
*     The header line should begin with the computer run number each
*     time to avoid confusion, and be bounded by quotation marks.  Other
*     entries must be separated either by commas or spaces (and
*     pre-checked before running).
*
*
*     The program is set to seek a minimum using up to MAX=750
*     iterations, begin with four parameters (NOP) and set a testing
*     point (NLOOP) equal to unity.
*/
      max = 750;
      nop = NUM_SHAPE_PARAMS;
      nloop = 1;

/*
*  Initialise miscellaneous local variables - J.Begg 24-May-2011
*/
      a = 0.0;

      nvals = params->nvals;
      clint = params->clint;
      dist = params->dist;
      thh = params->thh;
      ltmin = params->ltmin;
      ifx = params->ifx;
      iry = params->iry;
      ns = params->ns;
      km = params->km;
      kdt = params->kdt;
      iprint = params->iprint;
      jprint = params->jprint;
      ishow = params->ishow;
      maxjb = params->maxjb;
      durn = params->durn;
      rate = params->rate;
      pd = params->pd;
      ps = params->ps;

/*
*  Check selected input parameters for limits
*/
      if (nvals > MAX_OBSERVATIONS) {
	printf("FATAL ERROR: NVALS value(%d) exceeds allowed max(%d)!\n", nvals, MAX_OBSERVATIONS);
	abort();
      }

      if (maxjb > 5000) {
	printf("FATAL ERROR: MAXJB value(%d) exceeds allowed max(%d)!\n", maxjb, 5000);
	abort();
      }

/*
*  Initialize model parameter values and step sizes
*/
      search_params start = {0};
      select_search_parameters(params, &start);

Loop_10:
      for (ig=0; ig < nop; ig++) {    		//  DO ig=1,nop
        f[ig] = start.f[ig];
        step[ig] = start.step[ig];
      }						//  END DO Loop_10
      clint = start.clint;

/*
*     Care is required to ensure that the group size data is submitted
*     in PRECISELY the same sequence as the corresponding R[IH] data,
*     and that non-overtaking cases where r=0 are altered to r=0.01.
*/
	
Loop_20:
      for (ih=0; ih < nvals; ih++) {		//  DO ih=1,nvals
        nsize[ih] = params->nsize[ih];
      }						//  END DO Loop_20

/*
*  Create the output file
*/
      output_results = fopen(outfile, "w");	// OPEN (UNIT=2,FILE=outfile,STATUS='NEW',IOSTAT=ios,ERR=Line_1930)
      if (output_results == NULL) {

	/* File could not be opened */

	perror("Could not open outfile");
	return;
      }

/*
*     The header line now begins the program output.
*/
	
      fprintf(output_results, "%s\n", header);

/*
*     The program now prints out a statement of the various data
*     supplied as input, as a check on input accuracy.
*/

      fprintf(output_results, "\nINPUT:\n\n");

      if (ifx == 0) {
	switch(iry) {
		case 0:	fprintf(output_results, " Line transect data, based on radial distances from observer\n");
			break;
		case 1:	fprintf(output_results, " Line transect data, based on radial distances and horizontal angles\n");
			break;
		case 2:	fprintf(output_results, " Line transect data, using precalculated perpendicular distances\n");
			break;
		default: fprintf(output_results, " Unkown line transect data (IRY=%d)\n", iry);
	}
	switch (ns) {
		case 1: fprintf(output_results, " Observations from one side of transect line only\n");
			break;
		case 2: fprintf(output_results, " Observations from both sides of transect line\n");
	}
      }
      else {
	fprintf(output_results, " Fixed observing point data\n");
      }

      if (km == 1) {
	fprintf(output_results, " Transect lengths in km\n");
      } else if (km == 2) {
	fprintf(output_results, " Detection distances and transect lengths in km\n");
      }

      if (kdt > 1) {
	fprintf(output_results, " Designated distance range: %5.1f - %5i\n", start.stt, kdt);
      }

/*
*     The program prints out the class interval width (CLINT) and
*     either the total transect length (DIST) or the total time
*     spent (DURN) at fixed points, then the class interval set
*     and the total transect length travelled.
*/

      if (ifx == 1) {
	fprintf(output_results, " Class interval width =%7.1f m.    Total time spent =%7.1f min.\n", clint, durn);
      }
      else {

	if (km == 0) {

	  /* Transect lengths have been expressed in metres. */

	  fprintf(output_results, " Class interval width =%7.1f m.   Total transect length (L) =%10.3f m.\n", clint, dist);
	}

	else {

	  /* Transect lengths have been expressed in kilometres. */

	  fprintf(output_results, " Class interval width =%7.1f m.   Total transect length (L) =%10.3f km.\n", clint, dist);

	  /* Convert from km to metres */

	  dist *= 1000;
        }
        fprintf(output_results," Total time spent =%7.1f min.\n", durn);
      }

      estdmax = f[3];

Line_90:
      fprintf(output_results," Overall population movement rate =%4.1f m/min.\n", rate);

      if (ltmin < 999) {
	fprintf(output_results, " Topography uneven; approximate minimum obscuring distance =%6.1f m.\n", ltmin);
      }
      else {
	fprintf(output_results, " Topography approximately level\n");
      }

      fprintf(output_results, "\nOUTPUT:\n\n");

/*
*     NCLASS is the number of distance classes in the selected range.
*/
    stt = start.stt;
    nclass = start.nclass;
    ngroups = start.ngroups;

/*
*     The original values of NUMA and NUMO are retained (as NUMOIN
*     and NUMAIN) so they can be printed in the output.
*/
    numoin = numo = start.numo;
    numain = numa = start.numa;

/*
 *     The frequency distribution of the original data is now saved,
 *     as VALT[IG].
 */

    for (ig=0; ig < nclass; ig++) {		//  DO ig=1,nclass
        valt[ig] = val[ig] = start.val[ig];
    }

/*
*     For line transect data, the movement-corrected overall distance travelled (LJ) is
*     now calculated, overriding the DIST value submitted originally.
*     DURN and RATE are both set at 0 to avoid later computation problems.
*     
*/
            if (ifx == 0)       {
                dist *= start.estj;
                durn = 0;
                rate = 0;
            }

/*
*     F[2] is now raised in value to approximate D2LJ*Ns*Pd in the case
*     of line transect data, or D2ut*Ps*Pd for fixed point data, in
*     order to prepare for comparisons with observed values within
*     subroutine GIVEF.
*/
      sns = ns;

      if (ifx == 0) {
         f[2] *= (sns*dist*pd)/1.e4;
         step[2] *= (sns*dist*pd)/1.e4;
      } else {
         f[2] *= (2.*ps*pd*rate*durn)/1.e4;
         step[2] *= (2.*ps*pd*rate*durn)/1.e4;
      }

/*
*     The initial values of 'a', 'b', 'D2L', and 'dmax' are retained
*     as FT to make possible reruns of calculations.
*/

Loop_210:
      for (ia=0; ia < nop; ia++) {	//  DO ia=1,nop
        ft[ia] = f[ia];
      }					//  END DO Loop_210

/*
*
*     The program calculates a topographical cover value
*     (TCOV).  If the topography is effectively level (LTMIN=999)
*     or the value of f[3] is less than the LTMIN value supplied,
*     then TCOV is set at zero and topography has no effect on
*     computation of the model.  Otherwise, if LTMIN
*     values are supplied, a TCOV value is calculated.
*     TCOV=ln(0.001)/RMAX-LTMIN) = -6.9078/(RMAX-LTMIN)).
*/
	
      ltmax = f[3];

      if ((ltmin < 999) && (ltmax > ltmin)) {
         tcov = -6.9078/(ltmax - ltmin);
       } else {
         tcov = 0.;
      }

/*
*     The stopping criterion (STOPC) is set at a suitable value,
*     based on the type of data supplied: radial, perpendicular
*     distance data to the limit of visibility, or perpendicular
*     distance data to a distance limit (KDT>1). Tested on data.
*/
	
      if ((iry == 1) || ((iry == 2) && (kdt <= 1))) {
        stopc = 0.0001;
      } else if ((iry == 2) && (kdt > 1)) {
        stopc = 0.00001;
      } else {
        stopc = 0.0001;
      }

/*
*     The program now increases the value of STOPC to allow for
*     highly variable data, by multiplying by NUMAIN/NVALS.
*/
      stopc *= numain;
	  stopc = stopc / nvals;

/*
*     If progress reports are required (IPRINT=1), the program
*     prints a heading for them.
*/
	
      if (iprint == 1) {
	fprintf(output_results,"\n Progress Report (every 1 to 3 function evaluations):\n\n");
	fprintf(output_results,"\n Sequence:  Evaln.no.  Min.difnce  Cnsp.cfnt.  Cover prpn.  Density x attributes   max.dist.\n\n");
      }

/*
*     The term 'APPROX' is used to test closeness to zero.
*/
	
      approx = 1.e-12;

/*
*     If no values have been submitted in the input, the program sets
*     the step variations at: a=1.0, b=0.5, and c= 2.0
*/

      if (fabs(a) < approx) {
        a = 1.0;
        b = 0.5;
        c = 2.0;
      }

/*
*     NAP is the number of parameters to be varied (i.e. with STEP
*     not equal to zero.
*
*/
      nap = 0;
      iflag = 0;
      kprint = 0;
      dcoeff = 0;

/*
*     If all STEP sizes have been set at zero, and MAXJB has not been
*     set at 1, then computation goes to Label 1470, calculates the
*     set of values resulting from the values of a, b, D2L and dmax
*     supplied, and ends.  Otherwise it uses NAP to indicate the number
*     of submitted parameters to be varied during subsequent iterations.
*/

Loop_350:
      for (i=0; i < nop; i++) {		//  DO i=1,nop
        if (fabs(step[i]) > approx) {
          nap++;
        }
      }					//  END DO Loop_350

      if ((nap == 0) && (maxjb != 1)) {
	kprint = 1;
	sns = (ns == 0) ? 2.0 : ns;

	if (km > 0) {
          estden = (1.e6*f[2])/(sns*dist*pd);
	} else if (ifx == 0) {
          estden = (1.e4*f[2])/(sns*dist*pd);
	} else {
          estden = (1.e4*f[2])/(2.*rate*durn*ps*pd);
	}

	goto Line_1470;
      }

/*
*     To enable parameter estimation using bootstrapping, the basic
*     MINIM routine is run a predetermined MAXJB times, beginning
*     with an initial run.  A number of functions are set at zero first.
*/

Line_370:
      mfail = 0;
      msfail = 0;
      tden = 0.0;
      ttrden = 0.0;
      tcoeff1 = 0.0;
      tcoeff2 = 0.0;
      tcoeff3 = 0.0;

/*
*     Seed the random number generator
*/
	
      iseed = 0;
Loop_370:
      for (in=0; in < nvals; in++) {		//  DO in = 1, nvals
	if (iseed > 150000) iseed /= 100;
        iseed += (params->r[in] * nsize[in] * 72493);
      }						//  END DO Loop_370
      srandnag(iseed);

/* ==============================================================
*  ======                                                  ======
*  ======             Start of the main loop               ======
*  ======                                                  ======
*  =============================================================*/

Loop_1410:
    for (bootstrap = 0; bootstrap < maxjb; bootstrap++) {	//  DO bootstrap=1, maxjb
	params->bootstrap = bootstrap;

/*
*     The flag variable MTEST is first set at zero.
*/
		  
	mtest = 0;

          if (bootstrap == 0) {

              if ((maxjb != 1) && (nap <= 0)) goto Line_1320;

          } else {
/*
*     To calculate a bootstrapped distribution, values of R(IN) and
*     the corresponding group NSIZE(IN) are to be chosen at random
*     with replacement, based on a randomly-selected value of IN,
*     which ranges between 1 and NVALS, the total number of groups of
*     animals detected.  Loop 510 selects these values.
*
*/
              resample (start.distance, nsize, nvals, resamp_dist, nbsz);
              freq_distrib(nclass, stt, clint, nvals, kdt, resamp_dist, (iry > 0), nbsz,
                           &numa, &numo, val, &ngroups);
          }
/*
*     The calculated set of values, VAL[IC], in each class interval
*     is now printed out for the set of data concerned if JPRINT=1.
*/

Line_620:
	if ((jprint == 1) && (maxjb != 1)) {

	  fprintf(output_results, "\n\nBootstrap Replicate No. =%4i      Individuals per class:\n\n", bootstrap+1);
	  for (i=0, j=0; i < nclass; i++) {
	    fprintf(output_results, "%6.0f", val[i]);
	    j++;
	    if (j == 10) {
	      fprintf(output_results, "\n");
	      j = 0;
	    }
	  }
        }

/*
*     The initial simplex of program MINIM is now set up.
*
*/

Loop_660:
	for (i=0; i < nop; i++) {	//  DO i=1,nop
	  g[0][i]=f[i];
	}				//  END DO Loop_660

	irow = 1;

Loop_690:
	for (i=0; i < nop; i++) {		//  DO i=1,nop
          if (fabs(step[i]) >= approx) {
	    for (j=0; j < nop; j++) {	//  DO j=1,nop
	      g[irow][j]=f[j];
	    }				// END DO j
	    g[irow][i] += step[i];
	    irow++;
	  }
	}					//  END DO Loop_690

	np1 = nap+1;
	neval = 0;

Loop_730:
	for (i=0; i < np1; i++) {		//  DO i=1,np1

	  for (j=0; j < nop; j++) {	//  DO j=1,nop
	    f[j] = g[i][j];
	  }				//  END DO j

	  h[i] = func;

	  givef (f, &h[i], &s, val, clint, /* pd, */ stt, tcov,
		 thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint, /* dmax, */
		 ltmax, ltmin, nclass, numa, numo, &msfail, maxjb,
		 mtest, false, &params->results);

	  neval++;

/*
*
*     All points in the initial simplex become output if IPRINT=1.
*/
	  if (iprint==1) {
	    DUMP_SS(h[i], f)
	  }

	}					//  END DO Loop_730

/*
*
*     Now follows the basic loop.  That is, given a simplex, it
*     determines a new simplex and tests for convergence as
*     required (following the flow chart given in Nelder and
*     Mead).
*
*     HMAX and HMIN are the maximum and minimum function values
*     of the current simplex.
*
*/
        bool continue_search = true;
        while (continue_search) {

            for (int loop = 0; loop < nloop; ++loop) {

                int hmax_index = 0;
                int hmin_index = 0;
                double hmax = h[0];
                double hmin = h[0];

/* Loop_750: */
	for (i=1; i < np1; i++) {	//  DO i=2,np1

	  if (h[i] > h[hmax_index]) {
	    hmax_index = i;
	    hmax = h[i];
	  }

	  if (h[i] < h[hmin_index]) {
	    hmin_index = i;
	    hmin = h[i];
	  }

	}				//  END DO Loop_750

/*
*     The centroid of all vertices, excluding the maximum, is
*     now found.
*/

/* Loop_790: */
	for (i=0; i < nop; i++) {	//  DO i=1,nop
	  pbar[i]=0.0;
	}				//  END DO Loop_790


/* Loop_800: */
	for (i=0; i < np1; i++) {		//  DO i=1,np1
	  if (i != hmax_index) {
Loop_810:   for (j=0; j < nop; j++) {	//  DO j=1,nop
	      pbar[j] += g[i][j]/(nap);
	    }				//  END DO Loop_810
	  }
	}					//  END DO Loop_800

/*
*     The program reflects the maximum through PBAR to PSTAR, and
*     evaluates the function at PSTAR (to give HSTAR).
*/

/* Loop_820: */
	for (i=0; i < nop; i++) {		//  DO i=1,nop
	  pstar[i] = a * (pbar[i]-g[hmax_index][i]) + pbar[i];
	}					//  END DO Loop_820

	hstar = 0.0;

	givef (pstar, &hstar, &dcoeff, val, clint, /* pd, */ stt,
	       tcov, thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint,
	       /* dmax, */ ltmax, ltmin, nclass, numa, numo,
	       &msfail, maxjb, mtest, false, &params->results);

/*
*     The next 5 statements test whether a progress report is
*     required and, if so, provide one.  This procedure occurs
*     frequently in the program.
*/

	neval++;

/*
*     If the number of function evaluations to date exceeds 750
*     (=MAX), the program prints out parameter values provided that
*     IPRINT and JPRINT have been set at 1.
*/
		  
	if ((neval > max) && (iprint == 1) && (jprint == 1)) {
	  DUMP_SS(hstar,pstar)
	}

                bool take_pstar = false;
                if (hstar < hmin) {
/*
*     If HSTAR is less than HMIN, PBAR is reflected through PSTAR
*     (to give PSTST) and the function is evaluated there (to
*     give HSTST).
*/

/* Loop_830: */
	for (i=0; i < nop; i++) {		//  DO i=1,nop
	  pstst[i] = c * (pstar[i]-pbar[i]) + pstar[i];
	}					//  END DO Loop_830

        hstst = 0.0;

	givef (pstst, &hstst, &dcoeff, val, clint, /* pd, */ stt,
	       tcov, thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint,
	       /* dmax, */ ltmax, ltmin, nclass, numa, numo,
	       &msfail, maxjb, mtest, false, &params->results);

/*
*     If IPRINT=1 the program prints out the progress of the
*     iteration.  This is not normally required.
*/
        neval++;

	if ((iprint == 1) && (jprint == 1)) {
	  DUMP_SS(hstst,pstst)
	}

/* Line_870: */
                    if (hstst >= hmin) {
                        take_pstar = true;
                    } else {

/*
*     If HSTST is less than HMIN, the maximum point of the current
*     simplex is replaced by PSTST and HMAX is replaced by HSTAR,
*     then a test is performed.
*/

/* Loop_880: */
	for (i=0; i < nop; i++) {		//  DO i=1,nop
	  g[hmax_index][i] = pstst[i];
	}					//  END DO Loop_880

                        h[hmax_index] = hstst;
                        continue;
                    }
                    
                } else {
/*
*
*     If HSTAR is not less than HMIN, the program tests is HSTAR
*     is greater than the function value at all vertices other
*     than the maximum one.
*/
/* Loop_900: */
	for (i=0; i < np1; i++) {		//  DO i=1,np1
                    if (i != hmax_index && hstar < h[i]) {
                        take_pstar = true;
                        break;
                    }
    }					//  END DO Loop_900
                }

                if (take_pstar) {
/* Loop_1030: */
                    for (i=0; i < nop; i++) {		//  DO i=1,nop
                        g[hmax_index][i] = pstar[i];
                    }					//  END DO Loop_1030
                    h[hmax_index] = hstar;
                    continue;
                }

/*
*     If it is less than at least one of these vertices, the
*     maximum point of the current simplex is replaced by PSTAR and
*     HMAX by HSTAR.  A test is then performed.
*
*     If HSTAR is greater than all the function values excluding
*     the maximum, the program tests if HSTAR is greater than HMAX.
*     If it is not, the maximum point of whichever simplex is
*     now in store (depending on whether HSTAR is greater or less
*     than HMAX) is replaced by PSTAR and HMAX by HSTAR, and the
*     contracted point PSTST and the function value there (HSTST)
*     are calculated.
*/

	if (hstar <= hmax) {
	  for (i=0; i < nop; i++) {		//  DO i=1,nop
	    g[hmax_index][i] = pstar[i];
	  }					//  END DO
	  hmax = hstar;
	  h[hmax_index] = hstar;
	}

/* Loop_930: */
	for (i=0; i < nop; i++) {		//  DO i=1,nop
	  pstst[i] = b*g[hmax_index][i] + (1.0-b)*pbar[i];
	}					//  END DO Loop_930

        givef (pstst, &hstst, &dcoeff, val, clint, /* pd, */ stt,
	       tcov, thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint,
	       /* dmax, */ ltmax, ltmin, nclass, numa, numo,
	       &msfail, maxjb, mtest, false, &params->results);

	neval++;
	if ((iprint == 1) && (jprint == 1)) {
//          WRITE (2,720) neval,hstst,(pstst[j],j=0,nop-1)	/* JMB: loop range has been set, was 1..nop */
	  DUMP_SS(hstst,pstst)
	}

/* Line_970: */

/*
*     If HSTST is not greater than HMAX, the maximum point is replaced by
*     PSTST and HMAX by HSTST.  A test is then applied.
*/

	if (hstst <= hmax) {
	  for (i=0; i < nop; i++) {		//  DO i=1,nop
	    g[hmax_index][i] = pstst[i];
	  }					//  END DO
	  h[hmax_index] = hstst;
	  continue;
	}
/*
*     If HSTST is greater than HMAX, each point in the current
*     simplex is replaced by a point midway between its current
*     position and the position of the minimum point of the
*     current simplex.  The function is evaluated at each new
*     vertex and the test performed.
*/

/* Loop_990: */
	for (i=0; i < np1; i++) {		//  DO i=1,np1,1
	  for (j=0; j < nop; j++) {	//  DO j=1,nop,1
	    g[i][j] = (g[i][j]+g[hmin_index][j])/2.0;
	  }				//  END DO
	}					//  END DO Loop_990

/*
*     The conditional loop Inner_1010 is the first of several
*     introduced to ensure the program uses supplied f values
*     wherever the stepsize is set at zero (to avoid a processor
*     error?).
*/

/* Loop_1020: */
	for (i=0; i < np1; i++) {		//  DO i=1,np1,1
	  for (j=0; j < nop; j++) {	//  DO j=1,nop,1
	    f[j] = g[i][j];
	    if (step[i] < approx) f[i] = ft[i];	/* Inner_1010 */
	  }                                         //  END DO Loop_1010

          givef (f, &h[i], &dcoeff, val, clint, /* pd, */ stt,
		 tcov, thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt,
		 kprint, /* dmax, */ ltmax, ltmin, nclass, numa, numo,
		 &msfail, maxjb, mtest, false, &params->results);

	  neval++;
	  if ((iprint == 1) && (jprint == 1)) {
	    DUMP_SS(h[i],f)
          }
	}					//  END DO Loop_1020

/*
* End of the basic loop
*/
            }

/*
*   Tests for Convergence -
*
*     The mean and standard deviation of the function values of the
*     current simplex are now calculated.
*/

	hstd = 0.0;
	hmean = 0.0;

Loop_1060:
	for (i=0; i < np1; i++) {		//  DO i=1,np1
	  hstd += h[i]*h[i];
	  hmean += h[i];
	}					//  END DO Loop_1060

	hmean = hmean/(np1);
	hstd = (hstd-(np1)*hmean*hmean) / (np1);
	hstd = sqrt(fabs(hstd));

/*
*     The parameter values (F) at the centroid of the current
*     simplex and the function value there (FUNC) are now
*     calculated.
*
*/

Loop_1080:
	for (i=0; i < nop; i++) {
	    if (step[i] < approx) {
		f[i] = ft[i];
	    }
	    else {
		f[i] = 0.0;
		for (j=0; j < np1; j++) {
		    f[i] += g[j][i];
		}
		f[i] = f[i]/np1;
	    }
	}

        givef (f, &func, &dcoeff, val, clint, /* pd, */ stt, tcov,
	       thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint, /* dmax, */
	       ltmax, ltmin, nclass, numa, numo, &msfail, maxjb,
	       mtest, false, &params->results);

	neval++;

/*
*     If the number of evaluations has exceeded the value of MAX
*     set (=750), the convergence process is judged not to have
*     succeeded in this case.  If so, this particular run of Loop
*     1410 is assumed to have yielded a 'no result'.  A counting
*     variable (MFAIL) is given a value of 1, to be subtracted from
*     the value of MAXJB later, and the next run through the loop
*     begins.
*/

	if (neval >= max) {

/*
*     MFAIL is a parameter which counts the number of times a series
*     of iterations failed to converge on a minimum within the set
*     maximum number of iterations.  Parameter values are set at
*     zero so that the failed bootstraps do not contribute to the
*     overall estimates.
*/
/* Line_1090: */
	  mfail++;
	  den[bootstrap] = 0;
	  coeff1[bootstrap] = 0;
	  coeff2[bootstrap] = 0;
	  dcoeff = 0;

	  if ((iprint != 1) && (jprint != 1)) goto Line_1380;

	  fprintf(output_results, " NUMBER OF FUNCTION EVALUATIONS EXCEEDS %4i\n", max);
	  fprintf(output_results, " STANDARD ERROR OF FUNCTION VALUES OF LAST SIMPLEX %13.6e\n", hstd);
	  fprintf(output_results, "  CENTROID OF LAST SIMPLEX  %13.5e%13.5e%13.5e%13.5e\n", f[0], f[1], f[2], f[3]);
#if (NUM_SHAPE_PARAMS != 4)
#error "Output format does not match NUM_SHAPE_PARAMS"
#endif
	  fprintf(output_results, "  FUNCTION VALUE AT CENTROID   %13.6e\n", func);

	  goto Line_1375;
	}

/* Line_1150: */

	if (iprint == 1) {
	  DUMP_SS(func,f)
	}

/*
*     If the standard deviation calculated above is not less than
*     the criterion set (STOPC) and IFLAG are set to zero and the
*     basic loop begins again.
*/
	  if (hstd >= stopc) {
	  iflag = 0;
	  continue;
	}

	if ((iprint == 1) && (jprint == 1)) {
	  fprintf(output_results," *\n  INITIAL EVIDENCE OF CONVERGENCE\n");
	  fprintf(output_results,"  CENTROID OF LAST SIMPLEX  %13.5e%13.5e%13.5e%13.5e\n", f[0],f[1],f[2],f[3]);
#if (NUM_SHAPE_PARAMS != 4)
#error "Output format does not match NUM_SHAPE_PARAMS"
#endif
	  fprintf(output_results,"  FUNCTION VALUE AT CENTROID   %13.6e\n", func);
	}

/*
*     If the standard deviation is less than the stopping
*     criterion, IFLAG is set =0 if there was no evidence of
*     convergence on the last test and =1 if there was evidence
*     of convergence.
*
*     If IFLAG=0, reset IFLAG=1,  the mean of the function
*     values of the current simplex are saved (as SAVEMN), and
*     computation goes back to the beginning of the basic loop.
*/

Line_1210:
	if (iflag <= 0) {
	  iflag = 1;
	  savemn = hmean;
	  continue;
	}

/*
*     If IFLAG=1, the program tests if the change in the mean is
*     less than the stopping criterion.  If it is, the process is
*     said to have converged.  If not, IFLAG is set at zero and
*     computation reverts to the start of the basic loop.  STOPC
*     and TEST values tested empirically.
*/

	if (hmean != 0) {
	  test = savemn / hmean;
	  if ((test <= 0.9999995) || (test >= 1.0000005)) {
	    iflag = 0;
	    continue;
	  }
	}

/*
*     None of the continuation criteria have been met. Conclude the
*     search of this bootstrap.
*/
            continue_search = false;
        }

/*
*     If JPRINT=1 the program prints out the results of each successful
*     convergence on a minimum.
*/

Line_1250:
	if ((jprint == 1) && (maxjb != 1)) {

	  fprintf(output_results, "\n\n Process converges on minimum after %4i function evaluations\n\n", neval);
	  fprintf(output_results," Minimum at    %13.6e %13.6e %13.6e %13.6e\n", f[0],f[1],f[2],f[3]);
#if (NUM_SHAPE_PARAMS != 4)
#error "Output format does not match NUM_SHAPE_PARAMS"
#endif
	  fprintf(output_results, "\n Minimum function value   %13.6e\n", func);
	  fprintf(output_results, "\n End of Search\n *************\n");
        }

/*
*     MTEST is set at 1 to flag that convergence has occurred.  This
*     is carried into Subroutine GIVEF to trigger computation of
*     MSFAIL where computation of S cannot occur.
*/

Line_1320:
	mtest = 1;

/*
*     Program execution returns to the subroutine to yield final
*     values of F[0], F[1] and F[2].
*
*/
	givef (f, &func, &dcoeff, val, clint, /* pd, */ stt, tcov,
	       thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint, /* dmax, */
	       ltmax, ltmin, nclass, numa, numo, &msfail, maxjb,
	       mtest, false, &params->results);


	kprint = 0;

/*
*     The estimated 'best fit' values of the parameters from the current
*     pass through Loop 1410 are now computed, tabulated and printed out.
*
*     This involves dividing by the same variables used earlier to
*     compare calculated with observed values in Subroutine GIVEF.
*
*     A density estimate (DEN) is calculated from D2L=F[2] by
*     correcting units to no./ha, and dividing by 2L in the case of
*     line transect data, and by 2Vt in the case of fixed-point data.
*/

	if (km > 0) {
	  den[bootstrap] = (1.e6*f[2])/(sns*dist*pd);
	}
	else {
	  if (ifx == 0) {
	    den[bootstrap] = (1.e4*f[2])/(sns*dist*pd);
	  } else {
	    den[bootstrap] = (1.e4*f[2])/(2.*rate*durn*ps*pd);
	  }
	}

/*
*     Values of the other parameter estimates are obtained by redefining
*     the estimates of F[0] and F[1], thus:
*/
/* Line_1370: */

	coeff1[bootstrap] = f[0];
	coeff2[bootstrap] = f[1];
	coeff3[bootstrap] = dcoeff;

/*
*      The density estimate is transformed to TRDEN by a logarithmic
*      transformation, then a running total TTRDEN is calculated.
*/

	trden[bootstrap] = log(1 + den[bootstrap]);
	ttrden += trden[bootstrap];

/*
*     Running totals of DEN, COEFF1 and COEFF2 are now made to
*     enable calculation of mean values after the loop ends.
*/

Line_1375:
	tden += den[bootstrap];
	tcoeff1 += coeff1[bootstrap];
	tcoeff2 += coeff2[bootstrap];
	tcoeff3 += coeff3[bootstrap];

	if (maxjb == 1) goto Line_1468;

/*
*     Loop 1410 now ends, returning calculations to the start until the
*     maximum preset number of bootstraps (MAXJB) value is reached.
*     G[I][J] values are reset to zero to enable a new series of
*     iterations in the basic loop to begin again.
*
*/

Line_1380:
	iflag = 0;
	kprint = 0;

/* Loop_1390: */
	for (i=0; i < np1; i++) {		//  DO i=1,np1
	  for (j=0; j < nop; j++) {	//  DO j=1,nop
	    g[i][j] = 0.0;
	  }				//  END DO
	}					//  END DO Loop_1390

/*
*
*     Other F values are reset to their original values.
*/

/* Loop_1400: */
	for (ie=0; ie < nop; ie++) {		//  DO ie=1,nop
	  f[ie] = ft[ie];
	}					//  END DO Loop_1400

/* ==============================================================
*  ======                                                  ======
*  ======               End of the main loop               ======
*  ======                                                  ======
*  =============================================================*/
        
      }								//  END DO Loop_1410

/*
*
*     Following completion of the runs through Loop 1410, the means and
*     standard errors of each of the parameters are now calculated,
*     based on the values of the three arrays DEN[bootstrap],
*     COEFF1[bootstrap] and COEFF2[bootstrap].
*
*     The overall means ESTDEN, COEFFNT1 and COEFFNT2 are calculated
*     first.  NUMEST is the number of parameter estimations made.
*/
	
      numest = maxjb - mfail;
      if (numest == 0) numest = 1;

/*
*     Each of the following variables is an estimate of a parameter
*     mean value, beginning with the population mean ESTDEN, then
*     the other parameters - conspicuousness coefficient COEFFNT1,
*     cover proportion or attenuation coefficient COEFFNT2 - and the
*     detectability coefficient COEFFNT3, in sequence.
*/

      estden = tden/numest;
      double coeffnt1 = tcoeff1/numest;
      double coeffnt2 = tcoeff2/numest;
      double coeffnt3 = (numest == msfail) ? tcoeff3 : (tcoeff3/(numest-msfail));

/*
*     A mean of the transformed values, TTRDENMN, is calculated too.
*/
      float ttrdenmn = ttrden/numest;

/*
*
*     The next step is to calculate the standard errors of each
*     parameter, provided that the number of analyses exceeds 1.
*     Each is the standard deviation of the parameter estimates.
*     If NUMEST is 0 or 1, standard error calculation is bypassed.
*/
      double sden = 0, scf1 = 0, scf2 = 0, scf3 = 0;

      if (numest > 1) {

      double dsum = 0.0;
      double cf1sum = 0.0;
      double cf2sum = 0.0;
      double cf3sum = 0.0;
      float tdsum = 0.0;

      for (bootstrap=0; bootstrap < maxjb; bootstrap++) {
	if ((den[bootstrap]) != 0.) {
	  double dendif = den[bootstrap]-estden;
	  dsum += dendif*dendif;
	}

        if (coeff1[bootstrap] != 0.0) {
	  double cf1dif = coeff1[bootstrap]-coeffnt1;
	  cf1sum += cf1dif*cf1dif;
        }

	if (coeff2[bootstrap] != 0.0) {
	  double cf2dif = coeff2[bootstrap]-coeffnt2;
	  cf2sum += cf2dif*cf2dif;
	}

	if (coeff3[bootstrap] != 0.0) {
	  double cf3dif = coeff3[bootstrap]-coeffnt3;
	  cf3sum += cf3dif*cf3dif;
	}

/*
*     A standard deviation of the transformed means, STRDEN, is now
*     computed.
*/
	if (trden[bootstrap] != 0.0) {
	  float tdendif = trden[bootstrap]-ttrdenmn;
	  tdsum += tdendif*tdendif;
	}
      }

      sden = sqrt(dsum/(numest-1));
      scf1 = sqrt(cf1sum/(numest-1));
      scf2 = sqrt(cf2sum/(numest-1));
      if ((numest-msfail-1) > 0) scf3 = sqrt(cf3sum/(numest-msfail-1));
      float strden = sqrt(tdsum/(numest-1));

/*
*     The p=0.05(2) distribution of t with sample number is now
*     approximated by a function t05.
*/

      double t05 = 1.97623 + 5.95434/pow(nvals, 1.4);

/*
*     The 95% lower (CL1) and upper (CL2) confidence limits are now
*     calculated from the transformed mean and standard deviation,
*     then back-transformed to the original units.
*/

      tcl1 = ttrdenmn - t05*strden;
      tcl2 = ttrdenmn + t05*strden;
      cl1 = exp(tcl1) - 1.0;
      cl2 = exp(tcl2) - 1.0;
      }

/*
*     If all the STEP values supplied have been set at zero, several
*     parameters need to be given values that will correspond to
*     those where at least one parameter was estimated.
*/

Line_1470:
      if (nap <= 0) {
	numest = 1;
	coeffnt1 = f[0];
	coeffnt2 = f[1];
      }

/*
*     The products of the program are now output.
*
*     The number of animal groups in the selected range (NGROUPS) is
*     first, followed by the estimated topographical cover (TCOV),
*     then the number of actual parameter estimations made (NUMEST).
*/

Line_1468:

//      WRITE (2,1471) ngroups
// 1471 FORMAT (/,' Number of Groups in Distance Range = ',i4)

      fprintf(output_results, "\n Number of Groups in Distance Range = %4i\n", ngroups);

      if (ifx == 0) {

//        WRITE (2,1472) numain
// 1472   FORMAT (' Number of Individuals Detected Ahead = ',i5)

	fprintf(output_results," Number of Individuals Detected Ahead = %5i\n", numain);
      } else {

//        WRITE (2,1473) numain
// 1473   FORMAT (' Number of Individuals Detected =',i5)

	fprintf(output_results, " Number of Individuals Detected = %5i\n", numain);
      }

      if (ifx == 0) {
	fprintf(output_results, " Number Overtaking (distance unmeasured) = %4i\n", numoin);
      }

/* Line_1475: */

      fprintf(output_results, " Height Difference from Eyelevel = %5.1f m\n", thh);

      if (ifx == 0) {
        fprintf(output_results, " Movement Correction Factor (J) = %6.3f\n", start.estj);
      }

      if (ifx == 0) {
	fprintf(output_results, " Adjusted Transect Length (LJ) = %11.3f m\n", dist);
      }

      fprintf(output_results, " Topographical Cover Value = %6.4f\n", tcov);

      fprintf(output_results, " Maximum Detection Distance (%s) = %7.1f m\n",
              dmaxIsPreset ? "preset" : "estimated",
              estdmax);

      if (maxjb == 1) {
	fprintf(output_results, "\n\n Process converges on minimum after %4i function evaluations\n", neval);

/*
*     Subroutine qsf is now called if MAXJB=1.
*/

	qsf (f, func, approx, dist, durn, rate, step, stopc,
	     &s, val, clint, stt, tcov, thh, /* vgh, */ /* imv, */ ishow, &kprint, ltmax,
	     ltmin, nclass, numa, numo, &msfail, mtest, &estden, &sden, &hstst,
	     pd, ps, ifx, iprint, iry, kdt, km, nap, /* neval,*/ nop, ns, /* nvals, */
	     np1, g, h, maxjb, /* jprint, */ &params->results);
	iqsf = true;
	goto Line_1920;
      }

      fprintf(output_results, " Number of Parameter Estimations =%4i\n", numest);

/*
*     Now follow general headings for the results table.
*
*/
      fprintf(output_results, "\n Estimated Parameter Values:\n\n");
      kprint = 1;

	
      fprintf(output_results,     " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " x    Parameter     x      Value       x  Standard Error  x       Unit       x\n");
      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " x ESTIMATED        x                  x                  x                  x\n");

/*
*     Density estimates are printed.
*/
	
      if (km == 0) {
	if (maxjb > 1) {
	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x    %10.3f    x  indivs/hectare  x\n", estden,sden);
	}
	else {
	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x (indeterminate)  x  indivs/hectare  x\n", estden);
	}
      }
      else {
	if (maxjb > 1) {
	  fprintf(output_results, " x DENSITY  (D)     x    %10.3f    x    %10.3f    x   indivs./sq.km. x\n", estden,sden);
	}
	else {
	  fprintf(output_results, " x DENSITY  (D)     x    %10.2f    x (indeterminate)  x   indivs./sq.km. x\n", estden);
	}
      }

/*
*     The conspicuousness coefficient is next printed.
*/

Line_1650:

      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " x Conspicuousness  x                  x                  x                  x\n");

      if (maxjb == 1) {

	fprintf(output_results,   " x Coefficient  (a) x    %10.4f    x (indeterminate)  x      metres      x\n", coeffnt1);
      }
      else {

	fprintf(output_results,   " x Coefficient  (a) x    %10.3f    x    %10.3f    x      metres      x\n", coeffnt1, scf1);
      }

/*
*     The second coefficient is either a cover proportion or a sound
*     attenuation coefficient, decided by the values of KDT.
*/

Line_1710:
      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
      fprintf(output_results,     " x                  x                  x                  x                  x\n");


      if (kdt == 1) {
	fprintf(output_results,   " x Attenuation      x                  x                  x                  x\n");
	if (maxjb == 1) {
	  fprintf(output_results, " x Coefficient  (b) x    %10.4f    x (indeterminate)  x    per  metre    x\n", coeffnt2);
	}
	else {
	  fprintf(output_results, " x Coefficient  (b) x    %10.4f    x    %10.4f    x    per metre     x\n", coeffnt2, scf2);
	}
      }
      else {
	fprintf(output_results,   " x Cover            x                  x                  x                  x\n");
	if (maxjb == 1) {
	  fprintf(output_results, " x Proportion   (c) x    %10.4f    x (indeterminate)  x                  x\n", coeffnt2);
	}
	else {
	  fprintf(output_results, " x Proportion   (c) x    %10.4f    x    %10.4f    x                  x\n", coeffnt2, scf2);
	}
      }

/*
*     The table is now ruled off and some other details printed.
*
*/

Line_1840:

      fprintf(output_results,     " x                  x                  x                  x                  x\n");
      fprintf(output_results,     " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n");


/*
*     The next two parameters are output only if at least one parameter
*     is being estimated.
*/
	
      if (nap >= 1) {

	fprintf(output_results, " 95%% confidence limits for density estimate:%10.3f %10.3f\n", cl1, cl2);
	fprintf(output_results, " Detectability Coefficient (S) =%9.2f, SE =%9.2f\n", coeffnt3, scf3);
      }

/*
*     Key model estimates are now output if ISHOW was originally set
*     at 1.  The best estimates of F[0], F[1] and F[2] must be entered
*     in the subroutine GIVEF to do this.  The totals in the class
*     intervals are those from the initial data, derived from the
*     saved array variable VALT[].  The output is then produced
*     within the subroutine.  The same is done with NUMO and NUMA.
*
*     f[0] and f[1] are left at their original values if all parameters
*     were preset (STEP=0).
*/

      if (nap > 0) {
	f[0] = coeffnt1;
	f[1] = coeffnt2;
      }

/* Loop_1870: */

      for (jv=0; jv < nclass; jv++) {		//  DO jv=1,nclass
        val[jv] = valt[jv];
      }						//  END DO Loop_1870

      numo = numoin;
      numa = numain;

/*
*    Calculated density values are now multiplied by observing situation
*    variables for a last comparison with observed values within
*    Subroutine GIVEF.
*/
/* Line_1890: */

      if (ifx == 0) {

         if (km == 0) {
            f[2] = (estden*dist*pd*sns) / 1.e4;
         } else {
            f[2] = (estden*dist*pd*sns) / 1.e6;
         }
      }
      else {
         f[2] = (estden*2.*ps*rate*durn*pd) / 1.e4;
      }


Line_1920:
      givef (f, &func, &dcoeff, val, clint, /* pd, */ stt, tcov,
	     thh, /* vgh, */ ifx, /* imv, */ iry, ishow, kdt, kprint, /* dmax, */
	     ltmax, ltmin, nclass, numa, numo, &msfail, maxjb,
	     mtest, iqsf, &params->results);

    fclose(output_results);	// CLOSE (unit=2)

/*
* And finally, write a file which tabulates observed and calculated frequencies.
*/
    double tr = clint*nclass + stt;
    if (f[3] > thh) {
        double model_limit = sqrt(f[3]*f[3] - thh*thh);
        if (kdt > 1) {
            model_limit = fmin(kdt, model_limit);
        }
        if (tr > model_limit) {
            nclass = (model_limit - stt + clint) / clint;
        }
    }

      if (kprint > 0) {
          FILE *output_graph = fopen(graphfile, "w");
          if (output_graph == NULL) {

              perror("Error opening graph file \'%s\': ");

          } else {

              /* Graph file is open */

              fprintf(output_graph, "      Midpt.   Calculated     Observed   \n");

              for (int jv = 0; jv < nclass; ++jv) {
                  fprintf(output_graph,
                          "     %7.1f      %7.2f       %6.1f    \n",
                          params->results.midpoints[jv],
                          params->results.calcn[jv],
                          params->results.obsdn[jv]);
              }
              fclose(output_graph);
          }

      }  /* end (kprint > 0) */

      params->results.num_intervals = nclass;
      params->results.km = km;
      params->results.estden = estden;
      params->results.sden = sden;
      params->complete = 1;
}
//      END SUBROUTINE calculate_density
