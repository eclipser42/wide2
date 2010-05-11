!     PROGRAM WildlifeDensity  
!
!     (File mnps2.f, Version 0.9.8)
!
!     This program is designed to return population density estimates
!     from 'distance' data collected using either line transect or fixed
!     point methods, together with estimates of other parameters of the
!     observing situation.  It is intended primarily for use in
!     ground surveys of mammal and bird populations in terrestrial
!     habitats, and with aerial surveys using helicopters.  It is
!     capable, with modification, of being used with aerial survey
!     data collected from fixed-wing aircraft.
!
!     WildlifeDensity works by mathematical modelling of the frequency
!     distributions of animal numbers detected under uniform observing 
!     conditions over a range of distances from either an observer
!     or a transect line.  Curve-fitting is achieved by using the 
!     simplex method to seek parameter values which minimize the ordinary
!     least squares fit between the observed and calculated numbers across
!     the range of detection distances.
!
!     Based on the 'best-fit' solution so obtained, the program returns
!     population density estimates, estimates of other parameters of the
!     model and related information on the system in which the original
!     observations were made. Standard errors of the parameters are
!     obtained by 'bootstrap' resampling of data from the original data 
!     set, unless the number of sets of iterations is set at 1, when the
!     program computes standard errors by a quadratic surface fitting
!     method (qsf).
!
!
!     The WildlifeDensity model was designed by David Morgan, Department
!     of Zoology, The University of Melbourne, Vic. 3010, Australia, and
!     engineered by James Clough.
!
!
!     ******************************************************************
!
!
!     File mnps2.f has been written in Fortran, and consists of a 
!     main subroutine, SUBROUTINE calculate_density, and five subsidiary
!     subroutines: SUBROUTINE givef, SUBROUTINE resample, SUBROUTINE 
!     srandnag (iseed), SUBROUTINE srandgen (result) and SUBROUTINE qsf.
!
!
!     The main subroutine determines parameter values for a least
!     squares fit between calculated and observed frequencies,
!     using the simplex method, a numerical method for minimizing an 
!     objective function in a many-dimensional space.  Its design is 
!     derived from Nelder & Mead (1965) Computer Journal 7:308-313, as 
!     used in the 'MINIM' program by D.E.Shaw, Divn. of Mathematical 
!     Statistics, CSIRO, Australia, though extensively modified since. 
!     The quadratic surface fitting method within it is based on a Hessian 
!     matrix, as set out in the appendix to the Nelder and Mead paper.
!
!     Subroutine GIVEF calculates the sum of squares by using the
!     mathematical model, and returns this to the main program.
!     Subroutines SRANDNAG and SRANDGEN supply a pseudo-random number used 
!     in the data resampling procedure in response to a seed value supplied
!     by the main program.
!
!
!     The program computes density estimates based on the numbers of
!     animals detected either at various horizontal radial distances (r)
!     from an observer or at various horizontal distances (y) perpendicular
!     to a transect line.  It accepts several alternative types of distance
!     data, as determined by appropriate values of the control parameters
!     IFX and IRY supplied with the data, viz:
!
!       . radial distance data supplied (IFX=0 and IRY=0);
!       . perpendicular distances calculated from radial distance
!           and horizontal detection angle data supplied
!           (IFX=0 and IRY=1);
!       . perpendicular distance data supplied (IFX=0 and IRY=2);
!       . fixed-observer distance data supplied (IFX=1 and IRY=0).
!
!
!     Data are entered into the program by means of an appropriately
!     formatted tab- or comma-separated dataset in a file named
!     '<filename>.dat'.  It outputs two files:  the first, named
!     '<filename>.results', prints out final values of the parameters
!     of the dataset supplied, together with some of the main attributes
!     of the observing situation;  the second, named '<filename>.graphData',
!     sets out the original frequency distribution and the calculated
!     frequency distribution of best fit in a column format suitable
!     as input to a graphing application such as MS 'Excel'.
!
!
!     The parameters supplied in the input via the data file are, in 
!     sequence:
!
!
!     A header line - to identify and label the data set.
!
!     NVALS - the total number of individual animal detections
!          in the sample submitted for modelling.
!
!     CLINT - the class interval to be used within the program,
!          chosen so that 80 x CLINT is at least equal to the 
!          maximum radial detection distance.
!
!     STT - the detection distance value (r or y) from which class
!          intervals begin (usually zero).
!
!     NUMA - the total number of individual animals detected
!          ahead of a moving observer on a transect.
!
!     NUMO - the number of individuals overtaking the observer
!          from behind during a transect, obtained by totalling
!          the number of individuals in R=0 observations.
!
!     DIST - the overall transect length (L), expressed either in
!          metres or kilometres, or the distance corrected for the 
!          effect of animal movement (LJ).  (If F(4) is set at
!          zero, only the overall transect length is needed and
!          the program calculates its own LJ value.)
!
!     THH - the vertical distance between observer eye level
!          and the median horizontal plane occupied by
!          the population, best represented by the root mean
!          square of the measured vertical height differences.
!
!     LTMIN - the approximate minimum distance at which an animal
!          may be obscured from the observer by topography. Set
!          at 999 if the topography is approximately level.
!
!     LTMAX - the p=0.001 upper-limit distance at which an animal
!          may be detected by the observer in uneven topography
!          in the absence of vegetation cover.  Set at 0 if this
!          distance is to be calculated from the data.
!
!     IFX - a parameter to modify the computations appropriately
!          where data are collected by an observer sampling from
!          a fixed point, viz:
!
!            = 0  handles data from line transects;
!            = 1  handles data from fixed point sampling.
!
!     IRY - a parameter to control computations to model different 
!           types of line transect data, viz:
!
!            = 0  models radial distance transect data N(r);
!            = 1  models perpendicular distance data N(y), using
!                   radial distance and angle data supplied;
!            = 2  models perpendicular distance data N(y), using
!                    perpendicular distance data supplied.
!
!     NS - the number of sides of the transect line used for
!          observations during a line transect survey, viz:
!
!            = 0  not transect data (fixed point), so not relevant;
!            = 1  handles data from only one side of the line;
!            = 2  handles data from both sides of the transect line.
!
!     KM - a parameter to modify the program if transect lengths and
!          detection distances are measured in kilometres, viz:
!
!            = 0  all distances are measured in metres;
!            = 1  transect lengths are measured in kilometres
!                   and detection distances in metres;
!            = 2  both transect lengths and detection distances
!                 are measured in kilometres.
!
!     KDT - a control parameter to indicate attributes of the
!          observational data supplied to the program, viz:
!
!            = 0  for visual and flushing data;
!            = 1  for auditory data or long distance (no vegn.) data;
!            > 1  the maximum class interval boundary distance (in m) 
!                 when visual data are from a limited distance range.
!
!     IPRINT - a parameter which directs the program to print out
!          progress sum of squares and parameter values to allow
!          perusal of the progress of the minimization process.
!
!            = 0  no progress evaluations;
!            = 1  If JPRINT=0, reports initial convergence and the
!                 function and parameter values of the initial
!                 simplex from each bootstrap iteration series;
!                 If JPRINT=1, directs printing of these values at
!                 individual steps through the function minimization 
!                 process, and informs the  user if the process
!                 failed to converge on a minimum within 750
!                 steps (the limit set).
!
!     JPRINT - a parameter to enable perusal of the original and
!           bootstrapped frequency distributions, viz:
!       
!            = 0  no frequency distributions or internal output;
!            = 1  prints the frequency distributions of the original
!                 and bootstrapped data and makes possible perusal
!                 of some of the function minimization process.
!
!     ISHOW - a parameter which prints out the values of the
!          following functions, in sequence, for each class interval
!          across the distance range for the final result: PRC,
!          PRR, QR, E, TOTE (total E), ED, EXPD, EXPDV and the
!          detectability coefficient S (see comments within
!          Subroutine GIVEF.
!
!            = 0  no output;
!            = 1  produces output.
!
!     MAXJB  -  the number of sets of iterations used to derive
!          the backstrapped series of model parameters and
!          calculate standard errors.  This is best set at between
!          500 and 2000.  Fewer sets saves computation time but
!          increases the size of confidence intervals of parameters.
!          If MAXJB is set at 1, no bootstrapping occurs; instead
!          the program calls Subroutine qsf and calculates variances
!          and standard errors using surface fitting algorithms.
!
!     VGH - the approximate average height of vegetation cover
!          in the animal's habitat in situations where the observer
!          is well above the plane of the population and most of
!          the line of detection is unobstructed.
!
!     DURN - the duration of a fixed-point census (in min.).
!
!     RATE - the overall mean rate of animal movement (in m/min)
!
!     PD - the proportion of the population observable at the
!          time of a census, usually taken =1.
!
!     PS - in fixed-point censuses, the proportion of the circle
!          scanned by the observer during the census (usually
!          either 0.5 or 1).
!
!
!     F  - the model's parameter values used in the search for the
!          minimum point.  Four initial F values are required
!          as input to the program:
!
!          F(1):  an estimate of the conspicuousness coefficient
!                 of the species under the conditions of the
!                 census;
!
!          F(2):  an estimate of either the overall lateral vegetation
!                 cover (c) between observer and animals, or of the
!                 attenuation coefficient (b) of animal sounds in the 
!                 observing situation under the census conditions;
!
!          F(3):  an estimate of the population density (D),
!                 expressed in number of individuals per hectare
!
!          F(4):  an estimate of the maximum distance (in m) from the
!                 observer at which species recognition is
!                 possible under the conditions of the census.
!
!          Those supplied with the data are the starting values
!          for the search; those printed on exit are the parameter
!          values which specify the minimum point.  Initial F()
!          values supplied are best based on existing knowledge.
!
!
!     STEP - the step sizes used initially in the function
!          minimization process to modify the F values. Typical
!          STEP values are:
!
!          STEP(1):  just under 50% of the initial F(1) value;
!
!          STEP(2):  equal to the initial F(2) value;
!
!          STEP(3):  half the initial F(3) value;
!
!          STEP(4):  usually set at zero.
!
!          Setting the STEP value at 0.0 for a parameter fixes an F()
!          value at that supplied to the program.  Where a data set
!          (NVALS) is small (say, < 50 detections), STEP(2) should be
!          set at 0.0 .  Where the set is very small (say, < 30 
!          detections), STEP(1) should also be set at 0.0 .
!
!
!    R() - the radial detection distances (r) originally measured
!          in the field, submitted in the order in which the data
!          were collected.  Or the perpendicular distances from the
!          transect line if these are precalculated and entered directly.
!
!    NSIZE() - the size of each cluster (group) of animals at the
!          moment of detection, listed in the data input in the same
!          order as the corresponding R() values.
!
!    ANGLE() - the horizontal angle between the direction of a transect
!          and the bearing to a cluster of animals at the moment of
!          detection, and measured in degrees.  These data are
!          required if computations are to be based on perpendicular
!          distance modelling; otherwise they are not needed.  
!          The program accepts negative angles in the input (e.g. 
!          for observation left of a transect line): these are
!          pooled with positive angles during computations.
!
!
!    In addition, a number of temporary variables are used within
!    the program.  These include:
!
!    NOP - the number of parameters potentially varied during
!          the minimization process.
!
!    FUNC - the sum of squares of the differences between
!          observed and expected values;
!
!    MAX - the maximum number of function evaluations to be
!          allowed (arbitrarily set at 750);
!
!    IMV - a parameter used to determine the method of
!          calculating probability used in Subroutine GIVEF.
!          It either:
!
!            = 0  uses the mean value of P(r) in an interval;
!            = 1  uses the median value of P(r) in the interval.
!
!          IMV is also set at 1 within the program if the
!          coefficient becomes negative or if the data are
!          visual observations. (The two approaches produce
!          almost identical results.)
!	  
!     H, HMIN, HMAX, HSTAR, HSTST - sum of squares totals used in the
!          search for a function minimum.
!
!     PSTAR(), PSTST() - temporary function values used in the search
!          for a function minimum.
!
!     Three coefficients used in the search for a minimum value, viz:
!
!          A - a reflection coefficient;
!          B - a contraction coefficient; and
!          C - an expansion coefficient.
!
!     STOPC - a stopping criterion used to initiate a test for
!          convergence on a minimum point.
!
!     NLOOP - convergence is tested for every NLOOP times the
!          process changes the simplex.  After initial
!          convergence, NLOOP further changes are allowed
!          before testing for final convergence.
!
!     MFAIL and MSFAIL - two variables used to count the numbers of
!           times when convergence fails (MFAIL) or when the
!           detectability coefficient (S) cannot be calculated (MSFAIL).
!           They are used to correct the numbers of counts used in
!           computing means and standard deviations of key parameters.
!
!     MTEST - a variable used to flag a successful run through Loop 1410,
!           and assist in the computation of MSFAIL at the end of
!           Subroutine GIVEF.
!
!     NCLASS - the number of distances classes in the range used for the
!           analysis.  It is either equal to DMAX-STT, or KDT-STT,
!           depending on whether or not KDT is greater than 1.
!    
!     NOTIN - a count variable which records the number of occasions
!           in which an observation falls outside a selected data
!           range, i.e. where distance is less than STT or more than KDT.
!
!     NOVTKS - the number of groups overtaking the observer during a
!           line transect.
!
!     NUMGRA - the number of individuals or groups detected ahead of 
!           observer during a transect.
!
!     NGROUPS - the total number of observations within the selected
!           range of distances (= NVALS - NOTIN).
!
!     D2L - the product of the population density and twice the
!           total transect length - a convenient variable in
!           calculations.
!
!     ESTDEN  - the estimated population density ('D').
!
!     P  - the conspicuousness coefficient ('a').
!
!     Q  - for visual data (where KDT=0 or >1), the mean vegetation
!          cover proportion (c) in the habitat between animal and
!          observer.
!        - for auditory data (where KDT=1), the conspicuousness 
!          coefficient ('b').
!
!     S  - a coefficient of detectability, usable in density estimation
!          by direct calculation.
!
!     TCOV - a topographical cover value, i.e. the estimated proportion 
!          of topographical cover in the line of sight between observer 
!          and animal.
!
!     TOPDD - the proportion of the population at a given radial 
!          distance that is unobscured by topographical cover (such as
!          hills and ridges).
!
!     DMAX - the maximum direct-line detection distance ('dmax'), 
!          either submitted to the program as F(4) or calculated.
!
!     ESTDMAX - the value of DMAX (= F(4)) calculated from the radial
!          distance data in the program input if F(4) is set at zero.
!
!     RESAMP_DIST() - a value of either R() or Y() reselected at
!          random from the original data as the key part of the 
!          bootstrapping process.
!
!     RMAX - the maximum horizontal detection distance ('rmax'),
!          calculated from DMAX and THH.
!
!     ERMAX - a computational upper-limit maximum detection distance 
!          used as a range limit in certain cases.
!    
!     ESTJ - an overall movement correction factor (J) for line
!          transects, calculated when transect duration and animal
!          movement rate data are entered, using an approximation for
!          the J(k) range 0 < k < 5.
!
!     OBSW - the overall observer rate of travel along the line
!          transects, calculated from the total distance travelled and
!          the time taken, and expressed in m/min.
!
!     RLMEAN - the mean of the natural logarithms of radial distances.
!
!     RLSD \\u2014 the standard deviation of the natural logarithms of the
!          radial distances in the data set.
!
!     Further explanatory notes are included within the program below.
!
!
!     This program uses double precision for most real numbers.
!
!
      SUBROUTINE calculate_density (params, header, outfile, graph_file)
!
!
      IMPLICIT NONE
!
!     This seqence derived type must be kept in complete agreement with
!     the C-language definition in mnps2.h
!
      TYPE CALC_PARAMS
      SEQUENCE        ! SEQUENCE indicates not to insert internal
                      ! padding for data alignment
      INTEGER nvals, numa, numo, ifx, iry, ns, km
      INTEGER kdt, iprint, jprint, ishow, maxjb
      DOUBLE PRECISION durn, rate, clint, stt, dist, thh, ltmin, ltmax
      DOUBLE PRECISION vgh, pd, ps
      DOUBLE PRECISION f(4), step(4), r(10000)
      INTEGER nsize(10000)
      DOUBLE PRECISION angle(10000)
      INTEGER complete, bootstrap
      DOUBLE PRECISION estden, sden
      END TYPE CALC_PARAMS

      TYPE (CALC_PARAMS) params
!
      CHARACTER*(*) header, outfile, graph_file
!
!     The variables declared below are in rough alphabetical order
!
      INTEGER nsize(10000), nbsz(10000)
      INTEGER bootstrap
      INTEGER i, ia, ic, ie, iflag, ifx, ig, ih, imax
      INTEGER imin, imv, in, ios, iprint, ir, irb, irow, iry, iseed
      INTEGER ishow, j, jbstp, jprint, js, jv, iqsf
      INTEGER k, kdt, km, kprint, loop, max
      INTEGER maxjb, mfail, msfail, mtest, nap, neval, notin, ngroups
      INTEGER nloop, nop, np1, ns, numa, numest, numo, nclass
      INTEGER novtks, numoin, numain, nvals, numgra
      DOUBLE PRECISION a, approx, b, c, cf1dif, cf1sum, cf2dif, cf2sum
      DOUBLE PRECISION cf3dif, cf3sum, clint, coeffnt1, coeffnt2
      DOUBLE PRECISION coeffnt3, dcoeff, dendif, dist, dsum
      DOUBLE PRECISION estden, estj, fnk, frst, func, hmax, hmean, hmin
      DOUBLE PRECISION hstar, hstd, hstst, durn, ltmin, ltmax, obsw
      DOUBLE PRECISION pd, ps, rate, savemn, scf1, scf2, scf3
      DOUBLE PRECISION sden, sns, stopc, stt, test, tcoeff1, tcoeff2
      DOUBLE PRECISION tcoeff3, tcov, tden, thh, vgh, dmax, t001
      REAL rltot, rlmean, rlsum, rlsd, rlf4, rdifsq, estdmax, ttrdenmn
      REAL ttrden, tdsum, tdendif, strden, t,tcl1,tcl2,cl1,cl2,rkdt,s
      REAL r(10000), resamp_dist(10000), y(10000)
      REAL angle(10000), trden(5000)
      DOUBLE PRECISION val(80), valt(80)
      DOUBLE PRECISION g(21,20), step(4), stept(4), f(4), ft(4)
      DOUBLE PRECISION h(21), pbar(20), pstar(20), pstst(20)
      DOUBLE PRECISION coeff1(5000), coeff2(5000), coeff3(5000)
      DOUBLE PRECISION den(5000)
!
!
!     The program accepts up to 10000 data values, each being the total
!     number of observations [N(r) or N(y)] within the class
!     intervals, beginning with that nearest r=0 or y=0 [R(1),
!     NSIZE(1) and ANGLE(1), if supplied].  If ANGLE() is not
!     supplied and IRY==2, computation still proceeds using N(y)
!     values entered in place of N(r) values.
!
!     This program is designed to receive a list-directed data file with
!     a name of the form <filename>.dat and either spaces or commas
!     separating the values.
!
!     The header line should begin with the computer run number each
!     time to avoid confusion, and be bounded by quotation marks.  Other
!     entries must be separated either by commas or spaces (and
!     pre-checked before running).
!
!
!     The program is set to seek a minimum using up to MAX=750
!     iterations, begin with four parameters (NOP) and set a testing
!     point (NLOOP) equal to unity.
!
      max=750
      nop=4
      nloop=1
!
      nvals=params%nvals
      clint=params%clint
      stt=params%stt
      numa=params%numa
      numo=params%numo
      dist=params%dist
      thh=params%thh
      ltmin=params%ltmin
      ltmax=params%ltmax
      ifx=params%ifx
      iry=params%iry
      ns=params%ns
      km=params%km
      kdt=params%kdt
      iprint=params%iprint
      jprint=params%jprint
      ishow=params%ishow
      maxjb=params%maxjb
      vgh=params%vgh
      durn=params%durn
      rate=params%rate
      pd=params%pd
      ps=params%ps
!
!
      Loop_10: DO ig=1,nop
        f(ig)=params%f(ig)
        step(ig)=params%step(ig)
      END DO Loop_10
!
!     Care is required to ensure that the group size data is submitted
!     in PRECISELY the same sequence as the corresponding R(IN) data,
!     and that non-overtaking cases where r=0 are altered to r=0.01.
!
      Loop_20: DO ih=1,nvals
        r(ih)=params%r(ih)
        nsize(ih)=params%nsize(ih)
      END DO Loop_20
!
!     The same requirement applies to data on observing angles.
!
      IF (iry.lt.2) THEN
        Loop_30: DO ih=1,nvals
          angle(ih)=params%angle(ih)
        END DO Loop_30
      END IF
!
!
      OPEN (UNIT=2,FILE=outfile,STATUS='NEW',IOSTAT=ios,ERR=1930)
      novtks=0
!
!
!
!     If no value of the maximum detection distance F(4) has been
!     entered (i.e. f(4)=0), and radial detection distances are 
!     provided in the data, an estimated maximum distance is
!     calculated based on the assumption that the logarithm of
!     the radial detection distance r is distributed according
!     to a normal distribution, with a maximum value at the
!     mean + (t.001 x s.d.), which matches observed values from
!     large data sets.  This is calculated, then a back-transformation
!     undertaken to give an F(4) value.
!
!     If f(4)=0 and only perpendicular distance data are supplied,
!     the logarithms of the calculated perpendicular distances are 
!     assumed to follow a half-normal distribution (a special case of 
!     the folded normal distribution with a mean of zero). Its maximum 
!     value is then estimated as (t.001 x s.d.), where s.d. is the 
!     standard deviation of the half-normal distribution.  If only 
!     perpendicular distances are supplied AND a maximum class boundary
!     distance has been set which is within the observed distribution, 
!     f(4) should NOT be set at 0 but given an approximate (maximum) 
!     value instead.
!
!
      estdmax = 0.0
!
!
      Outer_40: IF (f(4).eq.0) THEN
        rltot = 0.0
        rlmean = 0.0
        rlsum = 0.0
        rdifsq = 0.0
        numgra = 0
        t001 = 3.2704 + 12.9288/nvals
!
        Inner_42: IF (iry.lt.2) THEN
!
!     This option handles all situations with radial data supplied. 
!     The first step is to calculate a logarithmic detection distance 
!     total RLTOT, then a logarithmic mean value RLMEAN, using data
!     from ahead of the observer only.
!
          Loop_40: DO ih=1,nvals
            Inner_40:  IF (r(ih).gt.0) THEN
              rltot = log(r(ih)+1) + rltot
              numgra = numgra + 1
            END IF Inner_40
          END DO Loop_40
!
!
          rlmean = rltot/numgra
!
!     A standard error of the logarithmic r (RLSD) is now calculated,
!     followed by the estimated logarithmic maximum distance RLF4,
!     which is then backtransformed to give an F(4) value.
!
        rlsum = 0.0
!        
          Outer_45: DO ih=1,nvals
            Inner_45: IF (r(ih).gt.0) THEN
              rdifsq = ((log(r(ih)+1) - rlmean)**2.)/(numgra-1)
              rlsum = rlsum + rdifsq
            END IF Inner_45 
          END DO Outer_45
!
!
          rlsd = sqrt(rlsum)
          rlf4 = rlmean + t001*rlsd
!
          f(4)= EXP(rlf4) - 1
!
        ELSE IF (iry.eq.2) THEN
!
!    This option is used if only perp. data are supplied.  The half-
!    normal distribution is first square root transformed, the
!    distance at the final step.
!
          Loop_43: DO ih=1,nvals 
            rdifsq = ((sqrt(r(ih)+1))**2.)/(nvals-1)
            rlsum = rlsum + rdifsq
          END DO Loop_43
!      
          rlsd = sqrt(0.36338*rlsum)
          rlf4 = rlsd*t001
!
          f(4) = (rlf4)**2 - 1
!
        ELSE
!
        END IF Inner_42
!        
        estdmax = f(4)
!
      END IF Outer_40
!
!
!     If a value of either the maximum detection distance F(4) 
!     or the maximum of the selected interval KDT has been
!     entered as more than 80 times the class interval, CLINT is
!     reset at (F(4)-STT or KDT-STT)/80 to avoid computation problems.
!
      Outer_48: IF ((kdt.le.1).and.(f(4).gt.(80*clint))) THEN
        	clint=(f(4)-stt)/80
!
      ELSE IF ((kdt.gt.1).and.(kdt.le.f(4).and.((kdt-stt).gt.(80*clint))
     &)) THEN
            clint=(kdt-stt)/80
!
      ELSE IF ((kdt.gt.1).and.(kdt.gt.f(4)).and.(f(4).gt.(80*clint))) TH
     &EN
            clint=(f(4)-stt)/80
!
      END IF Outer_48
!
!
!     The header line now begins the program output.
!
      WRITE (2,50) header
   50 FORMAT (a)
!
!
!     The program now prints out a statement of the various data
!     supplied as input, as a check on input accuracy.
!
      WRITE (2,51)
   51 FORMAT (/'INPUT:')
!
      IF ((ifx.eq.0) .and. (iry.eq.0)) THEN
        WRITE (2,52)
   52   FORMAT (/,' Line transect data, based on radial distances from o  
     &bserver')
      ELSE IF ((ifx.eq.0) .and. (iry.eq.1)) THEN
        WRITE (2,53)
   53   FORMAT (/,' Line transect data, based on radial distances and ho  
     &rizontal angles')
      ELSE IF ((ifx.eq.0) .and. (iry.eq.2)) THEN
        WRITE (2,54)
   54  FORMAT (/,' Line transect data, using precalculated perpendicular 
     & distances')
      ELSE
        WRITE (2,55)
   55   FORMAT (/' Fixed observing point data')
      END IF
!
      IF ((ifx.eq.0).and.(ns.eq.1)) THEN
        WRITE (2,56)
   56   FORMAT (' Observations from one side of transect line only')
      ELSE IF ((ifx.eq.0).and.(ns.eq.2)) THEN
        WRITE (2,57)
   57   FORMAT (' Observations from both sides of transect line')
      ELSE
      END IF
!
      IF (km.eq.1) THEN
        WRITE (2,58)
   58   FORMAT (' Transect lengths in km')
       ELSE IF (km.eq.2) THEN
        WRITE (2,59)
   59  FORMAT (' Detection distances and transect lengths in km')
       ELSE
      END IF
!
      IF (kdt.gt.1) THEN
        WRITE (2,60) stt,kdt
   60   FORMAT (' Designated distance range: ',f5.1,' -',i5) 
       ELSE
      END IF
!
!     The program prints out the class interval width (CLINT) and
!     either the total transect length (DIST) or the total time
!     spent (DURN) at fixed points, then the class interval set
!     and the total transect length travelled.
!
      IF (ifx.eq.1) THEN
        WRITE (2,65) clint,durn
   65   FORMAT (' Class interval width =',f7.1,
     &  ' m.    Total time spent =',f7.1,' min.')
        GO TO 90
      END IF
!
!
      IF (km.eq.0) THEN
        WRITE (2,70) clint,dist
   70   FORMAT (' Class interval width =',f7.1,
     &  ' m.   Total transect length (L) =',f10.3,' m.')
      ELSE
        WRITE (2,80) clint,dist
   80   FORMAT (' Class interval width =',f7.1,
     &  ' m.   Total transect length (L) =',f10.3,' km.')
      END IF
!
!
!     If transect lengths have been expressed in kilometres,
!     transect length is are converted to metres.
!
      Outer_90: IF (km.gt.0) THEN
        dist=dist*1000
        Inner_90: IF ((km.eq.2).and.(f(4).gt.0)) THEN
            f(4)=1000*f(4)
            estdmax=f(4)
        END IF Inner_90
      END IF Outer_90
!
!
!     The original values of NUMA and NUMO are retained (as NUMOIN
!     and NUMAIN) so they can be printed in the output.  
!
   90 numoin=numo
      numain=numa
!
        WRITE (2,91) durn
   91   FORMAT (' Total time spent =',f7.1,' min.')    
        WRITE (2,92) rate
   92   FORMAT (' Overall population movement rate =',f4.1,' m/min.')
!
      IF (ltmin.lt.999) THEN
        WRITE (2,93) ltmin
   93   FORMAT (' Topography uneven; approximate minimum obscuring dista
     &nce =',f6.1,' m.')
      ELSE
        WRITE (2,94)
   94   FORMAT (' Topography approximately level')
      END IF
!
      WRITE (2,95)
   95 FORMAT (/,'OUTPUT:',/)
!
!
!     If detection distances were entered in kilometres (km=2),
!     then these distances are also converted to metres.
!
      Outer_96: IF (km.eq.2) THEN
        Inner_96: DO ih=1,nvals
          r(ih)=1000*r(ih)
        END DO Inner_96 
      END IF Outer_96       
!
!
!     The program now calculates the mean overall observer movement
!     rate (w) as OBSW=DIST/DURN (DIST in m and DURN in min), provided
!     that an DURN value has been entered and the data are from line
!     transects (IFX=0).  Also, k is defined as k=u/w (FNK=RATE/OBSW).
!     If not, this entire step is bypassed.
!
      Outer_100: IF ((durn.ge.0) .and. (ifx.eq.0)) THEN
        obsw=dist/durn
        fnk=rate/obsw
!
!     Assuming that the actual distance (and not LJ) has been
!     entered, the movement correction factor (J) is calculated
!     using an approximation.  Different approximations are used
!     if k=u/w is less than or greater than 5.
!
        Inner_100: IF (fnk.eq.0) THEN
           estj=1
         ELSE IF (fnk.gt.0 .and. fnk.le.1) THEN
           estj = 1.000 + (0.00495562*fnk) + (0.0995941*fnk**2)
     &  + (0.0324447*fnk**3)
         ELSE IF (fnk.gt.1 .and. fnk.le.5) THEN
           estj = 1.0051 - (0.0212978*fnk) + (0.0002868*fnk**2) 
     &  + (0.279106*fnk**3) - (0.12982*fnk**4)
     &  + (0.0234994*fnk**5) - (0.00153274*fnk**6) 
         ELSE
           estj=0.8183*fnk
        END IF Inner_100
!
!
!     The movement-corrected overall distance travelled (LJ) is
!     calculated, overriding the DIST value submitted originally.
!     DURN and RATE are both set at 0 to avoid later computation problems
!     before the IF . . THEN loop ends.
!
        dist = estj*dist
        durn = 0
        rate = 0
      END IF Outer_100 
!
!
!     F(3) is now raised in value to approximate D2LJ*Ns*Pd in the case
!     of line transect data, or D2ut*Ps*Pd for fixed point data, in
!     order to prepare for comparisons with observed values within
!     Subroutine GIVEF.  
!
        sns=float(ns)
!
      IF (ifx.eq.0) THEN
         f(3) = (sns*dist*pd*f(3))/1.e4
         step(3) = (sns*dist*pd*step(3))/1.e4
        ELSE
         f(3)=(2.*ps*pd*rate*durn*f(3))/1.e4
         step(3)=(2.*ps*pd*rate*durn*step(3))/1.e4
      END IF
!
!
!     If calculations are to be based on perpendicular distances (y)
!     from the transect line, and the data supplied are radial
!     distances and angles, perpendicular distances are now calculated
!     for the data entered initially, this action being
!     prompted by IRY having a value less than 2.  0.001 is added
!     to raise 0 values to distinguish them from overtakes. If perp.
!     distance data as such were supplied (IRY=2), this step is
!     bypassed and the distance data are recognized as N(y).
!     At this step, any angle data supplied as negative numbers
!     (E.g. from the left of a transect line) are also converted to
!     positive and pooled with the remainder.  If r=0 and the data 
!     are line transect (IFX=0), the values are treated as
!     'overtakes', y(in) is made zero, and the number of groups 
!     overtaking (NOVTKS) is counted.
!
      IF (iry.eq.2) GO TO 190
      IF (iry.eq.0) GO TO 200
!
      Loop_150: DO in=1, nvals      
        IF ((r(in).eq.0).and.(ifx.eq.0)) THEN          
          y(in)=0.0
          novtks=novtks+1
         ELSE
          y(in)=ABS(r(in)*sin((angle(in)*3.14159265)/180.))+0.001
        END IF
      END DO Loop_150
!
      GO TO 210
!
!
!     If perp. distance data were entered as r values, they
!     are renamed as y values at this stage, unless r=0 
!     when they are 'overtakes' and omitted.  Negative y values
!     submitted to the program (as negative r value) are
!     converted to positive and pooled with the rest.
!     The number of groups overtaking (NOVTKS) is not counted,
!     because they are actually y values, along the transect.
!
  190 Loop_190: DO in=1,nvals
        y(in)=ABS(r(in))
      END DO Loop_190
!
      GO TO 210
!
!
!     For IFX=0, the number of overtakes (NOVTKS) 
!     is counted.  If IFX=1, they are not counted as overtakes.
!
  200 Loop_200: DO in=1,nvals
        IF ((r(in).eq.0).and.(ifx.eq.0)) THEN
          novtks=novtks+1
        END IF
      END DO Loop_200
!
!     The initial values of 'a', 'b', 'D2L, and 'dmax' are retained
!     as FT and the corresponding steps as STEPT to make possible
!     reruns of calculations.
!
  210 Loop_210: DO ia=1,nop
        ft(ia)=f(ia)
        stept(ia)=step(ia)
      END DO Loop_210
!
!
!     The program now calculates a topographical cover value
!     (TCOV).  If the topography is effectively level (LTMIN=999)
!     or the value of f(4) is less than the LTMIN value supplied, 
!     then TCOV is set at zero and topography has no effect on
!     computation of the model.  Otherwise, if LTMAX and LTMIN
!     values are supplied, a TCOV value is calculated.
!
      ltmax=f(4)
!      
      IF ((ltmin.lt.999) .and. (ltmax.gt.ltmin)) THEN
         tcov = 1-(exp(log(1/float(nvals))/(ltmax-ltmin)))
       ELSE
         tcov = 0.
      END IF
!
!     The stopping criterion (STOPC) is set at a suitable value,
!     based on the type of data supplied: radial, perpendicular
!     distance data to the limit of visibility, or perpendicular
!     distance data to a distance limit (KDT>1). Tested on data.
!
      IF ((iry.eq.1) .or. ((iry.eq.2).and.(kdt.le.1))) THEN
        stopc=0.0001
       ELSE IF ((iry.eq.2).and.(kdt.gt.1)) THEN
        stopc=0.00001
       ELSE
        stopc=0.005
      END IF
!
!     The program now increases the value of STOPC to allow for
!     highly variable data, by multiplying by NUMAIN/NVALS.
!
        stopc=stopc*numain/nvals
!
!     If progress reports are required (IPRINT=1), the program
!     prints a heading for them.
!
      IF (iprint.eq.1) THEN
        WRITE (2,280)
  280 FORMAT (/' Progress Report (every 1 to 3 function evaluations):'/
     &/' Sequence:  Evaln.no.  Min.difnce  Cnsp.cfnt.  Cover prpn.  Dens 
     &ity x attributes   max.dist.',/)
      END IF
!
!     The term 'APPROX' is used to test closeness to zero.
!
        approx=1.e-15
!
!     If no values have been submitted in the input, the program sets 
!     the step variations at: a=1.0, b=0.5, and c= 2.0
!
      IF (abs(a) .lt. approx) THEN
        a=1.0
        b=0.5
        c=2.0
      END IF
!
!     NAP is the number of parameters to be varied (i.e. with STEP
!     not equal to zero.
!
!      
      nap=0
      loop=0
      iflag=0
      kprint=0
      dcoeff=0
!
!
!     NCLASS is the number of distance classes in the selected range.
!     0.49 is added to avoid counting errors due to 'chopping'.
!
      nclass=int(((f(4)-stt)/clint)+0.49)
      IF (nclass.gt.80) nclass=80
!
!
!     If all STEP sizes have been set at zero, and MAXJB has not been
!     set at 1, then computation goes
!     to Label 1470, calculates the set of values resulting from
!     the values of a, b, D2L and dmax supplied, and ends.  Otherwise
!     it uses NAP to indicate the number of submitted parameters
!     to be varied during subsequent iterations.
!
      Loop_350: DO i=1,nop
        Inner_350: IF (abs(step(i)).gt.approx) THEN
          nap=nap+1
        END IF Inner_350
      END DO Loop_350
!
        IF ((nap.ge.1).or.((nap.eq.0).and.(maxjb.eq.1))) THEN
          GO TO 370
        END IF
! 
          kprint=1
          sns=float(ns)
!          
        IF (sns.eq.0.) THEN
          sns=2.
        END IF
!
      IF (km.gt.0) THEN
          estden=(1.e6*f(3))/(sns*dist*pd)
        ELSE IF (ifx.eq.0) THEN
          estden=(1.e4*f(3))/(sns*dist*pd)
        ELSE 
          estden=(1.e4*f(3))/(2.*rate*durn*ps*pd)
      END IF
!     
      GO TO 1470
!
!
!     To enable parameter estimation using bootstrapping, the basic
!     MINIM routine is run a predetermined MAXJB times, beginning
!     with an initial run.  A number of functions are set at zero first.
!
  370 jbstp=0
      mfail=0
      msfail=0
      tden=0.0
      ttrden=0.0
      tcoeff1=0.0
      tcoeff2=0.0
      tcoeff3=0.0
!
!
!     The program now sets IMV=1 as a default value, unless the
!     expected frequency distribution is likely to have a sharp peak,
!     which may happen if the data are radial, NCLASS has a low
!     value (say, <20) and the data range is not truncated (KDT is not
!     greater than 1).  [The possibility that Q is negative is dealt
!     with in Subroutine GIVEF.]
!
      IF ((iry.eq.0).and.(kdt.le.1).and.(nclass.lt.20)) THEN
         imv=0
       ELSE
         imv=1
      END IF
!
!
!     Seed the random number generator
!
      iseed=0
!      
      Loop_370: DO in = 1, nvals
         IF (iseed.gt.150000) iseed = iseed / 100
         iseed = iseed + int((r(in) * nsize(in) * 72493))
      END DO Loop_370
!      
      CALL srandnag (iseed)
!
!
!     Loop 1410 now begins.
!
!
      Loop_1410: DO bootstrap=1, maxjb
         params%bootstrap = bootstrap
!
!     The flag variable MTEST is first set at zero.
!
        mtest=0
!
!     The program now computes the first distribution of the numbers
!     detected within each class interval, based on the detection
!     distances [(R(IN)], the numbers in each group [(NSIZE(IN)] and
!     the class interval (CLINT) preset in the input.
!
!     In subsequent runs through Loop 1410, bootstrapping applies
!     (JBSTP=1) and a bootstrapped distribution is used instead.
!
!     The number detected within each class [VAL(IC)], is the sum of the
!     numbers in each class in which the R(IN) or Y(IN) values fall.
!     Calculating the various VAL(IC) values first requires
!     finding which R(IN) or Y(IN) values fall within the interval
!     concerned, then adding all the NSIZE(IN) values which fall
!     within that class.  This will be done for each class interval
!     in turn, beginning with the calculation of VAL(1) for the
!     nearest class to r=0 or STT or y=0 or STT.  If a minimum value in the
!     range (STT) has been specified, or a maximum value for r or y of KDT 
!     (>1), the program also computes the number of data clusters (NOTIN) 
!     below STT and above KDT and subtracts it from NVALS to give the 
!     correct magnitude of NVALS for use in later calculations.
!
!
!     Either the initial r class totals, VALT(), are calculated...
!
!     An alternative computation works with perpendicular distance
!     values:
!
        IF (iry.gt.0) GO TO 430
!
!
      frst=stt
!
!     Computation goes to Step 480 if bootstrapping has begun (JBSTP=1).
!
        IF (jbstp.eq.1) THEN
          GO TO 480
        END IF
!
        numo=0
        numa=0
        ngroups=0
!
!
      Loop_410: DO ic=1,nclass
        val(ic)=0.0
!
!     Numbers ahead (NUMA), overtaking (NUMO) and groups are totalled.
!     The frequency in the class, VAL(IC), is accumulated too.
!
          outer411: DO ir=1,nvals
!
             initr: IF ( (kdt > 1) .AND. (r(ir) > kdt) )  THEN
                CYCLE outer411
             ELSE IF ( (r(ir) > frst) .AND. (r(ir) <= (frst+clint)) )
     &          THEN
                 numa=numa+nsize(ir)
                 ngroups=ngroups+1
                 val(ic)=val(ic)+nsize(ir)
             ELSE IF ( (ic == 1) .AND. (r(ir) == 0.) ) THEN
                 numo=numo+nsize(ir)
                 ngroups=ngroups+1
             ELSE
             END IF initr
!
          END DO outer411
!
          frst=frst+clint
!
      END DO Loop_410
!
!
!     The frequency distribution of the original data is saved,
!     as VALT(IG).
!
      Loop_420: DO ig=1,nclass
          valt(ig)=val(ig)
      END DO Loop_420
!
!
      IF ((maxjb.ne.1).and.(nap.le.0)) THEN 
          GO TO 1320
        ELSE
      END IF
!
!
        jbstp=1
        GO TO 620
!
!
!     Or the initial y class totals, VALT(), are calculated...
!
!     Absolute values of Y() are used because data from the
!     two sides of a transect line are pooled.
!
!
  430   frst=stt
!
!     Computation goes to Step 480 if bootstrapping has begun (JBSTP=1).
!
        IF (jbstp.eq.1) THEN
          GO TO 480
        END IF
!
        numo=0
        numa=0
        ngroups=0
!
      Loop_460: DO ic=1,nclass
          val(ic)=0.0
!
!     If a group is in the included data set, the number in each
!     class is totalled in a series of passes through the values 
!     supplied.
!
        outer2: DO ir=1,nvals
!
            inity: IF ( (kdt > 1) .AND. (y(ir) > kdt) )  THEN
                CYCLE outer2
            ELSE IF ( (y(ir) > frst) .AND. (y(ir) <= (frst+clint)) )
     &          THEN
                 numa=numa+nsize(ir)
                 ngroups=ngroups+1
                 val(ic)=val(ic)+nsize(ir)
             ELSE IF ( (ic == 1) .AND. (y(ir) == 0.) ) THEN  
                 numo=numo+nsize(ir)
                 ngroups=ngroups+1
             ELSE
             END IF inity
!
         END DO outer2
!
          frst=frst+clint
!
      END DO Loop_460
!
!
!     The frequency distribution of the original data is now saved,
!     as VALT(IG).
!
      Loop_470: DO ig=1,nclass
        valt(ig)=val(ig)
      END DO Loop_470
!
        IF ((maxjb.ne.1).and.(nap.le.0)) THEN
          GO TO 1320
         ELSE
        END IF  
!
!
        jbstp=1
        GO TO 620
!
!
!     To calculate a bootstrapped distribution, values of R(IN) and
!     the corresponding group NSIZE(IN) are to be chosen at random
!     with replacement, based on a randomly-selected value of IN,
!     which ranges between 1 and NVALS, the total number of groups of
!     animals detected.  Loop 510 selects these values.
!
!
  480    IF (iry.gt.0) THEN
            CALL resample (y, nsize, nvals, resamp_dist, nbsz)
         ELSE
            CALL resample (r, nsize, nvals, resamp_dist, nbsz)
         END IF
!
!
!     There should now be a new set of (N=NVALS) R (or Y) and NSIZE
!     values to use in putting together a new frequency distribution.
!
!     The calculation path differs according to whether radial or
!     perpendicular distance data are being handled.
!
        IF (iry.gt.0) GO TO 570
!
!
!     Either: radial distance data are handled....
!
      frst=stt
!
!
        numo=0
        numa=0
        ngroups=0
!
!
       Loop_550: DO ic=1,nclass
         val(ic)=0.0
!
!     The first class for r begins just above the frst value so that
!     0's are not included because they are 'overtakes'.
!
          outer3: DO irb=1,nvals
          bstrpr: IF ( (kdt > 1) .AND. (resamp_dist(irb) > kdt) )  THEN
                CYCLE outer3
             ELSE IF ( (resamp_dist(irb) > frst) .AND. (resamp_dist(irb)
     &           <= (frst+clint)) ) THEN
                 numa=numa+nbsz(irb)
                 ngroups=ngroups+1
                 val(ic)=val(ic)+nbsz(irb)
             ELSE IF ( (ic == 1) .AND. (resamp_dist(irb) == 0.) ) THEN  
                 numo=numo+nbsz(irb)
                 ngroups=ngroups+1
             ELSE
             END IF bstrpr
           END DO outer3
!
          frst=frst+clint
!
       END DO Loop_550
!
!
!     RESAMP_DIST(IRB) values need to be reassigned at this point
!     or some values will be carried into subsequent loop
!     iterations.
!
      Loop_560: DO js=1,nvals
          nbsz(js)=0
      END DO Loop_560
!
        GO TO 620
!
!
!     Or: perpendicular distance data are handled.
!
  570 frst=stt
!
!
        numo=0
        numa=0
        ngroups=0
!
!
      Loop_600: DO ic=1,nclass
          val(ic)=0
!
!     Numbers in each class are totalled for the groups included
!     in the data set.
!
        outer4: DO irb=1,nvals
!
          bstrpy: IF ( (kdt > 1) .AND. (abs((resamp_dist(irb))) > kdt) )
     &          THEN
                CYCLE outer4
             ELSE IF ( (abs(resamp_dist(irb)) > frst) .AND. 
     &         (abs(resamp_dist(irb)) <= (frst+clint)) ) THEN
                 numa=numa+nbsz(irb)
                 ngroups=ngroups+1
                 val(ic)=val(ic)+nbsz(irb)
             ELSE IF ( (ic == 1) .AND. (resamp_dist(irb) == 0.) ) THEN  
                 numo=numo+nbsz(irb)
                 ngroups=ngroups+1
             ELSE
           END IF bstrpy
!
         END DO outer4
!
          frst=frst+clint
!
        END DO Loop_600
!
!
!     RESAMP_DIST(IRB) values need to be reassigned at this point
!     or some values will be carried into subsequent loop
!     iterations.
!
      Loop_610:  DO js=1,nvals
          nbsz(js)=0
      END DO Loop_610
!
!
!     The calculated set of values, VAL(IC), in each class interval
!     is now printed out for the set of data concerned if JPRINT=1.
!
!
  620   IF ((jprint.eq.1).and.(maxjb.ne.1)) THEN
          WRITE (2,640) bootstrap
  640     FORMAT (///'Bootstrap Replicate No. =',i4,'      Individuals p 
     &er class:',/)
          WRITE (2,650) (val(i),i=1,nclass)
  650     FORMAT (10f6.0,/)
        END IF
!
!
!     The initial simplex of program MINIM is now set up.
!
!
      Loop_660: DO i=1,nop
        g(1,i)=f(i)
      END DO Loop_660
!
        irow=2
!
        Loop_690: DO i=1,nop
          Inner_690: IF (abs(step(i)).ge.approx) THEN
!
            Loop_700: DO j=1,nop
              g(irow,j)=f(j)
            END DO Loop_700
!
            g(irow,i)=g(irow,i)+step(i)
            irow=irow+1
          END IF Inner_690
        END DO Loop_690
!        
  690 CONTINUE
!
        np1=nap+1
        neval=0
!
      Loop_730: DO i=1,np1
        DO j=1,nop
           f(j)=g(i,j)
        END DO
!  
        h(i)=func
!
        CALL givef (f, h(i), s, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
     & mtest, iqsf, graph_file)
!      
          neval=neval+1
!
!
!     All points in the initial simplex become output if IPRINT=1.
!
        IF (iprint.eq.1) THEN
          WRITE (2,720) neval,h(i),(f(j),j=1,nop)
  720     FORMAT (/3x,i4,4x,e13.6,8(1x,e13.6)/24x,8(1x,e13.6)/24x,4(1x,
     &     e13.6))
        END IF
      END DO Loop_730
!
!
!     Now follows the basic loop.  That is, given a simplex, it
!     determines a new simplex and tests for convergence as
!     required (following the flow chart given in Nelder and
!     Mead).
!
!     HMAX and HMIN are the maximum and minimum function values
!     of the current simplex.
!
!
  740   loop=loop+1
        imax=1
        imin=1
        hmax=h(1)
        hmin=h(1)
!
!
      Loop_750: DO i=2,np1
!
         Inner_760: IF (h(i).gt.h(imax)) THEN
           imax=i
           hmax=h(i)
         END IF Inner_760
!
         Inner_770: IF (h(i).lt.h(imin)) THEN
           imin=i
           hmin=h(i)
         END IF Inner_770
!
      END DO Loop_750
!
!
!     The centroid of all vertices, excluding the maximum, is
!     now found.
!
      Loop_790: DO i=1,nop
        pbar(i)=0.0
      END DO Loop_790
!
      Loop_800: DO i=1,np1
        Inner_800: IF (i.ne.imax) THEN
          Loop_810: DO j=1,nop
            pbar(j)=pbar(j)+g(i,j)/FLOAT(nap)
          END DO Loop_810
        END IF Inner_800
      END DO Loop_800
!
!
!     The program reflects the maximum through PBAR to PSTAR, and
!     evaluates the function at PSTAR (to give HSTAR).
!
      Loop_820: DO i=1,nop
        pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)
      END DO Loop_820
!
        hstar=0.
!
        CALL givef (pstar, hstar, dcoeff, val, clint, pd, stt,
     &   tcov, thh, vgh, ifx, imv, iry, ishow, kdt, kprint,
     &   dmax, ltmax, ltmin, nclass, numa, numo, 
     &   msfail, maxjb, mtest, iqsf, graph_file)
!
!     The next 5 statements test whether a progress report is
!     required and, if so, provide one.  This procedure occurs
!     frequently in the program.
!
        neval=neval+1
!
!     If the number of function evaluations to date exceeds 750
!     (=MAX), the program prints out parameter values provided that
!     IPRINT and JPRINT have been set at 1.
!
        IF ((neval.gt.max).and.(iprint.eq.1).and.(jprint.eq.1)) THEN
          j=neval/iprint
          k=neval-j*iprint
          IF (k.le.0) WRITE (2,720) neval,hstar,(pstar(j),j=1,nop)
        END IF
!
        IF (hstar.lt.hmin) GO TO 830
        GO TO 900
!
!     If HSTAR is less than HMIN, PBAR is reflected through PSTAR
!     (to give PSTST) and the function is evaluated there (to
!     give HSTST).
!
  830 Loop_830: DO i=1,nop
        pstst(i)=c*(pstar(i)-pbar(i))+pstar(i)
      END DO Loop_830
! 
        hstst=0.
!
        CALL givef (pstst, hstst, dcoeff, val, clint, pd, stt,
     &   tcov, thh, vgh, ifx, imv, iry, ishow, kdt, kprint,
     &   dmax, ltmax, ltmin, nclass, numa, numo, 
     &   msfail, maxjb, mtest, iqsf, graph_file)
!
!     If IPRINT=1 the program prints out the progress of the
!     iteration.  This is not normally required.
!
        neval=neval+1
        IF ((iprint.eq.1).and.(jprint.eq.1)) GO TO 850
        GO TO 870
  850   j=neval/iprint
        k=neval-j*iprint
!
        IF (k.le.0) THEN
        WRITE (2,860) neval,hstst,(pstst(j),j=1,nop)
  860     FORMAT (/3x,i4,4x,e13.6,8(1x,e13.6)/24x,8(1x,e13.6)/24x,4(1x,
     &     e13.6))
        END IF
!
  870   IF (hstst.lt.hmin) GO TO 880
        GO TO 1030
!
!
!     If HSTST is less than HMIN, the maximum point of the current
!     simplex is replaced by PSTST and HMAX is replaced by HSTAR,
!     then a test is performed.
!
  880 Loop_880: DO i=1,nop
        g(imax,i)=pstst(i)
      END DO Loop_880
!
        h(imax)=hstst
        GO TO 1050
!
!
!     If HSTAR is not less than HMIN, the program tests is HSTAR
!     is greater than the function value at all vertices other
!     than the maximum one.
!
  900 Loop_900: DO i=1,np1
        Inner_900: IF (i.ne.imax) THEN
          IF (hstar.lt.h(i)) GO TO 1030
        END IF Inner_900
      END DO Loop_900
!
!
!     If it is less than at least one of these vertices, the
!     maximum point of the current simplex is replaced by PSTAR and
!     HMAX by HSTAR.  A test is then performed.
!
!     If HSTAR is greater than all the function values excluding
!     the maximum, the program tests if HSTAR is greater than HMAX.
!     If it is not, the maximum point of whichever simplex is
!     now in store (depending on whether HSTAR is greater or less
!     than HMAX) is replaced by PSTAR and HMAX by HSTAR, and the
!     contracted point PSTST and the function value there (HSTST)
!     are calculated.
!
        IF (hstar.gt.hmax) GO TO 930
!
      Loop_920: DO i=1,nop
        g(imax,i)=pstar(i)
      END DO Loop_920
!
        hmax=hstar
        h(imax)=hstar
!
  930 Loop_930: DO i=1,nop
        pstst(i)=b*g(imax,i)+(1.0-b)*pbar(i)
      END DO Loop_930
!
        CALL givef (pstst, hstst, dcoeff, val, clint, pd, stt,
     &   tcov, thh, vgh, ifx, imv, iry, ishow, kdt, kprint,
     &   dmax, ltmax, ltmin, nclass, numa, numo, 
     &   msfail, maxjb, mtest, iqsf, graph_file)
!
        neval=neval+1
!
        IF ((iprint.eq.1).and.(jprint.eq.1)) GO TO 950
        GO TO 970
  950   j=neval/iprint
        k=neval-j*iprint
!
        IF (k.le.0) THEN
        WRITE (2,720) neval,hstst,(pstst(j),j=1,nop)
        END IF
!
  970   IF (hstst.gt.hmax) GO TO 990
!
!
!     If HSTST is less than HMAX, the maximum point is replaced by
!     PSTST and HMAX by HSTST.  A test is then applied.
!
      Loop_980: DO i=1,nop
        g(imax,i)=pstst(i)
      END DO Loop_980
!
        h(imax)=hstst
        GO TO 1050
!
!
!     If HSTST is not less than HMAX, each point in the current
!     simplex is replaced by a point midway between its current
!     position and the position of the minimum point of the
!     current simplex.  The function is evaluated at each new
!     vertex and the test performed.
!
  990 Loop_990: DO i=1,np1,1
        Loop_1000: DO j=1,nop,1
          g(i,j)=(g(i,j)+g(imin,j))/2.0
        END DO Loop_1000
      END DO Loop_990
!
!
!     The conditional loop Inner_1010 is the first of several 
!     introduced to ensure the program uses supplied f values 
!     wherever the stepsize is set at zero (to avoid a processor 
!     error?).
!
      Loop_1020:  DO i=1,np1,1
        Loop_1010:  DO j=1,nop,1
          f(j)=g(i,j)
            Inner_1010: IF (step(i).lt.approx) THEN
              f(i)=ft(i)
            END IF Inner_1010
!
        END DO Loop_1010  
!         
          CALL givef (f, h(i), dcoeff, val, clint, pd, stt,
     &     tcov, thh, vgh, ifx, imv, iry, ishow, kdt,
     &     kprint, dmax, ltmax, ltmin, nclass, numa, numo, 
     &     msfail, maxjb, mtest, iqsf, graph_file)
!
          neval=neval+1
!
        Outer_1020: IF ((iprint.eq.1).and.(jprint.eq.1)) THEN
            j=neval/iprint
            k=neval-j*iprint
          Inner_1020: IF (k.lt.0) THEN
            WRITE (2,720) neval,h(i),(f(j),j=1,nop)
          END IF Inner_1020
        END IF Outer_1020
      END DO Loop_1020
!
        GO TO 1050
!
 1030 Loop_1030: DO i=1,nop
        g(imax,i)=pstar(i)
      END DO Loop_1030
!
        h(imax)=hstar
!
!     If LOOP=NLOOP, tests for convergence begin.  Otherwise
!     computation goes back to the beginning of the basic loop.
!
 1050 IF (loop.ne.nloop) GO TO 740
!
!
!   Tests for Convergence -
!
!     The mean and standard deviation of the function values of the
!     current simplex are now calculated.
!
        hstd=0.0
        hmean=0.0
!
      Loop_1060: DO i=1,np1
        hstd=hstd+h(i)*h(i)
        hmean=hmean+h(i)
      END DO Loop_1060
!
        hmean=hmean/FLOAT(np1)
        hstd=(hstd-FLOAT(np1)*hmean*hmean)/FLOAT(np1)
        hstd=SQRT(ABS(hstd))
!
!
!     The parameter values (F) at the centroid of the current
!     simplex and the function value there (FUNC) are now
!     calculated.
!
      Loop_1080: DO i=1,nop,1
          f(i)=0.0
!          
        Loop_1070: DO j=1,np1,1
          f(i)=f(i)+g(j,i)
            Inner_1070: IF (step(i).lt.approx) THEN
              f(i)=ft(i)
            END IF Inner_1070
        END DO Loop_1070
!        
          f(i)=f(i)/FLOAT(np1)
            Inner_1080: IF (step(i).lt.approx) THEN
              f(i)=ft(i)
            END IF Inner_1080      
!
      END DO Loop_1080
! 
        CALL givef (f, func, dcoeff, val, clint, pd, stt, tcov,
     &    thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     &    ltmax, ltmin, nclass, numa, numo, msfail, maxjb,  
     &    mtest, iqsf, graph_file)
!
         neval=neval+1
!
!     If the number of evaluations has exceeded the value of MAX
!     set (=750), the convergence process is judged not to have
!     succeeded in this case.  If so, this particular run of Loop
!     1410 is assumed to have yielded a 'no result'.  A counting
!     variable (MFAIL) is given a value of 1, to be subtracted from
!     the value of MAXJB later, and the next run through the loop
!     begins.
!
        IF (neval.gt.max) GO TO 1090
        GO TO 1150
!
!     MFAIL is a parameter which counts the number of times a series
!     of iterations failed to converge on a minimum within the set
!     maximum number of iterations.  Parameter values are set at
!     zero so that the failed bootstraps do not contribute to the
!     overall estimates.
!
 1090   mfail=mfail+1
          den(bootstrap)=0
          coeff1(bootstrap)=0
          coeff2(bootstrap)=0
          dcoeff=0
!
        IF ((iprint.eq.1).and.(jprint.eq.1)) GO TO 1100
        GO TO 1380
 1100   WRITE (2,1110) max
 1110   FORMAT (' NUMBER OF FUNCTION EVALUATIONS EXCEEDS ',i4)
        WRITE (2,1120) hstd
 1120   FORMAT (' STANDARD ERROR OF FUNCTION VALUES OF LAST SIMPLEX ',
     &   e13.6)
        WRITE (2,1130) (f(i),i=1,nop)
 1130   FORMAT ('  CENTROID OF LAST SIMPLEX  ',8e13.5,(/28x,8e13.5))
        WRITE (2,1140) func
 1140   FORMAT ('  FUNCTION VALUE AT CENTROID   ',e13.6)
!
        GO TO 1375
!
!
 1150   IF (iprint.eq.1) THEN
          WRITE (2,1155) neval,func,(f(j),j=1,nop)
 1155     FORMAT (/3x,i4,4x,e13.6,8(1x,e13.6)/24x,8(1x,e13.6)/24x,4(1x,
     &     e13.6))
        END IF
!
       IF (hstd.lt.stopc) GO TO 1160
!
!     If the standard deviation calculated above is not less than
!     the criterion set (STOPC), IFLAG and LOOP are set to zero
!     and the basic loop begins again.
!
        iflag=0
        loop=0
        GO TO 740
!
 1160   IF ((iprint.eq.1).and.(jprint.eq.1)) GO TO 1170
        GO TO 1210
 1170   WRITE (2,1180)
 1180   FORMAT (' *'/'  INITIAL EVIDENCE OF CONVERGENCE')
        WRITE (2,1190) (f(i),i=1,nop)
 1190   FORMAT ('  CENTROID OF LAST SIMPLEX  ',8e13.5,(/28x,8e13.5))
        WRITE (2,1200) func
 1200   FORMAT ('  FUNCTION VALUE AT CENTROID   ',e13.6)
!
!     If the standard deviation is less than the stopping
!     criterion, IFLAG is set =0 if there was no evidence of
!     convergence on the last test and =1 if there was evidence
!     of convergence.
!
!     If IFLAG=0, reset IFLAG=1,  the mean of the function
!     values of the current simplex are saved (as SAVEMN), and
!     computation goes back to the beginning of the basic loop.
!
 1210  IF (iflag.le.0) THEN
        iflag=1
        savemn=hmean
        loop=0
        GO TO 740
       END IF
!
!     If IFLAG=1, the program tests if the change in the mean is
!     less than the stopping criterion.  If it is, the process is
!     said to have converged.  If not, IFLAG and LOOP are both set
!     at zero and computation reverts to the start of the
!     basic loop.  STOPC and TEST values tested empirically.
!
        IF (hmean.eq.0) GO TO 1250
        test=savemn/hmean
      IF ((test.gt.0.9999995).and.(test.lt.1.0000005)) GO TO 1250
        iflag=0
        loop=0
        GO TO 740
!
!     If JPRINT=1 the program prints out the results of each successful
!     convergence on a minimum.
!
 1250   IF ((jprint.eq.1).and.(maxjb.ne.1)) THEN
          WRITE (2, 1270) neval
 1270     FORMAT (/,' Process converges on minimum after ', i4,
     &' function evaluations'/)
          WRITE (2, 1280) (f(i), i=1, nop)
 1280     FORMAT (' Minimum at   ',4(1x,e13.6))
          WRITE (2, 1290) func
 1290     FORMAT (/' Minimum function value   ',e13.6)
          WRITE (2, 1300)
 1300     FORMAT (/' End of Search'/1x,13('*'))
        END IF
!
!
!     MTEST is set at 1 to flag that convergence has occurred.  This
!     is carried into Subroutine GIVEF to trigger computation of
!     MSFAIL where computation of S cannot occur.
!
 1320   mtest=1
!
!     Program execution returns to the subroutine to yield final 
!     values of F(1), F(2) and F(3).
!
!
        CALL givef (f, func, dcoeff, val, clint, pd, stt, tcov,
     &    thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     &    ltmax, ltmin, nclass, numa, numo, msfail, maxjb,  
     &    mtest, iqsf, jprint, graph_file)
!
!
        kprint=0
!
!
!     The estimated 'best fit' values of the parameters from the current
!     pass through Loop 1410 are now computed, tabulated and printed out.
!
!     This involves dividing by the same variables used earlier to 
!     compare calculated with observed values in Subroutine GIVEF.
!
!     A density estimate (DEN) is calculated from D2L=F(3) by
!     correcting units to no./ha, and dividing by 2L in the case of
!     line transect data, and by 2Vt in the case of fixed-point data.
!
        IF (km.gt.0) THEN
          den(bootstrap)=(1.e6*f(3))/(sns*dist*pd)
          GO TO 1370
        END IF
!
        IF (ifx.eq.0) THEN
          den(bootstrap)=(1.e4*f(3))/(sns*dist*pd)
         ELSE
          den(bootstrap)=(1.e4*f(3))/(2.*rate*durn*ps*pd)
        END IF
!
!     Values of the other parameter estimates are obtained by redefining
!     the estimates of F(1) and F(2), thus:
!
 1370   coeff1(bootstrap)=f(1)
        coeff2(bootstrap)=f(2)
        coeff3(bootstrap)=dcoeff
!
!      The density estimate is transformed to TRDEN by a logarithmic 
!      transformation, then a running total TTRDEN is calculated.
!
        trden(bootstrap)=log(1 + den(bootstrap))
        ttrden=ttrden+trden(bootstrap)
!
!     Running totals of DEN, COEFF1 and COEFF2 are now made to
!     enable calculation of mean values after the loop ends.
!
 1375   tden=tden+den(bootstrap)
        tcoeff1=tcoeff1+coeff1(bootstrap)
        tcoeff2=tcoeff2+coeff2(bootstrap)
        tcoeff3=tcoeff3+coeff3(bootstrap)
!
!
!
      IF (maxjb.eq.1) THEN
         GO TO 1468
      END IF
!
!
!     Loop 1410 now ends, returning calculations to the start until the
!     maximum preset number of bootstraps (MAXJB) value is reached.
!     LOOP is reset to zero to enable a new series of iterations in the
!     basic loop to begin again, as are G(I,J) values.
!
!
 1380   loop=0
        iflag=0
        kprint=0
!
      Loop_1390: DO i=1,np1
        Loop_1395: DO j=1,nop
          g(i,j)=0.0
        END DO Loop_1395
      END DO Loop_1390
!
!     If the number of detections is less than 80, the conspicuousness
!     coefficient step has been set greater than zero, and the cover
!     proportion step size is not zero, then the program saves
!     the 'a' value calculated from the original data by altering the
!     values of ft(1) and stept(1) in subsequent runs, fixing the
!     'a' value to that calculated from the original data. 
!     NAP is reduced by 1 because there is one variable fewer to alter. 
!
      IF ((bootstrap.eq.1) .and. (nvals.lt.80) .and. (step(1).ne.0)  
     &.and. (step(2).ne.0) .and. (ishow.ne.1)) THEN
        ft(1)=f(1)
        stept(1)=0.0
        nap=nap-1
       ELSE
      END IF
!
!
!     Other F and STEP values are reset to their original values.
!
      Loop_1400: DO ie=1,nop
        f(ie)=ft(ie)
        step(ie)=stept(ie)
      END DO Loop_1400
!
!     The main loop 1410 now ends.
!
 1410 END DO Loop_1410
!
!
!     Following completion of the runs through Loop 1410, the means and
!     standard errors of each of the parameters are now calculated,
!     based on the values of the three arrays DEN(bootstrap),
!     COEFF1(bootstrap) and COEFF2(bootstrap).
!
!     The overall means ESTDEN, COEFFNT1 and COEFFNT2 are calculated
!     first.  NUMEST is the number of parameter estimations made.
!
      numest=maxjb-mfail
!
      IF (numest.eq.0) THEN
         numest=1
      END IF
!
!     Each of the following variables is an estimate of a parameter
!     mean value, beginning with the population mean ESTDEN, then
!     the other parameters - conspicuousness coefficient COEFFNT1,
!     cover proportion or attenuation coefficient COEFFNT2 - and the 
!     detectability coefficient COEFFNT3, in sequence.
!
      estden=tden/numest
      coeffnt1=tcoeff1/numest
      coeffnt2=tcoeff2/numest
!
!     A mean of the transformed values, TTRDENMN, is calculated too.
!
      ttrdenmn=ttrden/numest
!
!     Should MSFAIL be identical to NUMEST, 1 is added to NUMEST to
!     prevent division by zero at the next step.
!
      IF (numest.eq.msfail) numest=numest+1
      coeffnt3=tcoeff3/(numest-msfail)
      IF (numest.eq.(msfail+1)) numest=numest-1
!
!
!     The next step is to calculate the standard errors of each
!     parameter, provided that the number of analyses exceeds 1.
!     Each is the standard deviation of the parameter estimates.
!     If NUMEST is 0 or 1, standard error calculation is bypassed.
!
      IF (numest.le.1.) GO TO 1470
      dsum=0.0
!      
      Loop_1430: DO bootstrap=1, maxjb
        IF ((den(bootstrap)).ne.0.) THEN
           dendif=(den(bootstrap)-estden)**2.
           dsum=dsum+dendif
        END IF
      END DO Loop_1430
! 
      sden = SQRT(dsum/(numest-1))
!
      cf1sum=0.0
!      
      Loop_1440: DO bootstrap=1, maxjb
        IF ((coeff1(bootstrap)).ne.0.) THEN
          cf1dif=(coeff1(bootstrap)-coeffnt1)**2.
          cf1sum=cf1sum+cf1dif
        END IF
      END DO Loop_1440
! 
      scf1 = SQRT(cf1sum/(numest-1))
!
      cf2sum=0.0
!      
      Loop_1450: DO bootstrap=1, maxjb
        IF ((coeff2(bootstrap)).ne.0.) THEN
          cf2dif=(coeff2(bootstrap)-coeffnt2)**2.
          cf2sum=cf2sum+cf2dif
        END IF
      END DO Loop_1450
! 
      scf2=sqrt(cf2sum/(numest-1))
      cf3sum=0.0
!      
      Loop_1460: DO bootstrap=1, maxjb
        IF (coeff3(bootstrap).ne.0) THEN
          cf3dif=(coeff3(bootstrap)-coeffnt3)**2.
          cf3sum=cf3sum+cf3dif
        END IF
      END DO Loop_1460
! 
      IF ((numest-msfail-1).le.0) GO TO 1462      
      scf3=sqrt(cf3sum/(numest-msfail-1))
!
!
!     A standard deviation of the transformed means, STRDEN, is now
!     computed.
!
 1462 tdsum=0.0
! 
      Loop_1465: DO bootstrap = 1, maxjb 
        IF (trden(bootstrap).ne.0) THEN
            tdendif=(trden(bootstrap)-ttrdenmn)**2
            tdsum=tdsum+tdendif
        END IF
      END DO Loop_1465
! 
      strden=sqrt(tdsum/(numest-1))
!
!
!     The p=0.05(2) distribution of t with sample number is now
!     approximated by a function T.
!
      t=(1/(0.093914*maxjb**2.09372))+2.03046
!
!
!     The 95% lower (CL1) and upper (CL2) confidence limits are now
!     calculated from the transformed mean and standard deviation,
!     then back-transformed to the original units.
!
      tcl1 = ttrdenmn - t*strden
      tcl2 = ttrdenmn + t*strden
      cl1 = exp(tcl1) - 1
      cl2 = exp(tcl2) - 1
!
!     If all the STEP values supplied have been set at zero, several
!     parameters need to be given valuues that will correspond to
!     those where at least one parameter was estimated.
!
 1470 IF (nap.le.0) THEN
        numest=1
        coeffnt1=f(1)
        coeffnt2=f(2)
      END IF
!
!
!     The products of the program are now output.
!
!     The number of animal groups in the selected range (NGROUPS) is 
!     first, followed by the estimated topographical cover (TCOV), 
!     then the number of actual parameter estimations made (NUMEST).
!
 1468 WRITE (2,1471) ngroups
 1471 FORMAT (/,' Number of Groups in Distance Range = ',i4)
!
      IF (ifx.eq.0) THEN
        WRITE (2,1472) numain
 1472   FORMAT (' Number of Individuals Detected Ahead = ',i5)
      ELSE
        WRITE (2,1473) numain
 1473   FORMAT (' Number of Individuals Detected =',i5)
      END IF  
!
      IF (ifx.eq.0) THEN
        WRITE (2,1474) numoin
 1474   FORMAT (' Number Overtaking (distance unmeasured) = ',i4)
      END IF
!
 1475 WRITE (2,1478) thh
 1478 FORMAT (' Height Difference from Eyelevel = ',f5.1,' m')
!
      IF (ifx.eq.0) THEN
        WRITE (2,1479) estj
 1479   FORMAT (' Movement Correction Factor (J) = ',f6.3)
 	  END IF
!
      IF (ifx.eq.0) THEN
        WRITE (2,1480) dist
 1480   FORMAT (' Adjusted Transect Length (LJ) = ',f11.3,' m')
      END IF
!
      WRITE (2,1481) tcov
 1481 FORMAT (' Topographical Cover Value = ',f6.4)
!
      IF (estdmax.eq.0) THEN
        WRITE (2,1482) f(4)
 1482   FORMAT (' Maximum Detection Distance (preset) =',f7.1,' m')
      ELSE
        WRITE (2,1483) estdmax
 1483   FORMAT (' Maximum Detection Distance (estimated) =',f7.1,' m')
      END IF
!
       IF (maxjb.eq.1) THEN
          WRITE (2, 1485) neval
 1485     FORMAT (///,' Process converges on minimum after ', i4,
     &' function evaluations')
!
!     Subroutine qsf is now called if MAXJB=1.
!
        CALL qsf (f, func, approx, dist, durn, rate, step, stopc,
     & s, val, clint, stt, tcov, thh, vgh, imv, ishow, kprint, ltmax,
     & ltmin, nclass, numa, numo, msfail, mtest, estden, sden, hstst,
     & pd, ps, ifx, iprint, iry, kdt, km, nap, neval, nop, ns, nvals,
     & np1, g, h, maxjb, iqsf, graph_file)
        GO TO 1920
       END IF
!
      WRITE (2,1490) numest
 1490 FORMAT (' Number of Parameter Estimations =',i4)
!
!
!     Now follow general headings for the results table.
!  
!
      WRITE (2,1500)
 1500 FORMAT (/' Estimated Parameter Values:'/)
      kprint=1
!
      WRITE (2,1520)
 1520 FORMAT (' ',77('x')/' x',4(18x,'x'))
      WRITE (2,1530)
 1530 FORMAT (' x    Parameter',5x,'x',6x,'Value',7x,'x',
     &'  Standard Error  ','x',7x,'Unit',7x,'x')
      WRITE (2,1540)
 1540 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/' x ESTIMATED'
     &,8x,'x',3(18x,'x'))
!
!     Density estimates are printed.
!
      IF (km.eq.0) THEN
 1550 IF (maxjb.gt.1) GO TO 1580
      WRITE (2,1570) estden
 1570 FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,
     &'x (indeterminate)  x  indivs./hectare  x')
      GO TO 1650
 1580 WRITE (2,1590) estden,sden
 1590 FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,'x',4x,f10.3,4x,
     &'x  indivs/hectare  x')
      GO TO 1650
      END IF
!
 1600 IF (maxjb.gt.1) GO TO 1630
      WRITE (2,1620) estden
 1620 FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.2,4x,
     &'x (indeterminate)  x   indivs./sq.km. x')
      GO TO 1650
 1630 WRITE (2,1640) estden,sden
 1640 FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,'x',4x,f10.3,4x,
     &'x   indivs./sq.km. x')
!
!     The conspicuousness coefficient is next printed.
!
 1650 WRITE (2,1660)
 1660 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Conspicuousness  x',3(18x,'x'))
      IF (maxjb.eq.1) THEN  
      WRITE (2,1680) coeffnt1
 1680 FORMAT (' x Coefficient  (a) x',4x,f10.4,4x,
     &'x (indeterminate)  x      metres      x')
      GO TO 1710
      ELSE
      WRITE (2,1700) coeffnt1,scf1
 1700 FORMAT (' x Coefficient  (a) x',4x,f10.3,4x,'x',4x,f10.3,4x,'x',
     &6x,'metres',6x,'x')
      END IF
!
!     The second coefficient is either a cover proportion or a sound
!     attenuation coefficient, decided by the values of KDT.
!
 1710 IF (kdt.eq.1) GO TO 1780
 1720 WRITE (2,1730)
 1730 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Cover      ',6x,'x',3(18x,'x'))
      IF (maxjb.eq.1) THEN  
      WRITE (2,1750) coeffnt2
 1750 FORMAT (' x Proportion   (c)  x',4x,f10.4,4x,
     &'x (indeterminate)  x                  x')
      GO TO 1840
      ELSE
      WRITE (2,1770) coeffnt2,scf2
 1770 FORMAT (' x Proportion  (c)  x',4x,f10.4,4x,'x',4x,f10.4,4x,'x',
     &4x,'         ',5x,'x')
      GO TO 1840
      END IF
!
 1780 WRITE (2,1790)
 1790 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Attenuation',6x,'x',3(18x,'x'))
      IF (maxjb.eq.1) THEN
      WRITE (2,1810) coeffnt2
 1810 FORMAT (' x Coefficient  (b) x',4x,f10.4,4x,
     &'x (indeterminate)  x    per  metre    x')
      GO TO 1840
      ELSE
      WRITE (2,1830) coeffnt2,scf2
 1830 FORMAT (' x Coefficient  (b) x',4x,f10.4,4x,'x',4x,f10.4,4x,'x',
     &4x,'per metre',5x,'x')
      END IF
!
!     The table is now ruled off and some other details printed.
!
!
 1840   WRITE (2,1850)
 1850   FORMAT (' x',4(18x,'x')/' ',77('x')/)
!
!
!     The next two parameters are output only if at least one parameter
!     is being estimated.
!
      IF (nap.ge.1) THEN
        WRITE (2,1865) cl1,cl2
 1865   FORMAT(1x,'95% confidence limits for density estimate:',f10.3,
     &' ',f10.3)
        WRITE (2,1860) coeffnt3,scf3
 1860   FORMAT (1x,'Detectability Coefficient (S) =',f9.2,', SE =',f9.2)
      END IF
!
!
!     Key model estimates are now output if ISHOW was originally set
!     at 1.  The best estimates of F(1), F(2) and F(3) must be entered
!     in the subroutine GIVEF to do this.  The totals in the class
!     intervals are those from the initial data, derived from the
!     saved array variable VALT().  The output is then produced
!     within the subroutine.  The same is done with NUMO and NUMA.
!
!     f(1) and f(2) are left at their original values if all parameters
!     were preset (STEP=0).
!
      IF (nap.le.0) THEN
        GO TO 1870
      END IF
!
      f(1)=coeffnt1
      f(2)=coeffnt2
!      
 1870 Loop_1870: DO jv=1,nclass
        val(jv)=valt(jv)
      END DO Loop_1870
!
        numo=numoin
        numa=numain
!
!    Calculated density values are now multiplied by observing situation
!    variables for a last comparison with observed values within
!    Subroutine GIVEF.
!
 1890 IF (ifx.eq.0) THEN
         IF (km.eq.0) THEN
            f(3)=(estden*dist*pd*sns)/1.e4
         ELSE
            f(3)=(estden*dist*pd*sns)/1.e6
         END IF
      ELSE
 1910    f(3)=(estden*2.*ps*rate*durn*pd)/1.e4
      END IF
!
!
 1920 CALL givef (f, func, dcoeff, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb,  
     & mtest, iqsf, graph_file)
!
!
      params%complete = 1
      params%estden = estden
      params%sden = sden
      CLOSE (unit=2)
      RETURN
!
!
 1930 WRITE (6,1940) outfile,ios
 1940 FORMAT (' Error opening ',a40,' - IOS = ',i6)
      RETURN
!
      END SUBROUTINE calculate_density
!
!     **************************************************************
!
!
!     This version of Subroutine GIVEF will handle ground survey
!     data and also aerial survey data for which there is complete
!     visibility ahead of the observer.  The subroutine
!     handles data divided into up to 80 classes, beginning its
!     computations at radial distances which exceed the maximum
!     recognition distance, and working inwards.
!
!     The subroutine calculates an 'expected' frequency in each
!     class, using parameter values supplied from the main
!     program.  It then compares each with the relevant
!     observed value supplied to the program, computes a sum
!     of squares of the differences, and returns this to the main
!     program as FUNC.
!
!
      SUBROUTINE givef (f, func, s, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb, 
     & mtest, iqsf, graph_file)
!
!     Double precision for all real numbers is set, together with
!     common values, dimensions and symbols for key variables.
!
!
      IMPLICIT NONE
!
!
      INTEGER iermax,ifx,imv,ios,iry,ishow,jj,jl,jrlow,jv,iqsf
      INTEGER kdt,kprint,l10,l20,limit,numa,numo,nclass,mtest
      INTEGER md,msfail,maxjb
      CHARACTER*(*) graph_file
      DOUBLE PRECISION aprexi,auc,cint,clint,d2l,dd,dds,ddsm,dh,dif
      DOUBLE PRECISION difsq,dint,dl,dmax,dnr,dnrl,dnrh,dvg,e,ed
      DOUBLE PRECISION corrn,ermax,expd,expdr,expdy,expdv,func
      DOUBLE PRECISION hcint,htot,ltmin,ltmax,obsd,p,pa,pad,pam,pd
      DOUBLE PRECISION pr,prc,prr,prmax,q,qdd,qdmax,qr
      DOUBLE PRECISION rlow,rmax,rr,s,ss,ssh,ssl,ssmax,stt,tcov
      DOUBLE PRECISION texpd,thh,topdd,topmax,tot,tote,tr,vegdd
      DOUBLE PRECISION vegmax,vgh,vh,visdd,vismax,vl,vlowm,wh
      DOUBLE PRECISION vhowh,wl,wm,wtot,yh,yl,yy,yyh,yyl,zh,zl
      DOUBLE PRECISION val(80), f(4)
      DOUBLE PRECISION rout(80),calcnr(80),obsdnr(80)
      DOUBLE PRECISION yout(80),calcny(80),obsdny(80)
!
!
!     TOT is the progressive value of the sum of squares of the
!     differences between observed and calculated values.  It is
!     initially set at zero.
!
      tot=0.
      wtot=0.
!
!     The program now sets upper and lower limits to the
!     parameter values supplied, usually by effectively giving
!     FUNC (through HTOT) a very high value if a parameter falls
!     outside predetermined limits.
!
      p=f(1)
!
!     HTOT is a variable used to deflect the search away from
!     highly improbable values of the parameters; HTOT=1 000 000.
!
      htot=0.
      dmax=f(4)
!
!
      q=f(2)
!
!     If Q is negative, the program sets IMV=1 and so uses the
!     median value of d in an interval as the basis of
!     computations, in order to avoid logarithms of negative
!     values appearing in Approximation 1 below.
!
        IF (q.lt.0) THEN
          imv=1
        END IF
!
!
!     F(3) is renamed D2L for its run through Subroutine GIVEF.
!
      d2l=f(3)
!
!     HTOT is set to a higher value if D2L<0.001 and MAXJB is not 1
!
      IF ((d2l.lt.0.001).and.(iqsf.ne.1)) THEN
          htot=1.e+6+htot
      END IF
!
!     HTOT is set high if a>400 or <0.01
!
      IF (((p.gt.400).or.(p.lt.0.01)).and.(iqsf.ne.1)) THEN
          htot=1.e+6+htot
      END IF
!
!
!     RMAX is set at zero if DMAX is equal to or less than THH
!
      IF (dmax.le.thh) THEN
           rmax=0.0
          GO TO 110
      END IF 
!
!     RMAX - the maximum horizontal recognition distance - is
!     computed from the direct line distance and height difference.
!
  100 rmax=dsqrt(dmax*dmax-thh*thh)
!
!     A range of 0 - 0.2 is set for the attenuation coefficient
!
  110 IF (((q.ge.0.2).or.(q.lt.0.)).and.(iqsf.ne.1)) THEN
          htot=1.e+6+htot
      END IF
!
!
!
!     The theoretical probability of detecting an individual at
!     the maximum recognition distance [PRMAX=Pr(rmax)] is now
!     calculated.
!
      f(4)=dmax
!
!     PA is the square of the conspicuousness coefficient.
!
!
      pa=p*p
!
!     Visibility and audibility will be affected by topographical
!     features in habitats where the ground is not level. The total
!     probability of visibility at DMAX (VISMAX) will be the product of
!     the probability VEGMAX that the animal is unobscured by vegetation
!     at distance DMAX, and the probability TOPMAX that it is unobscured
!     by topography.  VEGMAX will be a function of dMAX (VEGMAX=(1-Q)**DMAX),
!     while TOPMAX is approximated by the function
!     TOPMAX=(1-T)**(LTMAX-DMAX)=(EXP(ln(0.001)/(LTMAX-LTMIN))**(LTMAX-DMAX),
!     and VISMAX=VEGMAX*TOPMAX.  A first step is to calculate TOPMAX,
!     provided that either LTMIN is not very large or RMAX is currently
!     less than LTMIN.
!
      IF ((ltmin.lt.999.).and.(rmax.gt.ltmin)) GO TO 150
      topmax=1.
      GO TO 160
  150 topmax=(exp(-6.9078/(ltmax-ltmin)))**(dmax-ltmin)
!
!     For observing situations (e.g. aerial survey) where there is
!     'ground' cover for only the first part of the direct-line distance
!     d between animal and observer (indicated by the vegetation height
!     (VGH) exceeding zero), this cover will obscure some animals.  The
!     proportion (VISMAX) visible at DMAX will be a function of d
!     (VISMAX=(1-c)**DMAX). The distance obscured (DVG) will have a
!     value equal to (cover height) x (distance d)/ (observer-animal
!     height difference).
!
  160 IF (kdt .eq. 1) GO TO 230
!
!     If the calculated cover proportion is less than zero, 
!     htot is set to a high value to deflect the search away 
!     from such unreal values.
!
  170 IF ((maxjb.ge.1).and.(q.lt.0.)) THEN
        htot=1.e+6+htot
       ELSE
      END IF
!
!     If visual data are supplied, and there is vegetation close
!     to the ground only, and the animals are within that
!     vegetation, the cover proportion correction applies only
!     to the proportion of the direct-line distance potentially
!     obscured by that ground vegetation.
!
  180 IF (kdt .eq. 1) GO TO 230
  190 IF (vgh.gt.0) THEN
        dvg=vgh*dmax/thh
        IF (dmax.le.dvg) GO TO 210
        vegmax=(1.-q)**dvg
        vismax=vegmax*topmax
        GO TO 220
      END IF
  210 vegmax=(1.-q)**dmax
      vismax=vegmax*topmax
  220 ddsm=dmax*dmax
      prmax=pa*vismax/ddsm
!
!     For visual data collected where there is cover, IMV is put
!     =1 to avoid the possibility of negative values being taken
!     to logarithms later in the program.
!
      imv=1
      GO TO 260
!
  230 qdmax=q*dmax
!
!     To prevent overflow during computations, an upper limit of
!     70 and a lower limit of -68 are set to QDMAX.
!
      IF (qdmax.lt.70.) GO TO 240
      qdmax=70.
      GO TO 250
  240 IF ((-qdmax).lt.68.) GO TO 250
      qdmax=-68.
  250 ssmax=dexp(qdmax)
      ddsm=dmax*dmax
      pam=pa/ddsm
      prmax=pam/ssmax
!
!
!     L10, the number of class intervals used for comparisons, is
!     defined.
!
  260 l10=nclass
      tote=0.
!
!     TR is the highest r value in the current frequency class;
!     because the first class computed is that furthest
!     from the observer or from the transect line - TR is initially
!     set at CLINT*NCLASS + STT, i.e. equal to or just below the
!     maximum recognition distance.
!
      tr=clint*nclass + stt
!
!     To compute 'expected' values within each class, each class
!     is subdivided into subclasses, each of width 800/L10, so that
!     CINT is small and L10 x L20=800.  Each subclass
!     corresponds to an observing arc of width CINT ('delta-r')
!     that sweeps forwards ahead of the observer.  L20 is the
!     number of classes in the inner loop.
!     0.5 is added to remove errors due to 'chopping'.
!
      l20=(800/l10)+0.5
!
!     CINT ('delta-r') is the width of each subclass; it is set
!     at the class width divided by the number of subclasses used.
!
      cint=clint/float(l20)
      hcint=cint/2.
!
!     A function ERMAX is defined to be an exact multiple of CLINT
!
      iermax=int((f(4))/clint)
      ermax=float(iermax)*clint
!
!     In the event that PRMAX > 1., HTOT is set to a high value.
!
      IF (prmax.lt.1.) GO TO 290
      htot=prmax*1.e06+htot
!
!     QR is the proportion of the number of animals originally in
!     a subclass still present after animals detected in arcs
!     further out have been disregarded.  It is initially unity.
!
  290 qr=1.
!
!     JRLOW - a control variable used to control the printout
!     of RLOW - is set initially at zero.
!
      jrlow=0
!
!     EXPD is the cumulative expected number within a class; it is
!     set at zero before computations begin; so is EXPDV.
!
      expd=0.
      expdv=0.
      s=0.
!
!     The program computes a correction, CORRN=(NUMA+NUMO)/NUMA, to
!     cater for data sets that contain overtakes.
!
      corrn=float(numa+numo)/float(numa)
!
!
!
!     Loop 1180, which calculates expected values and compares them
!     with observed values, now begins .....
!
!
      Loop_1180:  DO jj=1,l10
!
!     MD is the hth class in the series, calculated in the reverse
!     order to JJ (going from MD=L10 downwards) because computations 
!     begin at RMAX and continue at progressively decreasing r values.
!
        md=nclass-jj+1
!
!     OBSD is the observed value for the hth arc, as supplied
!     in the input to the program.
!     The observed frequency value in the hth class is corrected to
!     to allow for animals that overtake the observer from behind.
!
        obsd=corrn*val(md)
!
!
!     Where the data are radial distance or fixed point data, EXPD
!     accumulates the expected frequency values within a class; it
!     is set at zero before computations for the class begin.
!
!     Where the data are perpendicular distance data (IRY>0), TEXPD
!     accumulates the expected frequency values within each class;
!     it is set at zero before computations for the class begin,
!     while EXPD is allowed to accumulate.
!
        IF (iry.gt.0) GO TO 310
  300   expd=0.
        GO TO 320
  310   texpd=0.
  320   IF (tote.lt.0.) GO TO 330
        GO TO 340
  330   tote=0.
!
!     WL is the lower boundary of the current class.
!
  340   wl=tr-clint
!
!     If the entire class lies above both the highest calculated
!     and highest preset recognition distances, the expected
!     value (EXPDV) is set at zero and the calculations in Loop
!     870 are bypassed.
!
        IF ((wl.gt.rmax).or.(wl.gt.ermax)) GO TO 350
        GO TO 360
  350   expdv=0.
        GO TO 900
!
!     Two alternative ways of calculating the probability of
!     detection are provided.  If the method is based on the area
!     under the probability curve (IMV=0), the area under
!     the curve from y=DH to y=(infinity) is calculated and
!     expressed as YYH.  If IMV=1, this calculation is bypassed.
!
  360   IF (imv.eq.1) GO TO 450
!
!     This method of computation uses one of two approximations to
!     the exponential integral (APREXI).  If bd < 1, Approximation
!     2 is used; if bd is equal to or greater than 1, Approximation
!     1 is used.
!
        zh=q*dh
!
!     To avoid overflow during computations, ZH is set at 68 if
!     ZH > 68.
!
        IF (zh.lt.68.) GO TO 380
        zh=68.
        GO TO 390
!
!     ZH is set at -70 if ZH < 70.
!
  380   IF ((-zh).lt.70.) GO TO 390
        zh=-70.
!
!     If ZH is zero, VHOWH (=VH/WH) is set at 0, and a few lines of
!     calculations are bypassed.
!
  390   IF (zh.eq.0) THEN
          vhowh=0.
          GO TO 430
        END IF  
!
!     The choice of approximation depends on the value of bd.
!
        IF (zh.ge.1.0) GO TO 420
!
!     Approximation 2 is used if bd<1.
!
        aprexi=-.577216+.999992*zh-.249910*zh**2+.055200*zh**3-.009760*
     &   zh**4+.0010792*zh**5-log(zh)
        GO TO 430
!
!     Approximation 1 is used if bd=1 or bd>1.
!
  420   vh=zh*zh+2.334733*zh+0.250621
        wh=zh*zh+3.330457*zh+1.681534
        vhowh=vh/wh
  430   ssh=dexp(zh)
        yh=pa/(dh*ssh)
!
!     YYH is the area under the detectability curve from
!     the current DH value to infinity.
!
        IF (zh.gt.0..and.zh.lt.1.0) GO TO 440
        yyh=yh*(1.-vhowh)
        GO TO 450
  440   yyh=yh-pa*q*aprexi
!
!
!
!     Loop 870, which calculates the expected number in each
!     arc within the class, and adds them together to produce
!     a progressive total (EXPD and TEXPD), now begins .....
!
!
  450 Loop_870: DO jl=1,l20
!
!     DNR is the central r value in an observing arc.
!     It is initially set half of CINT below TR
!
        IF (jl.eq.1) THEN
          dnr=tr-hcint
!
!     DNRH is the highest r value in an observing arc.  It is
!     initially equal to the original TR.
!
          dnrh=tr
        END IF
!
!     DH is the direct-line distance to the outer edge of the
!     observing arc.
!
          dh=dsqrt(thh*thh+dnrh*dnrh)
!
!     DNRL is the lowest r value in each arc; it is calculated
!     simply by subtracting the arc width from the highest r
!     value in the arc.
!
          dnrl=dnrh-cint
!
!     If, in the last class evaluated, DNRL comes to a negative
!     value, it is then set at zero.
!
        IF (dnrl.ge.0) GO TO 470
          dnrl=0.
!
!     DL is the direct-line distance from the observer to the
!     inner edge of the current observing arc.
!
  470     dl=dsqrt(thh*thh+dnrl*dnrl)
!
!     DINT is the difference between DL and DH, and thus the width
!     of the area under the probability density curve between DH
!     and DL.
!
          dint=dh-dl
!
!     If IMV has been set at zero, calculation of the probability
!     of detection now follows, using the same approximation to the
!     exponential integral as previously.  If IMV has been set at
!     1, computation moves to address 570.
!
          IF (imv.eq.1) GO TO 570
          zl=q*dl
          IF (zl.lt.68.) GO TO 490
          zl=68.
          GO TO 500
  490     IF ((-zl).lt.70.) GO TO 500
          zl=-70.
  500     IF (zl.ne.0) GO TO 520
          vlowm=0.
          GO TO 540
!
!     The choice of approximation depends on the value of bd.
!
  520     IF (zl.ge.1.0) GO TO 530
!
!     Approximation 2 is used if bd<1.
!
          aprexi=-.577216+.999992*zl-.249910*zl**2+.055200*zl**3-
     &     .009760*zl**4+.0010792*zl**5-log(zl)
          GO TO 540
!
!     Approximation 1 is used if bd=1 or bd>1.
!
  530     vl=zl*zl+2.334733*zl+0.250621
          wm=zl*zl+3.330457*zl+1.681534
          vlowm=vl/wm
  540     ssl=dexp(zl)
          yl=pa/(dl*ssl)
!
!     YYL is the area under the detectability curve from
!     the current DL value to infinity.
!
          IF (zl.gt.0..and.zl.lt.1.0) GO TO 550
          yyl=yl*(1.-vlowm)
          GO TO 560
  550     yyl=yl-pa*q*aprexi
!
!     AUC is the area under the detectability curve between
!     d=DL and d=DH.
!
  560     auc=abs(yyh-yyl)
!
!     The mean corrected probability [PR=P(r)] of detecting an
!     individual animal present between d=DL and d=DH is given by
!     dividing the area under the curve by the arc width, then
!     subtracting Pr(rmax).  It is the mean height of the curve.
!
          pr=auc/dint
          GO TO 680
!
!     Where IMV was set at 1, computation of the probability
!     is based on its median value in the arc rather than
!     the mean height of its density curve.
!
!     DD is the direct-line distance to the centre of the arc.
!
  570     dd=dsqrt(thh*thh+dnr*dnr)
!
!     Visibility and audibility will be affected by topographical
!     features in habitats where the ground is not level. The total
!     probability of visibility at DD (VISDD) will be the product of 
!     the probability VEGDD that the animal is unobscured by 
!     vegetation at distance DD, and the probability TOPDD that it 
!     is unobscured by topography there.  VEGDD will be a function of 
!     d (VISDD=(1-Q)**DD), while TOPDD is approximated by the function
!     TOPDD=(1-TCOV)**(DD-LTMIN)=(EXP(ln(0.001)/(LTMAX-LTMIN)))**(DD-LTMIN).
!     and VISDD=VEGDD*TOPDD.  A first step is to calculate TOPDD,
!     provided that LTMIN is not a very large value or d is
!     currently less than LTMIN.  [If LTMIN is set at 999 or higher,
!     topography is assumed not to affect detectability, so TOPDD is set
!     at 1 and TCOV is zero.  If DD is less than LTMIN, topography is
!     also assumed not to affect detectability, so TOPDD is again set at
!     1.] 
!
          IF (ltmin.lt.999 .and. dd.ge.ltmin) THEN
             topdd = (1 - tcov)**(dd-ltmin)
          ELSE
             topdd = 1.
          END IF
!
!     For observing situations where there is 'ground'
!     cover for only the first part of the direct-line distance
!     d between animal and observer (indicated by the vegetation
!     height (VGH) exceeding zero), this cover will obscure
!     some animals.  The proportion visible at d (VISDD) will
!     be a function of d (VISDD=(1-Q)**DD). The distance obscured
!     (DVG) will have a value equal to (cover height) x (distance d)/
!     (observer-animal height difference).
!
  600     IF (kdt .eq. 1) GO TO 650
!
  610     IF (vgh.gt.0) THEN
            dvg=vgh*dd/thh
            IF (dd.le.dvg) GO TO 630
            vegdd=(1.-q)**dvg
            visdd=vegdd*topdd
            GO TO 640
          END IF
!
  630     vegdd=(1.-q)**dd
          visdd=vegdd*topdd
  640     dds=dd*dd
          pr=pa*visdd/dds
          GO TO 680
!
  650     qdd=q*dd
          IF (qdd.lt.68.) GO TO 660
          qdd=68.
          GO TO 670
  660     IF ((-qdd).lt.70.) GO TO 670
          qdd=-70.
  670     ss=dexp(qdd)
          dds=dd*dd
          pad=pa/dds
          pr=pad/ss
!
!     PRC [g(r)=P(r)-P(rmax)] is the height of the detectability curve at DD,
!     corrected by subtracting PRMAX.
!
  680     prc=pr-prmax
!
!     Because 0<PRC<1, an upper limit of 1 and a lower limit of
!     0 are set to the probability function.
!
        IF (prc.gt.1) THEN
           prc=1
         ELSE IF (prc.lt.0) THEN
            prc=0
         ELSE
        END IF
!
!     PRR is the product of PRC=g(r) and the arc of width CINT=\u2206r.
!
          prr=prc*cint
!
!     If the median r value in an arc is above both the calculated
!     and preset highest r values, PRR is set at zero and QR at 1
!
          IF (dnr.ge.ermax.and.dnr.ge.rmax) GO TO 710
!
!     If the median r value is below either the calculated 
!     or observed r value, PRR is allowed to
!     have a value other than zero.
!
          IF (dnr.le.rmax.or.dnr.le.ermax) GO TO 720
  710     prr=0.
!
!     E [= \u2206r.g(r).Q(r)] is calculated by multiplying
!     together the probability of detection in the arc [PRC=g(r)]
!     and the probability [Q(r)] that an individual is still
!     available for detection in that arc. TOTE [=\u2211(\u2206r.g(r).Q(r))] is
!     the accumulating probability density function.
!
  720     e=prr*qr
          tote=tote+e
!
!     The subroutine now calculates the number of detections (ED)
!     expected for radial, fixed-point and perpendicular distance data.
!     For radial data, computation is based on the assumption
!     that the arc of radius r and width \u2206r sweeps over 
!     a series of plots of unit area.  The total number of detections 
!     expected in such a plot will be [plot area] x [apparent density]
!     x [probability of detection].
!
!     With perpendicular distance data, ED is the expected number
!     detected at a distance r from the observer in the strips
!     CINT units wide at a perpendicular distance y=r from the
!     transect line.  Because there are one or two such strips, each
!     L times the area of the individual plot at distance r,
!     the expected no. = \u2206r x [DxNsxLJ=d2l] x \u2206r x g(r) x Q(r).
!     NS has already been converted to a floating-point number (SNS).
!
      IF ((ifx.eq.0).and.(iry.ge.1))  THEN
          ed=d2l*cint*e
!
!     With radial distance data, each observing arc sweeps out an
!     area L units long and 2r units wide, within which there are
!     2r/CINT x L individual plots at radial distance r.  Thus the
!     expected total number detected  =  [number expected in one
!     plot]x[number of plots]=[DxNsxLJ=d2l] x r x \u2206r x g(r) x Q(r).
!
        ELSE IF ((ifx.eq.0).and.(iry.lt.1)) THEN
          ed=d2l*e*dnrh
!
!     With fixed-point data, the expected total number detected
!     = [2ut(=D2LxPd)] x Ps x r x \u2206r x g(r) x Q(r).
!
        ELSE
          ed=d2l*e*dnrh
!
      END IF
!
!     The total number expected in a class (EXPD) 
!     is obtained by progressively adding together
!     the ED values from each arc with radii between DNRH and DNRL.
!
  790     expd=ed+expd
!
!     In the case of perpendicular distance data, all the EXPD
!     values within a class are added together to give a TEXPD
!     value for the class (while, for perpendicular distance data, EXPD 
!     continues to accumulate from class to class).
!
          IF (iry.lt.1) GO TO 810
  800     texpd=expd+texpd
!
!     In preparation for the next arc in the range, DNRL is
!     now renamed DNRH (because the lowest value in one arc
!     becomes the highest in the next, going inward), while DL
!     similarly becomes DH.
!
  810     dnrh=dnrl
          dh=dl
!
!     Similarly, in the case of radial distance data, YYL
!     becomes YYH (bypassed if IMV=1).
!
          IF (imv.eq.1) GO TO 830
          yyh=yyl
!
!     The DNR value for the next arc is DNR minus the
!     arc width (CINT).
!
  830     dnr=dnr-cint
!
!     Where PRR is negative, QR is kept at 1, because no animals
!     are in fact removed from the sample at such distances from
!     the observer.
!
        IF (prr.lt.0) THEN
          qr=1.
        END IF
!
!     The proportion still present in the next arc will be
!     the proportion present in the present arc less the
!     proportion expected to have been detected already.
!
!
          qr=(1.-prr)*qr
!
!     If QR has fallen to a value of less than 0.001 - one
!     animal in a thousand - the program is set to print out the
!     approximate value of rmin where this happens.
!
        Outer_850: IF (qr.le.0.001) THEN
          Inner_850: IF (jrlow.le.0) THEN
!
!     RLOW is the inner boundary of the arc by which 99.9% of the
!     detections are expected to have been made.
!
            rlow=dnrl-cint
            jrlow=1
          END IF Inner_850
        END IF Outer_850
!
!
!     Loop 870 now ends.
!
      END DO Loop_870
!
!
!     Loop 1180 is completed by determining the expected value
!     (EXPDV) in the class; this is arrived at differently for
!     radial or fixed-point and perpendicular distance data. In
!     both cases, negative EXPD values (EXPDN) are then added in.
!
        IF (iry.gt.0) GO TO 890
  880   expdv=expd
        GO TO 900
  890   expdv=texpd
!
!     Calculated values are printed out once the program has
!     converged on a minimum, and KPRINT has been set at 1.
!
  900   IF ((wl.gt.rmax) .or. ((kdt.gt.1).and.(wl.gt.kdt))) GO TO 1080
!
!     If D2L has a trial value of zero, calculation of S at this stage
!     is bypassed and a record of this non-computation retained
!     as the temporary variable MSFAIL.  This is done only when the
!     program has converged on a minimum (MTEST=1) and in the final
!     pass through Loop 1180 (JJ=L10). 
!
        IF (d2l.gt.0) THEN
          s=(expdv/d2l) + s
        ELSE IF ((mtest.eq.1).and.(jj.eq.l10)) THEN
          msfail=msfail+1
        END IF
!
!     Final values of key internal functions are now output if
!     KPRINT=1; otherwise this step is bypassed.
!
  930   IF (kprint.eq.0) GO TO 1080
!
!     For output purposes only, EXPDV is redefined as EXPDR in
!     the case of radial distance data, and as EXPDY in the case
!     of perpendicular distance data.  Negative values of EXPDV
!     are printed as '0.0' in the output because the 'observations'
!     are unreal.   
!
		IF (iry.gt.0) GO TO 1000
  950   rr=tr-(clint/2.)
        expdr=expdv
        IF (expdr.lt.0) THEN
          expdr=0.
        END IF
        IF (ishow.le.0) GO TO 990
        IF (jj.eq.1) THEN
        WRITE (2,975) 
  975   FORMAT (//" Columns, in order: P(r), Q(r), pdf, E{N(r)}, E{N(y)}
     &"/)
        END IF
        WRITE (2,980) rr,expdr,obsd
  980   FORMAT ('  r=',f8.1,5x,'Calc.N(r)=',f9.2,5x,'Obsd.N(r)=',f9.1)
  990   rout(jj)=rr
        calcnr(jj)=expdr
        obsdnr(jj)=obsd
        GO TO 1050
!
 1000   yy=tr-(clint/2.)
        expdy=expdv
        IF (expdy.lt.0) THEN
          expdy=0.
        END IF
        IF (ishow.le.0) GO TO 1040
        IF (jj.eq.1) THEN
        WRITE (2,1025) 
 1025   FORMAT (/" Columns in order: P(r), Q(r), pdf, E{N(r)}, E{N(y)}"/
     &)
        END IF
        WRITE (2,1030) yy,expdy,obsd
 1030   FORMAT ('  y=',f8.1,5x,'Calc.N(y)=',f9.2,5x,'Obsd.N(y)=',f9.1)
 1040   yout(jj)=yy
        calcny(jj)=expdy
        obsdny(jj)=obsd
!
!
!      If the control variable ISHOW has been set at 1, a variety 
!      of intermediate variables is output.     
!
 1050   IF (ishow.eq.1) THEN
          WRITE (2,1070) prc,qr,tote,expd,expdv
 1070     FORMAT (1x,5(f12.6,3x),/)
        END IF
!
!
!     The difference between each observed (OBSD) and expected
!     value (EXPDV) for the class is now calculated.
!
!
 1080   dif=obsd-expdv
!
!
!      If an upper limit of KDT (>1) has been set for the input
!      data range, the computed difference between observed and
!      computed values is reset at zero.
!      
        IF ((kdt.gt.1).and.(wl.gt.kdt)) dif=0
!
!
!     Each difference between the expected (EXPDV) and observed
!     (OBSD) value is squared and a running sum of squares total
!     (TOT) computed.  The final overall sum (FUNC) - which
!     includes HTOT values where parameter values are inappropriate
!     - is then calculated at Step 1300 and returned to the 
!     main program.
!
 1090     difsq=dif*dif
          tot=tot+difsq
!
!     TR is reduced by CLINT before the subroutine goes to the
!     next class inward.
!
        tr=tr-clint
!
!
!     Loop 1180 now ends.
!
 1180 END DO Loop_1180
!
!
!     Tne next two outputs occur only if the range of distance
!     values includes the 99.9% r value and/or an r=0 value.
!
!     If the search for a minimum is complete, and KPRINT=1,
!     the program now prints out a value for RLOW and a file
!     which tabulates observed and calculated frequencies.
!
      IF (kprint.le.0) THEN
         GO TO 1300
      END IF
!
      IF ((jrlow.gt.0) .and. (ishow.gt.0)) THEN
        WRITE (2,1200) rlow
 1200   FORMAT (/' 99.9% r value (rmin) =',f7.2,' m ')
      END IF
!
!     If the range of distance values includes a g(0) value, the
!     program estimates detectability at g(y)=0.
!
      IF (tr.le.0.1) THEN
        WRITE (2,1205) tote
 1205   FORMAT (' Est.Detectability at g(y=0): ',f4.2)
      END IF
!
      OPEN (UNIT=3,FILE=graph_file,STATUS='NEW',IOSTAT=ios,
     & ERR=1320)
!
      IF (iry.gt.0) GO TO 1250
!      
 1210 WRITE (3,1220)
 1220 FORMAT (3x,'   Midpt.   Calculated     Observed   ')
! 
      Loop_1240: DO jj=1,l10
        jv=l10-jj+1
!        
        IF (kdt.gt.1) THEN
          limit = int((kdt-stt) / clint)
         ELSE
          limit = int((dmax-stt) / clint)
        END IF
!      
      IF (jj.gt.limit) GO TO 1290
      WRITE (3,1230) rout(jv),calcnr(jv),obsdnr(jv)
 1230 FORMAT (5x,f7.1,6x,f7.2,7x,f6.1,4x)
 1240 END DO Loop_1240
! 
      GO TO 1290
!      
 1250 WRITE (3,1260)
 1260 FORMAT (3x,'   Midpt.   Calculated     Observed   ')
! 
      Loop_1280: DO jj=1,l10
        jv=l10-jj+1
!        
        IF (kdt.gt.1) THEN
          limit = int((kdt-stt) / clint)
         ELSE
          limit = int((dmax-stt) / clint)
        END IF
!        
      IF (jj.gt.limit) GO TO 1290
      WRITE (3,1270) yout(jv),calcny(jv),obsdny(jv)
 1270   FORMAT (5x,f7.1,6x,f7.2,7x,f6.1,4x)
      END DO Loop_1280
! 
 1290 CLOSE (unit=3)
!
!     The totalled sum of squares of the differences between
!     observed and calculated values is now defined as FUNC.
!
 1300 func=tot+htot
      GO TO 1350
!
!     The next two lines print out an error message if needed.
!
 1320 WRITE (6,1330) graph_file,ios
 1330 FORMAT (' Error opening ',a40,' - IOS = ',i6)
!
!
 1350 IF (kprint.eq.1) THEN
        WRITE (2,1370) func
 1370  FORMAT (1x,'Final Difference at Minimum = ',f15.6/)
      END IF
!
      RETURN
      END SUBROUTINE givef
!
!
      SUBROUTINE resample (orig_dist, orig_size, nvals,
     &     resamp_dist, resamp_size)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: orig_size(10000), nvals
      INTEGER, INTENT(OUT) :: resamp_size(10000)
      INTEGER :: ix, n
      REAL, INTENT(IN) :: orig_dist(10000) 
      REAL, INTENT(OUT) :: resamp_dist(10000)
      DOUBLE PRECISION :: x
      DO ix = 1, nvals
         CALL randgen (x)
         n = INT(x*nvals + 1)
         resamp_dist(ix) = orig_dist(n)
         resamp_size(ix) = orig_size(n)
      END DO
      RETURN
      END SUBROUTINE resample
!
!
      SUBROUTINE srandnag (iseed)
!
!  This subroutine sets the integer seed to be used with the
!  companion RANDNAG function to the value of ISEED.  A flag is
!  set to indicate that the sequence of pseudo-random numbers
!  for the specified seed should start from the beginning.
!
!
      IMPLICIT NONE
      INTEGER iseed
      DOUBLE PRECISION jseed,ifrst
      COMMON /seed/ jseed,ifrst
!
      jseed=iseed
      ifrst=0
!
      END SUBROUTINE srandnag
!
!
      SUBROUTINE randgen (result)
!
!
!  This function returns a pseudo-random number for each invocation.
!  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
!  standard number generator whose Pascal code appears in the
!  article:
!     Park, Steven K. and Miller, Keith W., "Random Number
!     Generators: Good Ones are Hard to Find", Communications of
!     the ACM, October, 1988.
!  It is used in the form recommended by the Numerical
!  Algorithms Group (NAG) for UNIX systems.
!
!
      IMPLICIT NONE
      INTEGER hvlue, lvlue, mobymp, modlus, momdmp, mplier, nextn, testv
      DOUBLE PRECISION ifrst, jseed, result

      PARAMETER  (mplier=16807,modlus=2147483647,mobymp=127773,
     & momdmp=2836)
!
      COMMON /seed/ jseed,ifrst
      SAVE  nextn
!
      IF (ifrst.eq.0) THEN
        nextn=jseed
        ifrst=1
      END IF
!
      hvlue = nextn / mobymp
      lvlue = MOD(nextn, mobymp)
      testv = mplier*lvlue - momdmp*hvlue
      IF (testv.gt.0) THEN
        nextn = testv
      ELSE
        nextn = testv + modlus
      END IF
      result = REAL(nextn) / REAL(modlus)
!
      END SUBROUTINE randgen
!
!
      BLOCK DATA randbd
      DOUBLE PRECISION jseed,ifrst
      COMMON /seed/ jseed,ifrst
!
      DATA jseed,ifrst/123456789,0/
!
      END BLOCK DATA randbd
!
!
!
!      SUBROUTINE qsf                                                                       
!                                                                        
!      This subroutine enables a surface fitting procedure to be used
!      when the number of model iterations (maxjb) is set at 1.  It  
!      estimates the variances of the parameters by fitting a quadratic 
!      surface in the region of the minimum (for estimation of the 
!      variance-covariance matrix).
! 
! 
      SUBROUTINE qsf (f, func, approx, dist, durn, rate, step, stopc,
     & s, val, clint, stt, tcov, thh, vgh, imv, ishow, kprint, ltmax,
     & ltmin, nclass, numa, numo, msfail, mtest, estden, sden, hstst,
     & pd, ps, ifx, iprint, iry, kdt, km, nap, neval, nop, ns, nvals,
     & np1, g, h, maxjb, iqsf, jprint, graph_file)
! 
!      Double precision is set for most real numbers, together with
!      common values, dimensions and symbols for key variables.
! 
! 
       IMPLICIT NONE
! 
!      
       INTEGER iprint, i, i1, i2, i3, i4, i5, j, j0, j1, k, np1, krun 
       INTEGER nop, nap, nless1, in, jk, jless1, l, kprint, ltmin
       INTEGER ltmax, maxjb, msfail, mtest, nclass, numa, numo, iqsf       
       INTEGER iplus1, klessi, nu, nl, ij, nb, ndf, neval, ns, jplus1
       INTEGER ifx, iry, nvals, kdt, mnpd,  iless1, imv, ishow, km
       INTEGER jprint
       DOUBLE PRECISION func, step(4), stopc, g(21,20), h(21)
       DOUBLE PRECISION pstar(20), aval(20), hstst, bmat(210), ao, dmax
       DOUBLE PRECISION temp, ymin, vc(210), var(4), vra, sa, sb, sd2l
       DOUBLE PRECISION sdmax, sd, vrb, vrd2l, vrdmax, den, dist, simp
       DOUBLE PRECISION approx, pmin(20), pbar(20), pstst(20), pd, t
       DOUBLE PRECISION df, vrs, clint, s, stt, tcov, thh, val, vgh
       DOUBLE PRECISION estden, sden, rate, durn, ps, test, cl1, cl2
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(4) :: f
       CHARACTER*(*) graph_file
! 
! 
        WRITE (2,10)                                                        
   10 FORMAT (//' Fitting of quadratic surface in region of minimum',//) 
! 
!                                                                        
!      The fitting of the quadratic surface follows the procedure
!      outlined by Nelder and Mead exactly and, where possible,
!      the notation in the comments corresponds to theirs also.
! 
      neval=0 
      dmax=f(4)
!                                                                        
!      Further function evaluations are counted in NEVAL.
!                                                                        
!      The final simplex is expanded to overcome rounding errors.
! 
!      KRUN is a control variable which allows an automatic rerun
!      through the program if the criterion TEST is negative during
!      quadratic surface fitting.  The program begins again using
!      the computed F(I) values and values of STEP set at one 
!      tenth the F(I) values.
! 
!      SIMP is set to a higher value than the stopping criterion STOPC
!      and IQSF is set at 1 to indicate that Subroutine QSF has begun.
!
      simp=2*stopc
      iqsf=1
!
      outer0: DO i=1,np1                                                  
   20   test=abs(h(i))-func
!
        inner01: IF (test.le.approx) THEN
!
          inner02: IF (krun.le.0) THEN
            krun=1
            WRITE (2,30) 
   30 FORMAT (' Rerun needed with new initial & step values ')
! 
            inner03: DO in=1,nop
              step(in)=f(in)/10.
            END DO inner03
!
           GO TO 400
         END IF inner02
! 
         WRITE (2,50)
   50    FORMAT (' Rerun data with class intervals altered ')
           GO TO 400
        END IF inner01
! 
        inner04: IF (test.lt.simp) THEN 
!
         inner05: DO j=1,nop
           pstst(j)=1000*(g(i,j)-f(j))+g(i,j) 
         END DO inner05
!
         CALL givef (pstst, h(i), s, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb,  
     & mtest, iqsf, graph_file)
! 
         neval=neval+1                                                     
         go to 20
        END IF inner04
!
        GO TO 60
      END DO outer0  
!
!
   60 ao=h(1)                                                           
!
!
!      The function values Y0(I) are calculated and stored in AVAL.
! 
       DO i=1,nap                                                    
         i1=i+1
!
           DO j=1,nop                                                    
           pstar(j)=(g(1,j)+g(i1,j))/2.0
           END DO
!
         CALL givef (pstar, aval(i), s, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
     & mtest, iqsf, graph_file) 
!
         neval=neval+1                                                     
       END DO 
!
!
!      The matrix B(I,J) is calculated, and the lower diagonal
!      section stored in the vector BMAT.
! 
       outer: DO i=1,nap                                                  
          i1=i-1                                                            
          i2=i+1  
!
         inner1:  DO j=1,i1 
          j1=j+1
!
           inner2: DO k=1,nop                                                    
             pstst(k)=(g(i2,k)+g(j1,k))/2.0
           END DO inner2
!
          CALL givef(pstst, hstst, s, val, clint, pd, stt, tcov,
     & thh, vgh, ifx, imv, iry, ishow, kdt, kprint, dmax,
     & ltmax, ltmin, nclass, numa, numo, msfail, maxjb,
     & mtest, iqsf, graph_file) 
!
          neval=neval+1                                                     
          l=i*(i-1)/2+j                                                     
          bmat(l)=2.0*(hstst+ao-aval(i)-aval(j))                            
         END DO inner1 
!
       END DO outer
! 
       l=0
!
       DO i=1,nap                                                   
         i1=i+1                                                            
         l=l+i                                                             
         bmat(l)=2.0*(h(i1)+ao-2.0*aval(i))                                
       END DO  
!
!                                                                        
!      The vector A(I) is calculated and stored in AVAL.
! 
       DO i=1,nap                                                   
         i1=i+1                                                            
         aval(i)=2.0*aval(i)-(h(i1)+3.0*ao)/2.0  
       END DO
!
!                                                                        
!      The matrix Q is calculated and stored in the matrix G.  
!      Considering the usual orientation of rows and columns, 
!      TRANS(Q) is stored in G.
! 
       DO i=1,nop                                                   
         pmin(i)=g(1,i)
       END DO
! 
       outer2: DO i=1,nap                                                   
         i1=i+1
!
         inner21: DO j=1,nop                                                   
           g(i1,j)=g(i1,j)-g(1,j) 
         END DO inner21
!
       END DO outer2
! 
       outer3: DO i=1,nap                                                   
         i1=i+1
!
         inner31: DO j=1,nop                                                   
           g(i,j)=g(i1,j) 
         END DO inner31
!
       END DO outer3 
!
!                                                                        
!      The matrix B is inverted, using the modified square root
!      method (see SAZONOV, Geodeziya i Aerofotosyemka, No.6, 1962).
! 
       np1=nap                                                           
       nless1=nap-1                                                      
       i1=1   
! 
       outer4: DO i=2,np1                                                   
          iless1=i-1                                                        
          i3=1
!
         inner41: DO j=2,iless1                                                
            i2=i1+j                                                           
            jless1=j-1
!
           inner42: DO j0=1,jless1                                               
             i4=i3+j0                                                          
             i5=i1+j0                                                          
             bmat(i2)=bmat(i2)-bmat(i4)*bmat(i5)
           END DO inner42
!
            i3=i3+j
         END DO inner41
!
         i3=0                                                              
         i5=i1+i 
!
         inner43: DO k=1,iless1                                                
           i2=i1+k                                                           
           i3=i3+k                                                           
           temp=bmat(i2)/bmat(i3)                                            
           bmat(i5)=bmat(i5)-bmat(i2)*temp                                   
           bmat(i2)=temp 
         END DO inner43
!
!        If the matrix B is not positive definite, MNPD is set =-1,
!        the estimation of variances ends, and the parameter
!        estimates are printed out.
! 
         inner44: IF (bmat(i5).le.0) THEN
           WRITE (2,100)                                                         
  100     FORMAT (' Matrix to be inverted not positive definite ')          
           mnpd=-1                                                          
           GO TO 205                                                         
          END IF inner44
!
           i1=i1+i 
       END DO outer4
! 
       i1=1  
!
!
       outer5: DO i=2,np1                                                                                                               
         iless1=i-1
!
         inner51: DO j=1,iless1                                                
           i2=i1+j                                                           
           jplus1=j+1
!
           inner52: DO k=jplus1,iless1                                           
             i3=j+k*(k-1)/2                                                    
             i4=i1+k                                                           
             bmat(i2)=bmat(i2)+bmat(i3)*bmat(i4)
           END DO inner52
!
           bmat(i2)=-bmat(i2)
         END DO inner51
!
         i1=i1+i 
       END DO outer5
! 
       i1=0  
! 
       DO i=1,np1                                                   
         i1=i1+i                                                           
         bmat(i1)=1.0/bmat(i1)
       END DO
! 
       outer6: DO i=1,nless1                                                
         iplus1=i+1
!
         inner61: DO k=iplus1,np1                                              
           i1=k*(k-1)/2                                                      
           i2=i1+k                                                           
           i3=i1+i                                                           
           temp=bmat(i2)*bmat(i3)
           klessi=k-i                                                        
           inner62: DO j=1,klessi                                                
             j0=i+j-1                                                          
             i4=j0*(j0-1)/2+i                                                  
             i5=i1+j0                                                          
             bmat(i4)=bmat(i4)+bmat(i5)*temp
            END DO inner62
!
         bmat(i3)=temp 
         END DO inner61
!
       END DO outer6
!
!
!      (B**-1)*A is calculated , and stored in H                         
! 
       outer7: DO i=1,nap                                                   
         h(i)=0.0
!
         inner71: DO j=1,nap
!
           inner72: IF (j.le.i) THEN                                                  
             ij=i*(i-1)/2+j                                                    
            ELSE 
             ij=j*(j-1)/2+i 
           END IF inner72
!
           h(i)=h(i)+bmat(ij)*aval(j)
         END DO inner71
!
       END DO outer7 
!
!                                                                        
!      The estimated minimum value (YMIN) and its position (PMIN)
!      are calculated and printed out if IPRINT=1.
! 
       ymin=0 
! 
       DO i=1,nap                                                   
         ymin=ymin+h(i)*aval(i) 
       END DO
! 
       ymin=ao-ymin 
! 
       outer8: DO i=1,nop                                                   
         pstst(i)=0
!
         inner81: DO j=1,nap                                                   
           pstst(i)=pstst(i)+h(j)*g(j,i) 
         END DO inner81
!
       END DO outer8
! 
       DO i=1,nop                                                        
         pmin(i)=pmin(i)-pstst(i)
       END DO
! 
       IF (iprint.eq.1) THEN
         WRITE (2,120) ymin, (pmin(i),i=1,nop)                                   
  120 FORMAT(' Minimum of fitted quadratic surface is ',e13.6,' at: ',/
     $8(1x,e13.6)/8(1x,e13.6)/4(1x,e13.6)//)                            
         WRITE (2,130) func, (f(i),i=1,nop )                                   
  130 FORMAT(/'Compare with minimum found by iteration ',e13.6,' at: ',  
     &/8(1x,e13.6)/8(1x,e13.6)/4(1x,e13.6)//)                           
         WRITE (2,140)                                                        
  140 FORMAT(/' If difference is large, information matrix and errors
     &are inaccurate'///)
      END IF
! 
!                                                                       
!      Q*(B**-1)*TRANSQ is calculated, and its lower diagonal
!      section stored in the vector VC.
! 
       outer9: DO i=1,nop 
!
         inner91: DO j=1,nap                                                   
           h(j)=0.0
!
           inner92: DO k=1,nap
!
             inner93: IF (k.le.j) THEN                                                 
               jk=j*(j-1)/2+k 
              ELSE                                                       
               jk=k*(k-1)/2+j
             END IF inner93
!
  150       h(j)=h(j)+bmat(jk)*g(k,i)
           END DO inner92                                                         
         END DO inner91
!
         inner94: DO j=i,nop                                                   
           ij=i*(i-1)/2+j                                                    
           vc(ij)=0.0 
!
           inner95: DO k=1,nap                                                   
             vc(ij)=vc(ij)+h(k)*g(k,j) 
!
           END DO inner95
         END DO inner94                                                         
       END DO outer9                                                          
!
!
!      The diagonal elements of VC are stored in VAR.
! 
       DO i=1,nop                                                   
         j=i*(i+1)/2                                                       
         var(i)=vc(j) 
           IF (step(i).lt.approx) THEN
             var(i)=0
           END IF
       END DO
!         
! 
!      The inverse of the information matrix is printed out.
! 
       IF (iprint.eq.1) THEN
         WRITE (2,160)                                                        
  160   FORMAT (' Inverse of information matrix '/) 
!
         DO i=1,nop                                                    
           nu=i*(i+1)/2                                                      
           nl=i*(i-1)/2+1                                                    
           WRITE (2,170) (vc(j),j=nl,nu)                                        
  170     FORMAT (8(1X,E13.6)/4X,8(1X,E13.6)/8X,4(1X,E13.6))
         END DO
!
       END IF
! 
         WRITE (2,190) neval
  190 FORMAT(/,i2,' additional evaluations have been used',/)
!
       IF (iprint.eq.1) THEN
         WRITE (2,200)
  200 FORMAT(' End  of  quadratic  surface  fitting ')
       END IF
!
!
!      The estimated 'best fit' values of the parameters, together
!      with their estimated standard errors, are now computed,
!      tabulated and printed out.  Also, KPRINT is set at 1.
! 
  205  kprint=1
!
!
!      A density estimate (DEN) is calculated from D2L=F(3) by
!      correcting units to no./ha, and dividing by 2L in the case of
!      line transect data, and by 2Vt in the case of fixed-point data.
! 
       IF ((km.gt.0).and.(ifx.eq.0)) THEN
           den=(1.e6*f(3))/(ns*dist*pd)
         ELSE IF ((km.eq.0).and.(ifx.eq.0)) THEN
           den=(1.e4*f(3))/(ns*dist*pd)
         ELSE
           den=1.e4*f(3)/(2.*rate*durn*ps*pd)
       END IF
!
!
!      If the variance/covariance matrix is not positive definite
!      (MNPD=-1), the number of degrees of freedom (NDF) is set =0
!
  210 IF (mnpd.lt.0) THEN
         ndf=0
         GO TO 220
       END IF
!
!
!      The next series of steps determines the number of degrees of
!      freedom used to compute the variance estimate.  NDF is set
!      at one less than the number of values (NVALS) supplied to
!      the program, less the number of parameters estimated by
!      the program itself. 
! 
  220  ndf=nclass-nap-1
       df=float(ndf)
! 
!      If NDF is zero, VRS is set to an arbitrarily high value.
! 
       IF (df.le.approx) THEN
         df=approx
       END IF
! 
       vrs=hstst/df
!
!
!      The standard errors are computed and the results printed out.
! 
       vra=abs(var(1)*2.*vrs)
       sa=sqrt(vra)
! 
       IF (var(2).gt.1.e+12) THEN
         var(2)=1.e+12
       END IF
! 
       IF (var(4).gt.1.e+12) THEN
         var(4)=1.e+12
       END IF
! 
       vrb=abs(var(2)*2.*vrs)
       sb=sqrt(vrb)
       vrd2l=abs(var(3)*2.*vrs)
       sd2l=sqrt(vrd2l)
       vrdmax=abs(var(4)*2.*vrs)
       sdmax=sqrt(vrdmax)
!
! 
       IF ((km.eq.1).and.(ifx.eq.0)) THEN
           sd=(1.e6*sd2l)/(ns*dist*pd)
         ELSE IF ((km.eq.0).and.(ifx.eq.0)) THEN
           sd=(1.e4*sd2l)/(ns*dist*pd)
         ELSE 
           sd=1.e4*sd2l/(2.*rate*durn*ps*pd)
       END IF
!
!
!     The program now determines confidence limits for DEN.
!     The p=0.05(2) distribution of t with sample number is first
!     approximated by a function T.
!
      t=(1/(0.093914*nclass**2.09372))+2.03046
!
!
!     The 95% lower (CL1) and upper (CL2) confidence limits are now
!     calculated from the mean and standard deviation.
!
      cl1 = den - t*sd
      cl2 = den + t*sd
!
!
!
       WRITE (2,230)
  230 FORMAT(//' Estimated parameter values: '//)
!
      WRITE (2,240)
  240 FORMAT (' ',77('x')/' x',4(18x,'x'))
      WRITE (2,250)
  250 FORMAT (' x    Parameter',5x,'x',6x,'Value',7x,'x',
     &'  Standard Error  ','x',7x,'Unit',7x,'x')
      WRITE (2,260)
  260 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/' x ESTIMATED
     &',7x,'x',3(18x,'x'))
!
      Outer10: IF (km.eq.0) THEN
        Inner101: IF (ndf.eq.0) THEN
          WRITE (2,270) den
  270     FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,
     &'x (indeterminate)  x  indivs./hectare  x')
        ELSE
          WRITE (2,280) den, sd
  280     FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,'x',4x,f10.3,4x,
     &'x  indivs/hectare  x')
        END IF Inner101
!
      ELSE IF (km.eq.1) THEN
        Inner102: IF (ndf.gt.0) THEN
          WRITE (2,290) den, sd
  290     FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.3,4x,'x',4x,f10.3,4x,
     &'x   indivs./sq.km. x')
        ELSE
          WRITE (2,300) den, sd
  300     FORMAT (' x DENSITY  (D)',5x,'x',4x,f10.2,4x,'x',4x,f10.3,4x,
     &'x (indeterminate)  x')
        END IF Inner102
      END IF Outer10
! 
      WRITE (2,310)
  310 FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Conspicuousness  x',3(18x,'x'))
!
!
      IF (ndf.eq.0) THEN
         WRITE (2,320) f(1)
  320   FORMAT (' x Coefficient  (a) x',4x,f10.4,4x,
     &'x (indeterminate)  x      metres      x')
      ELSE
         WRITE (2,330) f(1),sa
  330   FORMAT (' x Coefficient  (a) x',4x,f10.3,4x,'x',4x,f10.3,4x,'x',
     &6x,'metres',6x,'x')
      END IF
!
!
      Outer11: IF (kdt.ne.1) THEN
!
        WRITE (2,340)
  340   FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Cover      ',6x,'x',3(18x,'x'))
!
        Inner111: IF (ndf.le.0) THEN
           WRITE (2,350) f(2)
  350      FORMAT (' x Proportion   (c)  x',4x,f10.4,4x,
     &'x (indeterminate)  x                  x')
        ELSE
           WRITE (2,360) f(2),sb
  360   FORMAT(' x Proportion  (c)  x',4x,f10.4,4x,'x',4x,f10.4,4x,'x',
     &4x,'         ',5x,'x')
        END IF Inner111
!
      ELSE 
!
        WRITE (2,370)
  370   FORMAT (' x',4(18x,'x')/' ',77('x')/' x',4(18x,'x')/
     &' x Attenuation',6x,'x',3(18x,'x'))
!
        Inner112: IF (ndf.le.0) THEN
           WRITE (2,380) f(2)
  380      FORMAT (' x Coefficient  (b) x',4x,f10.4,4x,
     &'x (indeterminate)  x    per  metre    x')
!
        ELSE
!
          WRITE (2,390) f(2),sb
  390   FORMAT (' x Coefficient  (b) x',4x,f10.4,4x,'x',4x,f10.4,4x,'x',
     &4x,'per metre',5x,'x')
         END IF Inner112
!
      END IF Outer11
! 
        WRITE (2,395)
  395   FORMAT (' x',4(18x,'x')/' ',77('x')///)
!
!
!     The next two parameters are output only if at least one parameter
!     is being estimated.
!
      IF (nap.ge.1) THEN
        WRITE (2,396) cl1,cl2
  396   FORMAT(1x,'95% confidence limits for density estimate:',f10.3,
     &' ',f10.3)
      END IF
!
!
      estden=den
      sden=sd
!
! 
  400 CONTINUE
! 
!
      END SUBROUTINE
!

