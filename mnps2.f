!     PROGRAM WildlifeDensity  (File mnps2.f, Version WD_0.2a5)
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
!     obtained by 'bootstrap' resampling of data from the original data set.
!
!
!     The WildlifeDensity model was designed by David Morgan, Department
!     of Zoology, The University of Melbourne, Vic. 3010, Australia.
!
!
!     ******************************************************************
!
!
!     File mnps2.f has been written in Fortran, and consists of a 
!     main subroutine, SUBROUTINE calculate_density, and two subsidiary
!     subroutines: SUBROUTINE GIVEF and SUBROUTINE SRANDNAG.
!
!
!     The main subroutine determines parameter values for a least
!     squares fit between calculated and observed frequencies,
!     using the simplex method.  Its design is derived from
!     Nelder & Mead, The Computer Journal, January, 1965, as used in
!     the 'MINIM' program by D.E.Shaw, Divn. of Mathematical Statistics,
!     CSIRO, Australia, though extensively modified since.
!
!     Subroutine GIVEF calculates the sum of squares by using the
!     mathematical model, and returns this to the main program.
!     Subroutine SRANDNAG supplies a pseudo-random number used in the
!     data resampling procedure in response to a seed value supplied
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
!          in the absence of vegetation cover.  Set at 0 if the
!          topography is approximately level, and calculated from
!          the data if current input formats are used.
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
!            = 1  for auditory data;
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
!                 of the function minimization process if IRPRINT=1.
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
!     NOP - the number of parameters potentially varied during
!          the minimization process.
!
!     FUNC - the sum of squares of the differences between
!          observed and expected values;
!
!     MAX - the maximum number of function evaluations to be
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
!     RLSD \u2014 the standard deviation of the natural logarithms of the
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
      INTEGER nvals, numa, numo, ifx, iry, ns, km, imv
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
      INTEGER i, ia, ib, ic, ie, iflag, ifx, ig, ih, imax
      INTEGER imin, imv, in, ios, iprint, ir, irb, irow, iry, iseed
      INTEGER ishow, j, jbstp, jk, jprint, jr, js, jv, jx
      INTEGER k, kdt, km, kprint, loop, max
      INTEGER maxjb, mfail, msfail, mtest, nap, neval, notin, ngroups
      INTEGER nloop, nop, np1, ns, numa, numest, numo, nvals, nclass
      INTEGER novtks, numoin, numain
      DOUBLE PRECISION a, approx, b, c, cf1dif, cf1sum, cf2dif, cf2sum
      DOUBLE PRECISION cf3dif, cf3sum, clint, coeffnt1, coeffnt2
      DOUBLE PRECISION coeffnt3, dcoeff, dendif, dist, dlim, dsum
      DOUBLE PRECISION estden, estj, fk, frst, func, hmax, hmean, hmin
      DOUBLE PRECISION hstar, hstd, hstst, durn, ltmin, ltmax, obsw
      DOUBLE PRECISION pd, ps, rate, rmax, savemn, scf1, scf2, scf3
      DOUBLE PRECISION sden, sns, stopc, stt, test, tcoeff1, tcoeff2
      DOUBLE PRECISION tcoeff3, tcov, tden, thh, vgh, x
      REAL rltot, rlmean, rlsum, rlsd, rlf4, rdifsq, estdmax
      REAL r(10000), resamp_dist(10000), y(10000)
      REAL angle(10000)
      DOUBLE PRECISION val(80), valt(80)
      DOUBLE PRECISION g(5, 4), step(4), stept(4), f(4), ft(4)
      DOUBLE PRECISION h(4), pbar(4), pstar(4), pstst(4)
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
      DO 10 ig=1,nop
        f(ig)=params%f(ig)
        step(ig)=params%step(ig)
   10 CONTINUE
!
!     Care is required to ensure that the group size data is submitted
!     in PRECISELY the same sequence as the corresponding R(IN) data,
!     and that non-overtaking cases where r=0 are altered to r=0.01.
!
      DO 20 ih=1,nvals
        r(ih)=params%r(ih)
        nsize(ih)=params%nsize(ih)
   20 CONTINUE
!
!     The same requirement applies to data on observing angles.
!
      IF (iry.eq.2) GO TO 40
      DO 30 ih=1,nvals
        angle(ih)=params%angle(ih)
   30 CONTINUE
!
!
   40 OPEN (UNIT=2,FILE=outfile,STATUS='NEW',IOSTAT=ios,ERR=1940)
      novtks=0
!
!
!
!     If no value of the maximum detection distance F(4) has been
!     entered (i.e. F(4)=0), this estimated maximum distance is
!     calculated based on the assumption that the logarithm of
!     the radial detection distance r is distributed according
!     to a normal distribution, with a maximum value at the
!     mean + (2.5 x s.d.), which matches observed values from
!     large data sets.  This is calculated, then a back-transformation
!     undertaken to give an F(4) value.
!
!     The first step is to calculate a logarithmic detection distance 
!     total RLTOT, then a logarithmic mean value RLMEAN.
!
      estdmax = 0.0
!
      IF (f(4).eq.0) THEN
        rltot = 0.0
        rlmean = 0.0
        rlsum = 0.0
        rdifsq = 0.0
!
      DO 45 ih=1,nvals
        rltot = log(r(ih)+1) + rltot
   45 CONTINUE
!
!
        rlmean = rltot/nvals
!
!     A standard error of the logarithmic r (RLSD) is calculated,
!     followed by the estimated logarithmic maximum distance RLF4,
!     which is then backtransformed to give an F(4) value.
!
        rlsum = 0.0
      DO 48 ih=1,nvals
        rdifsq = (log(r(ih)+1) - rlmean)**2.
        rlsum = rlsum + rdifsq
   48 CONTINUE
!
!
        rlsd = sqrt(rlsum/(nvals-1))
        rlf4 = rlmean + 2.5*rlsd
!
        f(4) = EXP(rlf4) - 1
        estdmax = f(4)
      END IF
!
!
!     If a value of the maximum detection distance F(4) has been
!     entered more than 80 times the class interval, CLINT is
!     reset at F(4)/80 to avoid computation problems.
!
      IF (f(4).gt.(80*clint)) clint=(f(4))/80
!
!
!     The header line now begins the program output.
!
      WRITE (2,50) header
   50 FORMAT (a80)
!
!
!     The program prints out the class interval width (CLINT) and
!     either the total transect length (DIST) or the total time
!     spent (DURN) at fixed points, then the class interval set
!     and the total transect length travelled.
!
      IF (ifx.eq.1) THEN
        WRITE (2,60) clint,durn
   60   FORMAT (/' Class Interval Width =',f7.1,
     &  ' m.    Total Time Spent =',f7.1,' min.')
        GO TO 90
      END IF
!
!
      IF (km.eq.0) THEN
        WRITE (2,70) clint,dist
   70   FORMAT (/' Class Interval Width =',f7.1,
     &  ' m.   Total Transect Length (L) =',f10.3,' m.')
      ELSE
        WRITE (2,80) clint,dist
   80   FORMAT (/' Class Interval Width =',f7.1,
     &  ' m.   Total Transect Length (L) =',f10.3,' km.')
      END IF
!
!
!     If transect lengths have been expressed in kilometres,
!     transect length is are converted to metres.
!
      IF (km.gt.0) THEN
        dist=dist*1000
      END IF
!
!
!     The original values of NUMA and NUMO are retained (as NUMOIN
!     and NUMAIN) so they can be printed in the output.  
!
   90 numoin=numo
      numain=numa
        WRITE (2,95) durn
   95   FORMAT (' Total Time Spent =',f7.1,' min.')

!
!
!     If detection distances were entered in kilometres (km=2),
!     then these distances are also converted to metres.
!
      IF (km.eq.2) THEN
        DO ih=1,nvals
          r(ih)=1000*r(ih)
        END DO  
      END IF        
!
!
!     The program now calculates the mean overall observer movement
!     rate (w) as OBSW=DIST/DURN (DIST in m and DURN in min), provided
!     that an DURN value has been entered and the data are from line
!     transects (IFX=0).  Also, k is defined as k=u/w (FK=RATE/OBSW).
!     If not, this entire step is bypassed.
!
      IF ((durn.gt.0) .and. (ifx.eq.0)) THEN
        obsw=dist/durn
        fk=rate/obsw
!
!     Assuming that the actual distance (and not LJ) has been
!     entered, the movement correction factor (J) is calculated
!     using an approximation.  Different approximations are used
!     if k=u/w is less than or greater than 5.
!
   98 IF (fk.eq.0) THEN
         estj=1
        ELSE IF (fk.gt.0 .and. fk.le.5) THEN
         estj = 1.0051 - (0.0212978*fk) + (0.0002868*fk**2) 
     &  + (0.279106*fk**3) - (0.12982*fk**4)
     &  + (0.0234994*fk**5) - (0.00153274*fk**6)        
        ELSE
         estj = 0.8183*fk
        END IF
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
      END IF     
!
!
!     F(3) is now raised in value to approximate D2LJ*NS*Pd in the case
!     of line transect data, or D2ut*PS*Pd for fixed point data.  
!
        sns=float(ns)
      IF (ifx.eq.0) THEN
         f(3) = (sns*dist*pd*f(3))/1.e4
         step(3) = (sns*dist*pd*step(3))/1.e4
        ELSE
         f(3)=(2.*ps*pd*rate*durn*f(3))/1.e4
         step(3)=(2.*ps*pd*rate*durn*step(3))/1.e4
!
      END IF
!
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
  150 IF (iry.eq.2) GO TO 190
      IF (iry.eq.0) GO TO 205
!
      DO 180 in=1, nvals
        IF (r(in).eq.0) THEN
          y(in)=0.0
          novtks=novtks+1
        ELSE
          y(in)=ABS(r(in)*sin((angle(in)*3.14159265)/180.))+0.001
        END IF
  180 CONTINUE
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
  190 DO 200 in=1,nvals
          y(in)=ABS(r(in))
  200 CONTINUE
!
      GO TO 210
!
!
!     For IFX=0, the number of overtakes (NOVTKS) 
!     is counted.  If IFX=1, they are not counted as overtakes.
!
  205 DO in=1,nvals
        IF ((r(in).eq.0) .and. (IFX.eq.0)) THEN
          novtks=novtks+1
        END IF
      END DO
!
!     The initial values of 'a', 'b', 'D2L, and 'dmax' are retained
!     as FT and the corresponding steps as STEPT to make possible
!     reruns of calculations.
!
  210 DO 220 ia=1,nop
        ft(ia)=f(ia)
        stept(ia)=step(ia)
  220 CONTINUE
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
      IF ((ltmin.lt.999) .and. (ltmax.gt.ltmin)) THEN
         tcov = 1-(exp(log(1/float(nvals))/(ltmax-ltmin)))
      ELSE
         tcov = 0.
      END IF
!
!     DMAX is given an upper limit (DLIM) which is 80 times the
!     interval width (CLINT).
!
        dlim=clint*80.
!
!     The stopping criterion (STOPC) is set at a suitable value,
!     based on the type of data supplied: radial, perpendicular
!     distance data to the limit of visibility, or perpendicular
!     distance data to a distance limit (KDT>1). Tested on data.
!
  240 IF ((iry.eq.1) .or. ((iry.eq.2).and.(kdt.le.1))) THEN
        stopc=0.001
      ELSE IF ((iry.eq.2).and.(kdt.gt.1)) THEN
        stopc=0.001
      ELSE
        stopc=0.01
      END IF
!
!
!     If progress reports are required (IPRINT=1), the program
!     prints a heading for them.
!
      IF (iprint.eq.1) THEN
        WRITE (2,280)
  280   FORMAT (/' PROGRESS REPORT EVERY 1 or 3 FUNCTION EVALUATIONS'/
     &/' Eval. No.    OLS Value ', 15x, 'Parameters')
      END IF
!
!     The term 'APPROX' is used to test closeness to zero.
!
        approx=1.e-20
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
!     If all STEP sizes have been set at zero, computation goes
!     to Label 1470, calculates the set of values resulting from
!     the values of a, b, D2L and dmax supplied, and ends.  Otherwise
!     it uses NAP to indicate the number of submitted parameters
!     to be varied during subsequent iterations.
!
        DO i=1,nop
          IF (abs(step(i)).gt.approx)   nap=nap+1
        END DO
!
        IF (nap.ge.1) THEN
          GO TO 370
        END IF
!
! 
          kprint=1
          sns=float(ns)
        IF (sns.eq.0.) THEN
          sns=2.
        END IF
!
      IF (km.gt.0) THEN
          estden=(1.e6*f(3))/(sns*dist)
        ELSE IF (ifx.eq.0) THEN
          estden=(1.e4*f(3))/(sns*dist)
        ELSE 
          estden=(1.e4*f(3))/(2.*rate*durn*ps)
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
      tcoeff1=0.0
      tcoeff2=0.0
      tcoeff3=0.0
!
!
!     NCLASS is the number of distance classes in the selected range.
!     0.49 is added to avoid counting errors due to 'chopping'.
!
      IF (kdt.le.1) THEN
        nclass=int(((f(4)-stt)/clint)+0.49)
      ELSE 
        nclass=int(((kdt-stt)/clint)+0.49)
      END IF
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
!     Loop 1410 now begins.
!
      DO 1410 bootstrap=1, maxjb
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
!     An alternative computation works with perpendicular distance
!     values (see below).
!
!
        IF (iry.gt.0) GO TO 430
!
!
!     Either the initial r class totals, VALT(), are calculated...
!
!
  380   frst=stt
!
!
      DO 410 ic=1,nclass
        val(ic)=0.0
!
!     Computation goes to Step 480 if bootstrapping has begun (JBSTP=1).
!
        IF (jbstp.eq.1) THEN
          GO TO 480
        END IF
!
!     A variable NOTIN counts the groups excluded from the data set.
!     The frequency in the class, VAL(IC), is accumulated too.
!
        notin=0
!
          DO 400 ir=1,nvals
             IF (kdt.gt.1 .and. r(ir).ge.kdt)  THEN
               notin=notin+1
             ELSE IF ((r(ir).ne.0) .and. (r(ir).le.stt)) THEN
               notin=notin+1
             ELSE IF ((r(ir)).gt.frst .and. (r(ir)).le.(frst+clint))
     &        THEN
               val(ic)=val(ic)+nsize(ir)
             END IF
!
  400     CONTINUE
!
          frst=frst+clint
  410 CONTINUE
!
		  ngroups=nvals-notin
!
!     The frequency distribution of the original data is saved,
!     as VALT(IG).
!
        DO 420 ig=1,nclass
          valt(ig)=val(ig)
  420   CONTINUE
!
        IF (nap.le.0) THEN 
          GO TO 1320
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
        DO 460 ic=1,nclass
          val(ic)=0.0
!
!     Computation goes to Step 480 if bootstrapping has begun (JBSTP=1).
!
        IF (jbstp.eq.1) THEN
          GO TO 480
        END IF
!
!     A variable NOTIN counts the groups excluded from the data set.
!     If a group is in the included data set, the number in each
!     class is totalled in a series of passes through the values 
!     supplied.
!
          notin=0
!
           DO 450 ir=1,nvals
             IF (kdt.gt.1 .and. abs(y(ir)).ge.kdt) THEN
               notin=notin+1
             ELSE IF ((y(ir).ne.0) .and. (abs(y(ir)).le.stt)) THEN
               notin=notin+1
             ELSE IF (ABS(y(ir)).gt.frst .and. ABS(y(ir)).le.
     &         (frst+clint)) THEN
               val(ic)=val(ic)+nsize(ir)
             END IF
!
  450     CONTINUE
!
          frst=frst+clint
  460   CONTINUE
		  ngroups=nvals-notin
!
!     The frequency distribution of the original data is now saved,
!     as VALT(IG).
!
        DO 470 ig=1,nclass
          valt(ig)=val(ig)
  470   CONTINUE
!
        IF (nap.le.0) THEN
          GO TO 1320
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
  480     iseed=0
!
      DO 510 jr=1,nvals
!
!     A seed value is first chosen, based on the (sequential) values of
!     R in the data input, the size of the data set provided (NVALS),
!     the stage reached in the main backstrapping loop (bootstrap),
!     and the relevant step (JR) in Loop 510.  This should minimize
!     the likelihood of common seed values being submitted in
!     different random number searches.
!     In this way, no two seed values are likely either to be
!     identical or to vary in a consistent way.  They are multiplied
!     by a relatively large number (72493) to ensure the product
!     is a correspondingly large number.
!
!     An upper limit is placed on ISEED to prevent integer overflow.
!
       IF (iseed.gt.150000) iseed=iseed/100
         iseed = iseed + (int((r(jr) * nsize(jr) * 72493)) / bootstrap)
!
!     A random number between 0 and 1 is now called from the Subroutine
!     SRANDNAG, multiplied by NVALS to give it a value between 0 and
!     NVALS, and converted to integer form, adding 1 to allow for
!     chopping.
!
          CALL srandnag (iseed)
          CALL randgen (x)
          jx = INT(x*nvals + 1)
!
!     What happens now depends on whether radial or perpendicular
!     distance measurements are being used.  !
!
          IF (iry.gt.0) GO TO 500
!
!     Either: radial distance computations are made.
!
!     The randomly-chosen RESAMP_DIST(JR), the JRth bootstrapped value
!     of R.  The corresponding value of NSIZE, NSIZE(JX), is now 
!     redefined as NBSZ(JR), the JRth bootstrapped value of NSIZE. Loop
!     510 then goes back to its beginning to select another R value,
!     and so on until all NVALS selections have been made.
!
  490      resamp_dist(jr)=r(jx)
           nbsz(jr)=nsize(jx) 
!
        GO TO 510
!
!     Or: perpendicular distance computations are made.
!
!     The randomly-chosen bootstrapped value of Y, Y(JX), is now
!     redefined as RESAMP_DIST(JR), the JRth bootstrapped value of Y. 
!     The corresponding value of NSIZE, NSIZE(JX), is now redefined
!     as NBSZ(JR), the JRth bootstrapped value of NSIZE. Loop
!     510 then goes back to its beginning to select another Y value,
!     and so on until all NVALS selections have been made.
!
  500      resamp_dist(jr)=y(jx)
           nbsz(jr)=nsize(jx) 
!
  510 CONTINUE
!
!
!     There should now be a new set of (N=NVALS) R (or Y) and NSIZE
!     values to use in putting together a new frequency distribution.
!
!     The calculation path differs according to whether radial or
!     perpendicular distance data are being handled.
!
!
        IF (iry.gt.0) GO TO 570
!
!
!     Either: radial distance data are handled....
!
  520     frst=stt
       DO 550 ic=1,nclass
         val(ic)=0.0
!
!     The first class for r begins just above the frst value so that
!     0's are not included because they are 'overtakes'; that variation
!     is not needed for y values.
!
!
!     A variable NOTIN counts the groups excluded from the data set.
!
        notin=0
        numo=0
        numa=0
!
          DO 540 irb=1,nvals
             IF (kdt.gt.1 .and. resamp_dist(irb).gt.kdt)  THEN
                  notin=notin+1
             ELSE IF ((resamp_dist(irb).ne.0) .and. (resamp_dist(irb)
     &         .le.stt)) THEN
                   notin=notin+1
             ELSE IF ((resamp_dist(irb).gt.frst)
     &         .and. (resamp_dist(irb).le.(frst+clint))) THEN
                  val(ic)=val(ic)+nbsz(irb)
             END IF
!
!     New values of NUMO and NUMA are now required because a different
!     selection of data values has been made.  This is done by
!     totalling the R()=0 and R()>0 values in the new sample to get
!     new NUMO and NUMA values respectively.
!
            IF (resamp_dist(irb).eq.0) THEN
                 numo=numo+nsize(irb)
            ELSE IF ((resamp_dist(irb).gt.stt) 
     &        .and. ((kdt.le.1) .or. (kdt.gt.1 
     &        .and. resamp_dist(irb).le.kdt)))  THEN
                 numa=numa+nsize(irb)
            END IF
!
  540     CONTINUE
!
          frst=frst+clint
!
  550  CONTINUE
!
          ngroups=nvals-notin
!
!     RESAMP_DIST(IRB) values need to be reassigned at this point
!     or some values will be carried into subsequent loop
!     iterations.
!
        DO 560 js=1,nvals
          resamp_dist(js)=0.0
          nbsz(js)=0
  560   CONTINUE
        GO TO 620
!
!
!     Or: perpendicular distance data are handled.
!
  570    frst=stt
        DO 600 ic=1,nclass
          val(ic)=0
!
!     A variable NOTIN counts the groups excluded from the data set.
!     Numbers in each class are totalled for the groups included
!     in the data set.
!
          notin=0
          numa=0
          numo=0
!
           DO 590 irb=1,nvals
             IF (kdt.gt.1 .and. abs(resamp_dist(irb)).gt.kdt)  THEN
               notin=notin+1
             ELSE IF ((abs(resamp_dist(irb)).ne.0) .and. 
     &         (abs(resamp_dist(irb)).le.stt)) THEN
               notin=notin+1
             ELSE IF ((abs(resamp_dist(irb)).gt.frst)
     &         .and. (abs(resamp_dist(irb)).le.(frst+clint))) THEN               
               val(ic)=val(ic)+nbsz(irb)
           END IF
!
!     New values of NUMO and NUMA are required because a different
!     selection of data values has been made.  This is done by
!     totalling the Y()=0 and Y()>0 values in the new sample to get
!     new NUMO and NUMA values respectively.
!
            IF (resamp_dist(irb).eq.0) THEN
                 numo=numo+nsize(irb)
            ELSE IF ((resamp_dist(irb).gt.stt) 
     &        .and. ((kdt.le.1) .or. (kdt.gt.1 
     &        .and. resamp_dist(irb).le.kdt)))  THEN
                 numa=numa+nsize(irb)
            END IF

!
  590     CONTINUE
!
          frst=frst+clint
!
  600   CONTINUE
!
          ngroups=nvals-notin
!
!     RESAMP_DIST(IRB) values need to be reassigned at this point
!     or some values will be carried into subsequent loop
!     iterations.
!
        DO 610 js=1,nvals
          resamp_dist(js)=0.0
          nbsz(js)=0
  610   CONTINUE
!
!
!     The calculated set of values, VAL(IC), in each class interval
!     is now printed out for the set of data concerned if JPRINT=1.
!
!
  620   IF (jprint.eq.1) THEN
          WRITE (2,640) bootstrap
  640     FORMAT (///'  Bootstrap Replicate No. =',i4/)
          WRITE (2,650) (val(i),i=1,nclass)
  650     FORMAT (8f10.4,//)
        END IF
!
!
!     The initial simplex of program MINIM is now set up.
!
!
  660   DO 670 i=1,nop
  670     g(1,i)=f(i)
        irow=2
        DO 690 i=1,nop
          IF (abs(step(i)).lt.approx) GO TO 690
            DO j=1,nop
              g(irow,j)=f(j)
            END DO
          g(irow,i)=g(irow,i)+step(i)
          irow=irow+1
  690   CONTINUE
        np1=nap+1
        neval=0
        DO 730 i=1,np1
          DO j=1,nop
           f(j)=g(i,j)
          END DO
        CALL givef (f, h(i), dcoeff, val, clint, pd, ps, stt,
     &     tcov, thh, vgh, sns, ifx, imv, iry, ishow, kdt,
     &     kprint, ltmax, ltmin, nclass, numa, numo, nvals, 
     &     ns,msfail, mtest, graph_file)

!
          neval=neval+1
!
!     All points in the initial simplex become output if IPRINT=1.
!
        IF (iprint.eq.1) THEN
          WRITE (2,720) neval,h(i),(f(j),j=1,nop)
  720     FORMAT (/3x,i4,4x,e13.6,8(1x,e13.6)/24x,8(1x,e13.6)/24x,4(1x,
     &     e13.6))
        END IF
  730  CONTINUE
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
!
  740   loop=loop+1
        imax=1
        imin=1
        hmax=h(1)
        hmin=h(1)
        DO 780 i=2,np1
          IF (h(i).gt.h(imax)) GO TO 750
          GO TO 760
  750     imax=i
          hmax=h(i)
  760     IF (h(i).lt.h(imin)) GO TO 770
          GO TO 780
  770     imin=i
          hmin=h(i)
  780   CONTINUE
!
!     The centroid of all vertices, excluding the maximum, is
!     now found.
!
        DO 790 i=1,nop
  790     pbar(i)=0.0
        DO 810 i=1,np1
          IF (i.eq.imax) GO TO 810
          DO 800 j=1,nop
  800       pbar(j)=pbar(j)+g(i,j)/FLOAT(nap)
  810   CONTINUE
!
!     The program reflects the maximum through PBAR to PSTAR, and
!     evaluates the function at PSTAR (to give HSTAR).
!
        DO 820 i=1,nop
  820     pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i)
        CALL givef (pstar, hstar, dcoeff, val, clint, pd, ps, stt,
     &   tcov, thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     &   ltmax, ltmin, nclass, numa, numo, nvals, ns, 
     &   msfail,mtest, graph_file)
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
  830   DO 840 i=1,nop
  840     pstst(i)=c*(pstar(i)-pbar(i))+pstar(i)
        CALL givef (pstst, hstst, dcoeff, val, clint, pd, ps, stt,
     &   tcov, thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     &   ltmax, ltmin, nclass, numa, numo, nvals, ns, 
     &   msfail, mtest, graph_file)
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
        WRITE (2,720) neval,hstst,(pstst(j),j=1,nop)
        END IF
!
  870   IF (hstst.lt.hmin) GO TO 880
        GO TO 1030
!
!     If HSTST is less than HMIN, the maximum point of the current
!     simplex is replaced by PSTST and HMAX is replaced by HSTAR,
!     then a test is performed.
!
  880   DO 890 i=1,nop
  890     g(imax,i)=pstst(i)
        h(imax)=hstst
        GO TO 1050
!
!     If HSTAR is not less than HMIN, the program tests is HSTAR
!     is greater than the function value at all vertices other
!     than the maximum one.
!
  900   DO 910 i=1,np1
          IF (i.eq.imax) GO TO 910
          IF (hstar.lt.h(i)) GO TO 1030
  910   CONTINUE
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
        DO 920 i=1,nop
  920     g(imax,i)=pstar(i)
        hmax=hstar
        h(imax)=hstar
  930   DO 940 i=1,nop
  940     pstst(i)=b*g(imax,i)+(1.0-b)*pbar(i)
        CALL givef (pstst, hstst, dcoeff, val, clint, pd, ps, stt,
     &   tcov, thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     &   ltmax, ltmin, nclass, numa, numo, nvals, ns, 
     &   msfail, mtest, graph_file)
!
        neval=neval+1
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
!     If HSTST is less than HMAX, the maximum point is replaced by
!     PSTST and HMAX by HSTST.  A test is then applied.
!
        DO 980 i=1,nop
  980     g(imax,i)=pstst(i)
        h(imax)=hstst
        GO TO 1050
!
!     If HSTST is not less than HMAX, each point in the current
!     simplex is replaced by a point midway between its current
!     position and the position of the minimum point of the
!     current simplex.  The function is evaluated at each new
!     vertex and the test performed.
!
  990   DO 1000 i=1,np1
          DO 1000 j=1,nop
 1000     g(i,j)=(g(i,j)+g(imin,j))/2.0
        DO 1020 i=1,np1
          DO 1010 j=1,nop
 1010       f(j)=g(i,j)
          CALL givef (f, h(i), dcoeff, val, clint, pd, ps, stt,
     &     tcov, thh, vgh, sns, ifx, imv, iry, ishow, kdt,
     &     kprint, ltmax, ltmin, nclass, numa, numo, nvals, 
     &     ns, msfail, mtest, graph_file)
!
          neval=neval+1
          IF ((iprint.eq.1).and.(jprint.eq.1)) THEN
            j=neval/iprint
            k=neval-j*iprint
!
            IF (k.lt.0) THEN
              WRITE (2,720) neval,h(i),(f(j),j=1,nop)
            END IF
          END IF
!
 1020   CONTINUE

        GO TO 1050
 1030   DO 1040 i=1,nop
 1040     g(imax,i)=pstar(i)
        h(imax)=hstar
!
!     If LOOP=NLOOP, tests for convergence begin.  Otherwise
!     computation goes back to the beginning of the basic loop.
!
 1050   IF (loop.ne.nloop) GO TO 740
!
!
!   Tests for Convergence -
!
!     The mean and standard deviation of the function values of the
!     current simplex are now calculated.
!
        hstd=0.0
        hmean=0.0
        DO 1060 i=1,np1
          hstd=hstd+h(i)*h(i)
 1060     hmean=hmean+h(i)
        hmean=hmean/FLOAT(np1)
        hstd=(hstd-FLOAT(np1)*hmean*hmean)/FLOAT(np1)
        hstd=SQRT(ABS(hstd))
!
!     The parameter values (F) at the centroid of the current
!     simplex and the function value there (FUNC) are now
!     calculated.
!
        DO 1080 i=1,nop
          f(i)=0.0
          DO 1070 j=1,np1
 1070       f(i)=f(i)+g(j,i)
          f(i)=f(i)/FLOAT(np1)
 1080   CONTINUE
        CALL givef (f, func, dcoeff, val, clint, pd, ps, stt, tcov,
     &    thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     &    ltmax, ltmin, nclass, numa, numo, nvals, ns, msfail, 
     &    mtest, graph_file)
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
!     maximum number of iterations.
!
 1090   mfail=mfail+1
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
        GO TO 1380
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
        IF (hmean.eq.0) GO TO 1240
        test=savemn/hmean
      IF ((test.gt.0.9999995).and.(test.lt.1.0000005)) GO TO 1240
        iflag=0
        loop=0
        GO TO 740
!
!     If a calculated value of F(4) [=dmax] is greater than the
!     upper limit (DLIM) set by the program, F(4) is set at
!     DLIM, STEP(4) is set at zero, and the search goes
!     back to the beginning.
!
 1240   IF (dlim.ge.f(4)) GO TO 1250
        f(4)=dlim
        step(4)=0.
        GO TO 740
!
!     If JPRINT=1 the program prints out the results of each successful
!     convergence on a minimum.
!
 1250   IF (jprint.eq.1) THEN
          WRITE (2, 1270) neval
 1270     FORMAT (5(/),' PROCESS CONVERGES ON MINIMUM AFTER ', i4,
     &' FUNCTION EVALUATIONS'///)
          WRITE (2, 1280) (f(i), i=1, nop)
 1280     FORMAT (' MINIMUM AT   ',4(1x,e13.6))
          WRITE (2, 1290) func
 1290     FORMAT (//' MINIMUM FUNCTION VALUE   ',e13.6)
          WRITE (2, 1300)
 1300     FORMAT (///' END  OF  SEARCH'/1x,15('*'))
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
        CALL givef (f, func, dcoeff, val, clint, pd, ps, stt, tcov,
     &    thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     &    ltmax, ltmin, nclass, numa, numo, nvals, ns, msfail, 
     &    mtest, graph_file)
!
!
        kprint=0
!
!
!     The estimated 'best fit' values of the parameters from the current
!     pass through Loop 1410 are now computed, tabulated and printed out.
!
!     A density estimate (DEN) is calculated from D2L=F(3) by
!     correcting units to no./ha, and dividing by 2L in the case of
!     line transect data, and by 2Vt in the case of fixed-point data.
!
        IF (km.gt.0) THEN
          den(bootstrap)=(1.e6*f(3))/(sns*dist)
          GO TO 1370
        END IF
!
        IF (ifx.eq.0) THEN
        den(bootstrap)=(1.e4*f(3))/(sns*dist)
        ELSE
        den(bootstrap)=(1.e4*f(3))/(2.*rate*durn*ps)
        END IF
!
!     Values of the other parameter estimates are obtained by redefining
!     the estimates of F(1) and F(2), thus:
!
 1370   coeff1(bootstrap)=f(1)
        coeff2(bootstrap)=f(2)
        coeff3(bootstrap)=dcoeff
!
!
!     Running totals of DEN, COEFF1 and COEFF2 are now made to
!     enable calculation of mean values after the loop ends.
!
        tden=tden+den(bootstrap)
        tcoeff1=tcoeff1+coeff1(bootstrap)
        tcoeff2=tcoeff2+coeff2(bootstrap)
        tcoeff3=tcoeff3+coeff3(bootstrap)
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
        DO 1390 i=1,np1
          DO 1390 j=1,nop
 1390     g(i,j)=0.0
!
!     F and STEP values are reset to their original values also.
!
        DO 1400 ie=1,nop
          f(ie)=ft(ie)
          step(ie)=stept(ie)
 1400   CONTINUE
!
!
 1410 CONTINUE
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
      DO 1430 bootstrap=1, maxjb
        IF ((den(bootstrap)) .ne. 0.) THEN
           dendif=(den(bootstrap)-estden)**2.
           dsum=dsum+dendif
        END IF
 1430 CONTINUE
      sden = SQRT(dsum/(numest-1))
!
      cf1sum=0.0
      DO 1440 bootstrap=1, maxjb
        IF ((coeff1(bootstrap)).eq.0.) GO TO 1440
        cf1dif=(coeff1(bootstrap)-coeffnt1)**2.
        cf1sum=cf1sum+cf1dif
 1440 CONTINUE
      scf1 = SQRT(cf1sum/(numest-1))
!
      cf2sum=0.0
      DO 1450 bootstrap=1, maxjb
        IF ((coeff2(bootstrap)).eq.0.) GO TO 1450
        cf2dif=(coeff2(bootstrap)-coeffnt2)**2.
        cf2sum=cf2sum+cf2dif
 1450 CONTINUE
      scf2=sqrt(cf2sum/(numest-1))
      cf3sum=0.0
      DO 1460 bootstrap=1, maxjb
        IF (coeff3(bootstrap).eq.0) GO TO 1460
        cf3dif=(coeff3(bootstrap)-coeffnt3)**2.
        cf3sum=cf3sum+cf3dif
 1460 CONTINUE
      IF ((numest-msfail-1).le.0) GO TO 1470
      scf3=sqrt(cf3sum/(numest-msfail-1))
!
!
!     The products of the program are now output.
!
!
!     The number of animal groups in the selected range (NGROUPS) is 
!     first, followed by the estimated topographical cover (TCOV), 
!     then the number of actual parameter estimations made (NUMEST).
!
 1470 WRITE (2,1475) ngroups
 1475 FORMAT (/' Number in Groups in Distance Range = ',i4)
!
      IF (ifx.eq.0) WRITE (2,1476) numain
 1476 FORMAT (' Number of Individuals Detected Ahead = ',i5)
!
      IF (ifx.eq.0) WRITE (2,1477) numoin
 1477 FORMAT (' Number Overtaking (distance unmeasured) = ',i4)
!
      WRITE (2,1478) thh
 1478 FORMAT (' Height Difference from Eyelevel = ',f5.1,' m')
!
      IF (ifx.eq.0) WRITE (2,1479) estj
 1479 FORMAT (' Movement Correction Factor (J) = ',f6.3)
!
      If (ifx.eq.0) WRITE (2,1480) dist
 1480 FORMAT (' Adjusted Transect Length (LJ) = ',f10.3,' m')
!
      WRITE (2,1481) tcov
 1481 FORMAT (' Topographical Cover Value = ',f6.4)
!
      WRITE (2,1482) estdmax
 1482 FORMAT (' Calculated Maximum Detection Distance =',f8.1,' m')
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
!     The table is now ruled off.
!
 1840 WRITE (2,1850)
 1850 FORMAT (' x',4(18x,'x')/' ',77('x')/)
!
      WRITE (2,1860) coeffnt3,scf3
 1860 FORMAT (1x,'Detectability Coefficient (S) =',f8.2,', SE =',f6.2)
!
!
!     Key model estimates are now output if ISHOW was originally set
!     at 1.  The best estimates of F(1), F(2) and F(3) must be entered
!     in the subroutine GIVEF to do this.  The totals in the class
!     intervals are those from the initial data, derived from the
!     saved array variable VALT().  The output is then produced
!     within the subroutine.  The same is done with NUMO and NUMA.
!
!
      f(1)=coeffnt1
      f(2)=coeffnt2
      DO 1870 jv=1,nclass
 1870   val(jv)=valt(jv)
!
        numo=numoin
        numa=numain
!
      IF (km.eq.0) GO TO 1890
      f(3)=(estden*dist*pd*sns)/1.e6
      GO TO 1920
!
 1890 IF (ifx.eq.0) THEN
         f(3)=(estden*dist*pd*sns)/1.e4
       ELSE
 1910    f(3)=estden*2.*ps*rate*durn*pd/1.e4
      END IF
!
!
 1920 CALL givef (f, func, dcoeff, val, clint, pd, ps, stt, tcov,
     & thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     & ltmax, ltmin, nclass, numa, numo, nvals, ns, msfail, 
     & mtest,graph_file)
!
!
      params%complete = 1
      params%estden = estden
      params%sden = sden
      CLOSE (unit=2)
      RETURN
!
!
 1930 FORMAT (' Error opening ',a40,' - IOS = ',i6)
      RETURN
 1940 WRITE (6,1930) outfile,ios
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
      SUBROUTINE givef (f, func, s, val, clint, pd, ps, stt, tcov,
     & thh, vgh, sns, ifx, imv, iry, ishow, kdt, kprint,
     & ltmax, ltmin, nclass, numa, numo, nvals, ns, msfail, 
     & mtest, graph_file)
!
!     Double precision for all real numbers is set, together with
!     common values, dimensions and symbols for key variables.
!
!
      IMPLICIT NONE
!
!
      INTEGER iermax,ifx,imv,ios,iry,ishow,jj,jl,jrlow,jv
      INTEGER kdt,kprint,l10,l20,limit
      INTEGER max,md,msfail,mtest,ns,numa,numo,nvals,nclass
      CHARACTER*(*) graph_file
      DOUBLE PRECISION aprexi,auc,cint,clint,d2l,dd,dds,ddsm,dh,dif
      DOUBLE PRECISION difsq,dint,dl,dmax,dnr,dnrl,dnrh,dvg,e,ed
      DOUBLE PRECISION corrn,ermax,err,expd,expdr,expdy,expdv,func
      DOUBLE PRECISION hcint,htot,ltmin,ltmax,obsd,p,pa,pad,pam,pd
      DOUBLE PRECISION pr,prc,prr,prmax,ps,q,qdd,qdmax,qmin,qr
      DOUBLE PRECISION rlow,rmax,rr,s,sns,ss,ssh,ssl,ssmax,stt,tcov
      DOUBLE PRECISION texpd,thh,topdd,topmax,tot,tote,tr,vegdd
      DOUBLE PRECISION vegmax,vgh,vh,visdd,vismax,vl,vlowm,w,wdifsq
      DOUBLE PRECISION wh,vhowh,wl,wm,wtot,yh,yl,yy,yyh,yyl,z,zh,zl
      DOUBLE PRECISION val(80)
      DOUBLE PRECISION rout(80),calcnr(80),obsdnr(80)
      DOUBLE PRECISION yout(80),calcny(80),obsdny(80)
      DOUBLE PRECISION g(5,4),step(4),stept(4),f(4)
      DOUBLE PRECISION h(4),pbar(4),pstar(4),pstst(4)
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
!     highly improbable values of the parameters.
!
      htot=0.
      dmax=f(4)
!
!     HTOT is set high (1 000 000) if DMAX < 1
!
      IF (dmax.lt.1.) THEN
        htot=1.e+6
      END IF
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
!     HTOT is set to a higher value if D2L<0.01
!
      IF (d2l.lt.0.001) THEN
        htot=1.e+6+htot
      END IF
!
!     HTOT is set high if a>400
!
      IF (p.gt.400) THEN
        htot=1.e+6+htot
      END IF
!
!     HTOT is set high if a<0.01
!
      IF (p.lt.0.01) THEN
        htot=abs((p-2.)*1.e06)+htot
      END IF
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
!     An upper limit of 0.4 is set for the attenuation coefficient
!
  110 IF (q.ge.0) THEN
      IF (q.lt.0.4) GO TO 140
      htot=1.e+6+htot
      GO TO 140
      END IF
!
!     A lower limit of -2./DMAX is set to the attenuation coeffnt.
!
      qmin=-2./dmax
      IF (q.gt.qmin) GO TO 140
      q=qmin
!
!     The theoretical probability of detecting an individual at
!     the maximum recognition distance [PRMAX=Pr(rmax)] is now
!     calculated.
!
!     PA is the square of the conspicuousness coefficient.
!
!
  140 pa=p*p
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
!     If the calculated cover proportion is less than zero, the
!     program adds an arbitrary 100 to the accumulating total
!     to deflect the search away from such unreal values.
!
  170 IF (q.ge.0.) GO TO 180
      htot=1.e+2+htot
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
!     Loop 1180, which calculates expected values and compares them
!     with observed values, now begins .....
!
!
      DO 1180 jj=1,l10
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
        IF (wl.gt.rmax.and.wl.gt.ermax) GO TO 350
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
!     Loop 870, which calculates the expected number in each
!     arc within the class, and adds them together to produce
!     a progressive total (EXPD and TEXPD), now begins .....
!
!
  450   DO 870 jl=1,l20
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
  500     IF (zl) 520,510,520
  510     vlowm=0.
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
!     PRC [=P(r)] is the height of the detectability curve at DD,
!     corrected by subtracting PRMAX.
!
  680     prc=pr-prmax
!
!     Because 0<PRC<1, an upper limit of 1 and a lower limit of
!     0 are set to the probability function.
!
          IF (prc.le.1.) GO TO 690
          prc=1.
!
  690     IF (prc.ge.0) GO TO 700
          prc=0
!
!     The probability (PRR=g(r)) of detection within an arc of width
!     CINT is different from that in an arc of unit width, and is
!     equal to the product of the corrected probability and the
!     arc width.
!
  700     prr=prc*cint
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
!     E [= g(r).P(r)] is a probability density function which
!     describes the probability that an animal is both present and
!     detected in an observing arc.  It is calculated by multiplying
!     together the probability of detection in the arc [PRC=g(r)]
!     and the probability [QR=P(r)] that an individual is still
!     available for detection in that arc. TOTE is the accumulating
!     total of the probability density function.
!
  720     e=prr*qr
          tote=tote+e
!
!     The subroutine now calculates the number of detections (ED)
!     expected for radial, fixed-point and perpendicular distance data.
!     For radial data, computation is based on the assumption
!     that the arc of
!     radius r and width delta-r sweeps over a series of plots of
!     unit area.  The total number of detections expected in such
!     a plot will be [plot area] x [apparent density] x [probability
!     of detection].
!
!     With perpendicular distance data, ED is the expected number
!     detected at a distance r from the observer in the strips
!     CINT units wide at a perpendicular distance y=r from the
!     transect line.  Because there are one or two such strips, each
!     L times the area of the individual plot at distance r,
!     the expected no. = [DxNSxLJ=d2l] x Pd x CINT x g(r) x P(r).
!     NS must be converted to a floating-point number (SNS) first.
!
      IF (ifx.eq.0 .and. iry.ge.1)  THEN
          ed=d2l*cint*e*pd
!
!     With radial distance data, each observing arc sweeps out an
!     area L units long and 2r units wide, within which there are
!     2r/CINT x L individual plots at radial distance r.  Thus the
!     expected total number detected  =  [number expected in one
!     plot]x[number of plots]=[DxNSxLJ=d2l] x Pd x r x g(r) x P(r).
!
        ELSE IF (ifx.eq.0 .and. iry.lt.1) THEN
          ed=d2l*e*dnrh*pd
!
!     With fixed-point data, the expected total number detected
!     = [2Vt(=D2L)] x Pd x Ps x r x g(r) x P(r).
!
        ELSE
          ed=d2l*pd*e*dnrh
!
      END IF
!
!     With radial distance data, the total number expected in a
!     class (EXPD) is obtained by progressively adding together
!     the ED values from each arc with radii between DNRH and DNRL.
!
  790     expd=ed+expd
!
!     In the case of perpendicular distance data, all the EXPD
!     values within a class are added together to give a TEXPD
!     value for the class (while EXPD continues to accumulate
!     from class to class).
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
          IF (qr.gt.0.001) GO TO 870
!
        IF (jrlow.le.0) THEN
!
!     RLOW is the inner boundary of the arc by which 99.9% of the
!     detections are expected to have been made.
!
          rlow=dnrl-cint
          jrlow=1
        END IF
!
!     Loop 870 now ends.
!
  870   CONTINUE
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
  975   FORMAT (/91HIn order: P(r), g(r), Q(r), g(r).Q(r), SUM[g(r).Q(r)
     &], n(r), E{N(r)}, E{N(y)}, SUM[d.c.(S)]/)   
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
 1025   FORMAT (/91HIn order: P(r), g(r), Q(r), P(r).g(r), SUM[P(r).g(r)
     &], n(r), E{N(r)}, E{N(y)}, SUM[d.c.(S)]/)   
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
          WRITE (2,1070) prc,prr,qr,e,tote,ed,expd,expdv,s
 1070     FORMAT (1x,9(f12.6,3x))
        END IF
!
!
!     The difference between each observed (OBSD) and expected
!     value (EXPDV) for the class is now calculated.
!
!
 1080   dif=obsd-expdv
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
!     Loop 1180 now ends.
!
!
 1180 CONTINUE
!
!
!     If the search for a minimum is complete, and KPRINT=1,
!     the program now prints out a value for RLOW and a file
!     which tabulates observed and calculated frequencies.
!
      IF (kprint.eq.0) GO TO 1300
!
      WRITE (2,1200) rlow
 1200 FORMAT (' 99.9% r value (rmin) =',f7.2,' m ')
      OPEN (UNIT=3,FILE=graph_file,STATUS='NEW',IOSTAT=ios,
     & ERR=1320)
      IF (iry.gt.0) GO TO 1250
 1210 WRITE (3,1220)
 1220 FORMAT (3x,'   Midpt.   Calculated     Observed   ')
      DO 1240 jj=1,l10
        jv=l10-jj+1
      IF (kdt.gt.1) THEN
         limit = int((kdt-stt) / clint)
      ELSE
         limit = int((dmax-stt) / clint)
      END IF
      IF (jj.gt.limit) GO TO 1290
      WRITE (3,1230) rout(jv),calcnr(jv),obsdnr(jv)
 1230 FORMAT (5x,f6.1,6x,f7.2,7x,f6.1,4x)
 1240 CONTINUE
      GO TO 1290
 1250 WRITE (3,1260)
 1260 FORMAT (3x,'   Midpt.   Calculated     Observed   ')
      DO 1280 jj=1,l10
        jv=l10-jj+1
      IF (kdt.gt.1) THEN
        limit = int((kdt-stt) / clint)
      ELSE
        limit = int((dmax-stt) / clint)
      END IF
      IF (jj.gt.limit) GO TO 1290
      WRITE (3,1270) yout(jv),calcny(jv),obsdny(jv)
 1270   FORMAT (5x,f6.1,6x,f7.2,7x,f6.1,4x)
 1280 CONTINUE
 1290 CLOSE (unit=3)
!
!     The totalled sum of squares of the differences between
!     observed and calculated values is now defined as FUNC.
!
 1300 func=tot+htot
      GO TO 1340
!
!     The next two lines print out an error message if needed.
!
 1320 WRITE (6,1330) graph_file,ios
 1330 FORMAT (' Error opening ',a40,' - IOS = ',i6)
!
!     If F(2) has been given a high negative value, FUNC is set to
!     a high value before returning to the main program.
!
 1340 IF (f(2).ge.qmin) GO TO 1350
      func=(abs(f(2)-qmin-1.))*func
!
 1350 IF (kprint.eq.1) THEN
        WRITE (2,1370) func
 1370   FORMAT (1x,'OLS Difference at Minimum = ',f15.6/)
      END IF
!
      RETURN
      END SUBROUTINE givef
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
      DOUBLE PRECISION jseed,ifrst
      COMMON /seed/ jseed,ifrst
!
      jseed=iseed
      ifrst=0
!
      END SUBROUTINE srandnag

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
      DOUBLE PRECISION result

      PARAMETER  (mplier=16807,modlus=2147483647,mobymp=127773,
     & momdmp=2836)
!
      DOUBLE PRECISION jseed,ifrst
      COMMON /seed/ jseed,ifrst
      INTEGER hvlue,lvlue,testv,nextn
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
      END BLOCKDATA randbd
