*     PROGRAM MNPS2
*
*     This program is designed to return population density estimates
*     from 'distance' data collected using either line transect or fixed
*     point methods, together with estimates of other parameters of the
*     observing situation.  It has been designed primarily for use in 
*     ground surveys of mammal and bird populations in terrestrial 
*     habitats, though a variant can work with aerial surveys using both
*     fixed-wing aircraft and helicopters.
*
*     [The program as written does not allow for normal aerial survey 
*     density estimation, or estimation where a segment of the viewing   
*     arc is missing; for these types of data, program MNPSAS should 
*     be used instead.]
*
*     MNPS2 works by modelling mathematically the frequency distributions 
*     of the numbers sighted under uniform observing conditions over a  
*     range of distances from either the observer or the transect line. 
*     Curve-fitting is achieved by using the simplex method to minimize 
*     the sum of squares of the differences between the observed and 
*     calculated numbers across the frequency distribution.
*                 
*     Based on the 'best-fit' solution so obtained, the program returns  
*     population density estimates, estimates of other parameters of the
*     model and related information on the system in which the original 
*     observations were made. Standard errors of the parameters are 
*     obtained by bootstrap sampling of data from the original data set.
*
*
*     MNPS2 has been written by David Morgan, Department of Zoology,
*     University of Melbourne, Vic. 3010, Australia.
*
*
*     ************************************************************
*
*
*     MNPS2 consists of a MAIN program and two subroutines,  
*     SUBROUTINE GIVEF and SUBROUTINE SRANDNAG.
*
*
*     The main program determines parameter values for a least
*     squares fit between calculated and observed frequencies,
*     using the simplex method in a way derived from 
*     Nelder & Mead, The Computer Journal, January, 1965.  [Original 
*     'MINIM' program by D.E.Shaw, Divn. of Mathematical Statistics, 
*     CSIRO, Australia.]
*
*     Subroutine GIVEF calculates the sum of squares by using the
*     mathematical model, and returns this to the main program.
*     Subroutine SRANDNAG supplies a pseudo-random number used in the
*     data resampling procedure in response to a seed value supplied 
*     by the main program.
*
*
*     The program accepts several alternative types of data and
*     provides three different ways of calculating expected values, 
*     viz:
*
*		. radial distance data supplied (IFX=0 and IRY=0);
*		. perpendicular distance data supplied (IFX=0 and IRY=1);
*		. fixed-point data supplied (IFX=1 and IRY=0);
*
*
*     Data are supplied to the program by means of a data file, 
*     'IN<num>.dat'.  The parameters supplied in the input via the 
*     data file are, in sequence:
*
*
*     A header line - to label and describe the data set.
*
*     NVALS - the total number of individual animal detections
*          in the sample submitted for modelling.
*
*     CLINT - the class interval to be used within the program, 
*          chosen so that 80 x CLINT is less than the observed
*          maximum detection distance.
*
*     STT - the detection distance value (r or y) from which class 
*          intervals begin (usually zero).
*
*     NUMA - the total number of individual animals detected 
*          ahead of a moving observer on a transect.
*
*     NUMO - the number of individuals overtaking the observer
*          from behind.
*
*     DIST - the transect length corrected for the effect of
*          animal movement (=LJ).
*
*     THH - the vertical distance between observer eye level
*          and the median horizontal plane occupied by
*          the population, best represented by the root mean
*          square of the measured vertical differences.
*
*     LTMIN - the approximate minimum distance at which an animal
*          may be obscured from the observer by topography.
*
*     LTMAX - the p=0.001 upper-limit distance at which an animal
*          may be detected by the observer in the absence of
*          vegetation cover.
*
*     IFX - a parameter to modify the computations appropriately
*          where data are collected by fixed observer sampling, viz:
*
*            =0  handles data from transects;
*            =1  handles data from fixed point sampling.
*
*     IRY - a parameter to control the computations within
*          the subroutine for transect data, viz:
*
*            =0  handles N(r) radial distance transect data;
*            =1  handles N(y) perpendicular distance data.
*
*     NS - the number of sides of the transect line used during
*          a survey, viz:
*            
*            =0  not transect data, so not relevant (fixed point);
*            =1  handles data from only one side of the line;
*            =2  handles data from both sides of the transect line.
*
*     KM - a parameter to modify the program if distance are
*          measured in kilometres, viz:
*
*            =0  distances are measured in metres;
*            =1  distances measured in kilometres.
*
*         The program as written assumes that all distances are 
*         in metres.
*
*    IMV - a parameter used to determine which method of 
*		  calculating probability is used in the subroutine.
*		  It either:
*			. calculates a mean probability value
*			      within each interval (IMV=0); or
*		        . calculates the median probability
*			      value within each interval (IMV=1), viz: 
*
*            =0  uses the mean value of P(r) in an interval;
*            =1  uses the median value of P(r).
*
*          [If perpendicular distances are being supplied as raw 
*          input data instead of providing radial distances and 
*          angles, put IMV=1 as well as IRY=1.  This removes the need
*          for data on angles, but also precludes the possibility
*          of calculating median probabilities when such data are
*          supplied.   (The program resets IMV=0 once the data have
*          been entered.)]
*
*          IMV is set at 1 within the program if the attenuation
*          coefficient becomes negative or if the data are
*          visual observations.
*
*     KDT - a control parameter to indicate the type of
*          observational data supplied to the program, viz:
*
*            = 0  for visual and flushing data;
*            = 1  for auditory data.
*
*          If KDT=0, the program sets IMV=1 within the
*          subroutine.
*
*     IPRINT - a parameter to control the output of progress 
*          evaluations from the main program, viz:
*
*            = 0  no progress evaluations;
*            = 1  reports initial convergence (with the parameter 
*                 and function values there), overrunning of MAX, 
*                 and a progress report on the minimization every 
*                 function evaluation;
*
*     JPRINT - a second parameter to control the progress of the
*          main program, viz.:
*      
*            = 0  no frequency distributions or internal summaries;
*            = 1  prints the frequency distributions of the original
*                 and backstrapped data and the final outcome of
*                 each data set function minimization;
*
*     ISHOW - a control parameter which outputs in sequence: PRC,
*          PRR, QR, E, TOTE (total E), ED, EXPD, EXPDV and the
*          detectability coefficient S.
*
*            = 0  no output;
*            = 1  produces output.
*
*     MAXJB  -  the number of sets of iterations used to derive
*          the backstrapped series of model parameters and
*          calculate standard errors.  This is best set at between
*          500 and 2000.
*
*     R3S - three times the estimated interquartile range of all
*          the deviations between the calculated and observed
*          values, used in the process of converging on an
*          internal minimum.  Its default value is 100 if no R3S 
*          value is supplied in the input.
*
*     VGH - the approximate average height of vegetation cover
*          in the animal's habitat in situations where the observer 
*          is well above the plane of the population and most of
*          the line of detection is unobstructed.
*
*     IT - the duration of a fixed-point census (in min.).
*
*     IV - the mean speed of animal movement (in m/min), used in
*          fixed observer censuses only.
*
*     PD - the proportion of the population observable at the
*          time of a census, taken as =1 as the default.
*
*     PS - in fixed-point censuses, the proportion of the circle
*          scanned by the observer during the census.
*
*
*     F  - the model parameter values used in the search for the
*          minimum point.  Four initial F values are required:
*
*          F(1):  an estimate of the conspicuousness coefficient
*                 of the species under the conditions of the
*                 census;
*
*          F(2):  an estimate of the attenuation coefficient of
*                 the model under the census conditions;
*
*          F(3):  an estimate of the population density (D),
*                 expressed in number of individuals per hectare;
*
*          F(4):  an estimate of the maximum distance (in m) from the
*                 observer at which species recognition is
*                 possible under the conditions of the census.
*
*          Those supplied with the data are the starting values
*          for the search; those printed on exit are the parameter 
*          values which specify the minimum point.  Initial F()
*          values supplied are best based on existing knowledge.
*
*
*     STEP - the step sizes used initially in the function 
*          minimization process to modify the F values. Typical
*          STEP values are:
*
*          STEP(1):  just under 50% of the initial F(1) value;
*
*          STEP(2):  equal to the initial F(2) value;
*
*          STEP(3):  half the initial F(3) value;
*
*          STEP(4):  usually set at zero.
*
*          Setting the STEP value at 0.0 for a parameter fixes the F()
*          value at that supplied to the program.  Where a data set (NVALS) 
*          is small (say, < 50 detections), STEP(2) should be set at 0.0. 
*          Where the set is very small (say, < 30 detections), STEP(1)
*          should also be set at 0.0 .
*
*
*    R() - the radial detection distances (r) originally measured
*          in the field, submitted in the order in which the data
*          were collected.
*
*    NSIZE() - the size of each cluster (group) of animals at the
*          moment of detection, submitted in the same order as the
*          corresponding R() values.
*
*    ANGLE() - the horizontal angle between the direction of a transect
*          and the bearing to a cluster of animals at the moment of
*          detection, and measured in degrees.  These data are
*          required if computations are to be based on perpendicular
*          distance modelling; otherwise they are not needed.
*
*
*    In addition, a number of temporary variables are used within
*    the program.  These include:
*
*     NOP - the number of parameters potentially varied during
*          the minimization process.
*
*     FUNC - the sum of squares of the differences between 
*		  observed and expected values;
*
*     MAX - the maximum number of function evaluations to be
*          allowed (set at 750);
*
*     Three coefficients used in the search for a minimum value, viz:
*
*          A - a reflection coefficient;
*          B - a contraction coefficient; and
*          C - an expansion coefficient.
*
*     STOPC - a stopping criterion used to initiate a test for
*          convergence on a minimum point.
*
*     NLOOP - convergence is tested for every NLOOP times the
*          process changes the simplex.  After initial
*          convergence, NLOOP further changes are allowed
*          before testing for final convergence.
*
*     KWT - a control parameter which directs calculation of
*          a biweighted least squares once the number of 
*          iterations exceeds a predetermined value.
*
*            = 0  for normal least squares computation,
*            = 1  for biweighted least squares.
*
*     LPRINT - a control parameter to remove the biweighting
*          procedure.
*
*            = 0  except for final computation of FUNC;
*            = 1  when FUNC comes from simplex method.
*
*     MFAIL and MSFAIL - two variables used to count the numbers of
*           times when convergence fails (MFAIL) or when the
*           detectability coefficient (S) cannot be calculated (MSFAIL).
*           They are used to correct the numbers of counts used in
*           computing means and standard deviations of key parameters.
*
*     MTEST - a variable used to flag a successful run through Loop 855,
*           and assist in the computation of MSFAIL at the end of
*           Subroutine GIVEF.
*
*     D2L - the product of the population density and twice the
*           total transect length - a convenient variable in
*           calculations.
*
*     ESTDEN  - the estimated population density ('D').
*
*     P  - the conspicuousness coefficient ('a').
*
*     Q  - for auditory data (where KDT=1), the 
*          conspicuousness coefficient ('b').
*  
*        - for visual data (where KDT=0), the mean vegetation 
*          cover proportion in the habitat between animal and
*          observer.
*
*     S  - a coefficient of detectability, usable in density estimation
*          by direct calculation.
*
*     TCOV - the estimated proportion of topographical cover in
*          the line of sight between observer and animal.
*
*     DMAX - the maximum direct-line detection distance ('dmax').
*
*     RMAX - the maximum horizontal detection distance ('rmax').
*
*     ERMAX - a maximum detection distance calculated by the 
*          program.
*
*
*     Further explanatory notes are included within the program below.
*
*
*     This program uses double precision for all real numbers.
*
*
*
*     PROGRAM MNPS2

      SUBROUTINE CALCULATE_DENSITY(PARAMS, HEADER, OUTFILE)

*     This seqence derived type must be kept in complete agreement with the
*     C-language definition in mnps2.h
      TYPE CALC_PARAMS
      SEQUENCE        ! SEQUENCE indicates not to insert internal
                      ! padding for data alignment
      INTEGER NVALS,NUMA,NUMO,LTMIN,LTMAX,IFX,IRY,NS,KM,IMV,KDT
      INTEGER IPRINT,JPRINT,ISHOW,MAXJB,IT,IV
      DOUBLE PRECISION CLINT,STT,DIST,THH,R3S,VGH,PD,PS,F(4)
      DOUBLE PRECISION STEP(4),R(5000)
      INTEGER NSIZE(5000)
      DOUBLE PRECISION ANGLE(5000)
      END TYPE CALC_PARAMS

      TYPE (CALC_PARAMS) PARAMS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

*     CHARACTER*80 HEADER,INFILE*40,OUTFILE*40,RUNID*10
      CHARACTER*(*) HEADER, OUTFILE

      INTEGER NSIZE(5000),NBSZ(5000),I,IB,IC,IMAX,IN,IRB
      INTEGER ISEED,IPRINT,IMV,IRY,ISHOW,IT,IV,J,JB,JBSTP,JK,JR,JX
      INTEGER KDT,KWT,LTMIN,LTMAX,MAX,MAXJB,NS,NUMA,NUMEST,NUMO,NVALS 
      DOUBLE PRECISION VALT(80)
      DIMENSION R(5000),BSTR(5000),Y(5000),BSTY(5000),
     & DEN(5000),COEFF1(5000),COEFF2(5000),COEFF3(5000),ANGLE(5000)
      DIMENSION G(5,4),STEP(4),STEPT(4),F(4),FT(4)
      DOUBLE PRECISION VAL(80),H(4),PBAR(4),PSTAR(4),PSTST(4)
*     REAL RAND
*
*  
*     The program accepts up to 5000 data values, each being the total
*     number of observations [N(r) or N(y)] within the class
*     intervals, beginning with that nearest r=0 or y=0 [R(1),
*     NSIZE(1) and ANGLE(1), if supplied].  If ANGLE() is not
*     supplied and IRY=0 or =2, computation still proceeds.
*
*
*     Some variables are common to the main program and the
*     subroutine, viz.
*
*     COMMON / common1 /VAL(80),H(4),PBAR(4),PSTAR(4),PSTST(4)
*     COMMON / common2 /HSTAR,HSTST,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH
*    &,S,SNS
*     COMMON / common3 /KM,IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT
*     COMMON / common4 /LTMAX,LTMIN,MAX,NUMA,NUMO,NVALS,NS,MSFAIL,MTEST

*     This program is designed to receive a list-directed data file with
*     a name constructed from RUNID (a character value not more than 10
*     characters long).  The filename constructed will be IN'runid'.DAT.
*     RUNID is also used to construct the output file name of the form
*     OUT'runnid'.DAT.
*
*     The header line should begin with the computer run number each time
*     to avoid confusion, and be bounded by quotation marks.  Other entries
*     must be separated either by commas or spaces (and pre-checked before
*     running).
*
*
*     The program is set to seek a minimum using up to MAX=750
*     iterations, begin with four parameters (NOP) and set a testing
*     point (NLOOP) equal to unity.
*
      MAX=750
      NOP=4
      NLOOP=1

      NVALS  = PARAMS.NVALS
      CLINT  = PARAMS.CLINT
      STT    = PARAMS.STT
      NUMA   = PARAMS.NUMA
      NUMO   = PARAMS.NUMO
      DIST   = PARAMS.DIST
      THH    = PARAMS.THH
      LTMIN  = PARAMS.LTMIN
      LTMAX  = PARAMS.LTMAX
      IFX    = PARAMS.IFX
      IRY    = PARAMS.IRY
      NS     = PARAMS.NS
      KM     = PARAMS.KM
      IMV    = PARAMS.IMV
      KDT    = PARAMS.KDT
      IPRINT = PARAMS.IPRINT
      JPRINT = PARAMS.JPRINT
      ISHOW  = PARAMS.ISHOW
      MAXJB  = PARAMS.MAXJB
      R3S    = PARAMS.R3S
      VGH    = PARAMS.VGH
      IT     = PARAMS.IT
      IV     = PARAMS.IV
      PD     = PARAMS.PD
      PS     = PARAMS.PS
*     RUNID  = PARAMS.RUNID

      DO 2 IG=1,NOP
      F(IG)    = PARAMS.F(IG)
      STEP(IG) = PARAMS.STEP(IG)
    2 CONTINUE

      DO 24 IH=1,5000
      R(IH)     = PARAMS.R(IH)
      NSIZE(IH) = PARAMS.NSIZE(IH)
      ANGLE(IH) = PARAMS.ANGLE(IH)
   24 CONTINUE

*     The next set of commands names the input and output files.
*
*     WRITE(6,281)
* 281 FORMAT(' Enter run identifier : ',$)
*     READ(5,282)RUNID
* 282 FORMAT(A10)
*     LR=INDEX(RUNID,' ')
*     IF(LR.EQ.0)THEN
*       LR=10
*     ELSE
*       LR=LR-1
*     ENDIF
*     INFILE='IN'//RUNID(1:LR)//'.DAT'
*     OUTFILE='OUT'//RUNID(1:LR)//'.LIS'
*
*
*     All data lines are read into the program together at this point.
*     Care in required to ensure that the group size data is submitted
*     in PRECISELY the same sequence as the corresponding R(IN) data, 
*     and omit overtaking cases where r=0.
*           
*
*     OPEN(UNIT=1,FILE=INFILE,STATUS='OLD',IOSTAT=IOS,ERR=283)
      OPEN(UNIT=2,FILE=OUTFILE,STATUS='NEW',IOSTAT=IOS,ERR=284)
*     READ(1,*) HEADER,NVALS,CLINT,STT,NUMA,NUMO,DIST,THH,LTMIN,
*    & LTMAX,IFX,IRY,NS,KM,IMV,KDT,IPRINT,JPRINT,ISHOW,MAXJB,R3S,
*    & VGH,IT,IV,PD,PS,(F(I),I=1,NOP),(STEP(I),I=1,NOP),
*    & (R(IN),IN=1,NVALS),(NSIZE(IN),IN=1,NVALS),
*    & (ANGLE(IN),IN=1,(NVALS*(IRY-IMV)))
*
*
*     F(3) is now raised in value to approximate D2LJ in the case
*     of line transect data, or D2ut for fixed point data.  D is
*     also altered from no./ha to no./sq.m.
*
      IF (IFX) 105,105,106
  105 F(3)=(2.*DIST*F(3))/1.E4
      STEP(3)=(2.*DIST*STEP(3))/1.E4
      GO TO 107
  106 F(3)=(2.*IV*IT*F(3))/1.E4
      STEP(3)=(2.*IV*IT*STEP(3))/1.E4
*
*
*     If a value of the maximum detection distance F(4) has been
*     entered which exceeds 80 times the class interval, F(4) is
*     reset at 80 times that interval to avoid computation problems.
*
  107 IF (F(4).GT.(80*CLINT)) F(4)=80*CLINT
*
*
*     The header line now begins the program output.
*
      WRITE(2,4) HEADER
    4 FORMAT(1X,A80)
*
*
*     IRY is reduced by 1 for convenience in Fortran computations.
*
      IRY=IRY-1
*
*     The program prints out the class interval width (CLINT) and
*     either the total transect length (DIST) or the total time
*     spent (IT) at fixed points.
*
      IF (IFX) 5,5,6
    6 WRITE(2,3) CLINT, IT
    3 FORMAT(/23H Class Interval Width =,F7.1,
     &25H m.    Total Time Spent =,I5,5H min.)
      GO TO 15
*
    5 IF (KM) 7,7,9
    7 WRITE(2,8) CLINT, DIST
    8 FORMAT(/23H Class Interval Width =,F7.1,
     &28H m.    Total Distance (xJ) =,F10.3,3H m.)
      GO TO 15
*
    9 WRITE(2,10) CLINT, DIST
   10 FORMAT(/23H Class Interval Width =,F7.1,
     &28H m.    Total Distance (xJ) =,F10.3,4H km.)
*
*
*     If calculations are to be based on perpendicular distances (y)
*     from the transect line, and the data supplied are radial
*     distances and angles, perpendicular distances are now calculated
*     for the data entered initially,  this action being
*     prompted by IRY having a zero value.  If perpendicular 
*     distance data as such were supplied (IMV=1), this step is
*     bypassed.
*
   15 IF (IRY) 11,12,11
   12 IF (IMV) 14,14,21
   14 DO 13 IN=1,NVALS
      Y(IN)=R(IN)*SIN((ANGLE(IN)*3.14159265)/180.)
   13 CONTINUE
      GO TO 11
*
*     If perp. values were supplied in the input, R is redefined as Y.
*
   21 DO 22 IN=1,NVALS
      Y(IN)=R(IN)
   22 CONTINUE
*
*     If perpendicular distance data are being used, and IMV was set
*     at 1 because perpendicular distance data are being supplied
*     as original data (not as radial distances and horizontal angles),
*     then IMV is reset at 0.
*
      IF ((IRY.EQ.0).AND.(IMV.EQ.1)) IMV=0
*
*     The initial values of 'a', 'b', 'D2L, and 'dmax' are retained
*     as FT and the corresponding steps as STEPT to make possible 
*     reruns of calculations.
*
   11 DO 16 IA=1,NOP
      FT(IA)=F(IA)
      STEPT(IA)=STEP(IA)
   16 CONTINUE
*
*     DMAX is given an upper limit (DLIM) which is 80 times the
*     interval width (CLINT).
*
      DLIM=CLINT*80.
*
*     If transect lengths have been expressed in kilometres (KM=1),
*     distance data are converted to metres.
*
      IF (KM) 18,18,17
   17 DIST=DIST*1000.
*
*     The stopping criterion (STOPC) is set at a suitable value.
*
   18 STOPC=0.00001
*
*     If progress reports are required (IPRINT=1), the program
*     prints a heading for them.
*
      IF (IPRINT) 30,30,19
   19 WRITE(2,20) IPRINT                                                  
   20 FORMAT (22H PROGRESS REPORT EVERY,I4,21H FUNCTION EVALUATIONS/  
     &/24H EVAL. NO.  FUNC. VALUE ,10X,10HPARAMETERS)                       
*                                                                       
*     The term 'APPROX' is used to test closeness to zero.
*
   30 APPROX = 1.E-20
*                                                                       
*     IF NO VALUES OF A,B AND C ARE INPUT , I.E. A IS SET = 0.0 ,    
*     THEN THE PROGRAM SETS A = 1.0 , B = 0.5 , C = 2.0				
*
      IF (ABS(A).LT.APPROX) GO TO 40                                          
  31  GO TO 50
  40  A=1.0
      B=0.5                                                                  
      C=2.0                                                                  
*                                                                       
*     NAP is the number of parameters to be varied (i.e. with STEP
*     not equal to zero.
*
*
   50 NAP=0                                                       
      LOOP=0                                                                 
      IFLAG=0                                                                
      KPRINT=0                                                               
      LPRINT=0
*
*
*     If all STEP sizes have been set at zero, computation goes
*     to line 850, calculates the set of values resulting from
*     the values of a, b, D2L and dmax supplied, and ends.  Otherwise
*     it uses NAP to indicate the number of submitted parameters
*     to be varied during subsequent iterations.
*
      DO 70 I=1,NOP                                                     
      IF (ABS(STEP(I)).GT.APPROX) GO TO 60                                   
   51 GO TO 70                                                              
   60 NAP=NAP+1                                                         
   70 CONTINUE                                                          
      IF (NAP) 80,80,89
   80 KPRINT=1
      SNS=FLOAT(NS)
      IF (SNS.EQ.0.) SNS=2.
      IF (KM) 82,82,81
   81 ESTDEN=(1.E6*F(3))/(2.*DIST)
      GO TO 89
   82 IF (IFX) 83,83,84
   83 ESTDEN=(1.E4*F(3))/(2.*DIST)
      GO TO 89
   84 ESTDEN=(1.E4*F(3))/(2.*IV*IT)
*
*
*     To enable parameter estimation using bootstrapping, the basic 
*     MINIM routine is run a predetermined MAXJB times, beginning
*     with an initial run.  A number of functions are set at zero first.
*
   89 JBSTP=0
      JB=0
      MFAIL=0
      MSFAIL=0
      TDEN=0.0
      TCOEFF1=0.0
      TCOEFF2=0.0
      TCOEFF3=0.0
*
   90 DO 855 JB=1,MAXJB
*     
*     The flag variable MTEST is first set at zero.
*
      MTEST=0
*
*     The program now computes the first distribution of the numbers
*     detected within each class interval, based on the detection
*     distances [(R(IN)], the numbers in each group [(NSIZE(IN)] and
*     the class interval (CLINT) preset in the input.
*
*     In subsequent runs through Loop 855, bootstrapping applies  
*     (JBSTP=1) and a bootstrapped distribution is used instead.
*
*
*     The number detected within each class [VAL(IC)], is the sum of the
*     numbers in each class in which the R(IN) or Y(IN) values fall.
*     Calculating the various VAL(IC) values first requires
*     finding which R(IN) or Y(IN) values fall within the interval
*     concerned, then adding all the NSIZE(IN) values which fall 
*     within that class.  This will be done for each class interval 
*     in turn, beginning with the calculation of VAL(1) for the
*     nearest class to r=0 or y=0.
*
*     An alternative computation works with perpendicular distance
*     values (see below).
*     
*
      IF (IRY) 32,33,32
*
*
*     Either the initial r class totals, VALT(), are calculated...
*
*
   32 FRST=STT
      DO 92 IC=1,80
      VAL(IC)=0.0
      IF (JBSTP.EQ.1) GO TO 91
*
      DO 93 IR=1,NVALS
      IF (R(IR).GT.FRST .AND. R(IR).LE.(FRST+CLINT)) 
     & GO TO 85
      GO TO 93
   85 VAL(IC)=VAL(IC)+NSIZE(IR)
   93 CONTINUE
      FRST=FRST+CLINT
   92 CONTINUE
*
*     The frequency distribution of the original data is saved, 
*     as VALT(IG).
*
      DO 75 IG=1,80
      VALT(IG)=VAL(IG)
   75 CONTINUE
      IF (NAP.LE.0) GO TO 1531
*
      JBSTP=1
      GO TO 98
* 
*
*     Or the initial y class totals, VALT(), are calculated...
*
*     Absolute values of Y() are used because data from the
*     two sides of a transect line are pooled.
*
*   
   33 FRST=STT
      DO 34 IC=1,80
      VAL(IC)=0.0
      IF (JBSTP.EQ.1) GO TO 91
*
      DO 35 IR=1,NVALS
      IF (ABS(Y(IR)).GE.FRST .AND. ABS(Y(IR)).LT.(FRST+CLINT)) 
     & GO TO 36
      GO TO 35
   36 VAL(IC)=VAL(IC)+NSIZE(IR)
   35 CONTINUE
      FRST=FRST+CLINT
   34 CONTINUE
*
*     The frequency distribution of the original data is now saved, 
*     as VALT(IG).
*
      DO 37 IG=1,80
      VALT(IG)=VAL(IG)
   37 CONTINUE
      IF (NAP.LE.0) GO TO 1531
*
      JBSTP=1
      GO TO 98  
*
* 
*     To calculate a bootstrapped distribution, values of R(IN) and  
*     the corresponding group NSIZE(IN) are to be chosen at random
*     with replacement, based on a randomly-selected value of IN, 
*     which ranges between 1 and NVALS, the total number of groups of 
*     animals detected.  Loop 95 selects these values.
*
*
   91 ISEED=0
      DO 95 JR=1,NVALS
*     
*     A seed value is first chosen, based on the (sequential) values of
*     R in the data input, the size of the data set provided (NVALS),
*     the stage reached in the main backstrapping loop (JB), 
*     and the relevant step (JR) in Loop 95.  This should minimize
*     the likelihood of common seed values being submitted in
*     different random number searches.
*     In this way, no two seed values are likely either to be
*     identical or to vary in a consistent way.  They are multiplied
*     by a relatively large number (72493) to ensure the product
*     is a correspondingly large number.
*
*     An upper limit is placed on ISEED to prevent integer overflow.
*
      IF (ISEED.GT.150000) ISEED=ISEED/100
      ISEED = ISEED + (R(JR)*724*JR)/JB
*
*     A random number between 0 and 1 is now called from the Subroutine
*     SRANDNAG, multiplied by NVALS to give it a value between 0 and NVALS,
*     and converted to integer form, adding 1 to allow for chopping. 
*
      CALL SRANDNAG(ISEED)
      X = RANDNAG()
      JX = INT(X*NVALS + 1)
*
*     What happens now depends on whether radial or perpendicular
*     distance measurements are being used.
*
*
      IF (IRY) 41,42,41
*    
*     Either: radial distance computations are made.
*
*     The randomly-chosen bootstrapped value of R, R(JX), is now
*     redefined as BSTR(JR), the JRth bootstrapped value of R.  The
*     corresponding value of NSIZE, NSIZE(JX), is now redefined
*     as NBSZ(JR), the JRth bootstrapped value of NSIZE. Loop
*     95 then goes back to its beginning to select another R value,
*     and so on until all NVALS selections have been made.
*
   41 BSTR(JR)=R(JX)
      NBSZ(JR)=NSIZE(JX)
      GO TO 95
*
*     Or: perpendicular distance computations are made.
*
*     The randomly-chosen bootstrapped value of Y, Y(JX), is now
*     redefined as BSTY(JR), the JRth bootstrapped value of Y.  The
*     corresponding value of NSIZE, NSIZE(JX), is now redefined
*     as NBSZ(JR), the JRth bootstrapped value of NSIZE. Loop
*     95 then goes back to its beginning to select another Y value,
*     and so on until all NVALS selections have been made.
*
   42 BSTY(JR)=Y(JX)
      NBSZ(JR)=NSIZE(JX)
*
   95 CONTINUE
*
*
*     There should now be a new set of (N=NVALS) R (or Y) and NSIZE  
*     values to use in putting together a new frequency distribution.
*
*     The calculation path differs according to whether radial or
*     perpendicular distance data are being handled.
*
*
      IF (IRY) 43,44,43
*
*
*     Either: radial distance data are handled.
*
   43 FRST=STT
      DO 96 IC=1,80
      VAL(IC)=0
      DO 97 IRB=1,NVALS
      IF ((BSTR(IRB).GT.FRST) .AND. (BSTR(IRB).LE.(FRST+CLINT))) 
     & GO TO 86
      GO TO 97
   86 VAL(IC)=VAL(IC)+NBSZ(IRB)
   97 CONTINUE
      FRST=FRST+CLINT
   96 CONTINUE
*
*     BSTR(IRB) values need to be reassigned at this point 
*     or some values will be carried into subsequent loop
*     iterations.
*
      DO 88 JS=1,NVALS
      BSTR(JS)=0.0
      NBSZ(JS)=0
   88 CONTINUE	
      GO TO 98
*
*
*     Or: perpendicular distance data are handled.
*
   44 FRST=STT
      DO 45 IC=1,80
      VAL(IC)=0
      DO 46 IRB=1,NVALS
      IF ((BSTY(IRB).GE.FRST) .AND. (BSTY(IRB).LT.(FRST+CLINT))) 
     & GO TO 47
      GO TO 46
   47 VAL(IC)=VAL(IC)+NBSZ(IRB)
   46 CONTINUE
      FRST=FRST+CLINT
   45 CONTINUE
*
*     BSTY(IRB) values need to be reassigned at this point 
*     or some values will be carried into subsequent loop
*     iterations.
*
      DO 48 JS=1,NVALS
      BSTY(JS)=0.0
      NBSZ(JS)=0
   48 CONTINUE	
*
*
*     The calculated set of values, VAL(IC), in each class interval 
*     is now printed out for the set of data concerned if JPRINT=1.
*
*
   98 IF(JPRINT) 103,103,102
  102 WRITE(2,87) JB
   87 FORMAT(///27H  Bootstrap Replicate No. = ,I4/)
      WRITE(2,99) (VAL(I),I=1,80)
   99 FORMAT(8F10.4,//)
*
*                                                                       
*                                                                       
*     The initial simplex of program MINIM is now set up.
*
*
  103 DO 100 I=1,NOP
  100 G(1,I)=F(I)                                                       
      IROW=2                                                            
      DO 130 I=1,NOP                                                    
      IF (ABS(STEP(I)).LT.APPROX)  GO TO 130                                   
  110 DO 120 J=1,NOP                                                    
  120 G(IROW,J)=F(J)                                                    
      G(IROW,I)=G(IROW,I)+STEP(I)                                       
      IROW=IROW+1                                                       
  130 CONTINUE                                                          
      NP1=NAP+1                                                         
      NEVAL=0                                                           
      DO 170 I=1,NP1                                                    
      DO 140 J=1,NOP                                                    
  140 F(J)=G(I,J)
      CALL GIVEF (F,H(I),VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
      NEVAL=NEVAL+1                                                     
*                                                                       
*     All points in the initial simplex become output if IPRINT=1.
*
      IF (IPRINT) 170,170,150
  150 WRITE(2,160) NEVAL,H(I),(F(J),J=1,NOP)
  160 FORMAT (/3X,I4,4X,E13.6,8(1X,E13.6)/24X,8(1X,E13.6)/24X,4(1X,
     &E13.6))
  170 CONTINUE
*
*                                                           
*                                                                       
*     Now follows the basic loop.  That is, given a simplex, it
*     determines a new simplex and tests for convergence as
*     required (following the flow chart given in Nelder and
*     Mead).
*
*     HMAX and HMIN are the maximum and minimum function values
*     of the current simplex.
*
*
  180 LOOP=LOOP+1                                                       
      IMAX=1                                                                
      IMIN=1                                                                
      HMAX=H(1)                                                             
      HMIN=H(1)                                                             
      DO 220 I=2,NP1                                                    
      IF (H(I).GT.H(IMAX)) GO TO 190                                         
  181 GO TO 200                                                          
  190 IMAX=I                                                            
      HMAX=H(I)                                                         
  200 IF (H(I).LT.H(IMIN)) GO TO 210
  201 GO TO 220                                                            
  210 IMIN = I                                                          
      HMIN=H(I)                                                         
  220 CONTINUE                                                          
*                                                                       
*     The centroid of all vertices, excluding the maximum, is
*     now found.
*
      DO 230 I=1,NOP                                                    
  230 PBAR(I)=0.0                                                       
      DO 260 I=1,NP1                                                    
      IF (I.EQ.IMAX) GO TO 260                                                
  240 DO 250 J=1,NOP                                                    
  250 PBAR(J)=PBAR(J)+G(I,J)/FLOAT(NAP)                                   
  260 CONTINUE                                                          
*                                                                       
*     The program reflects the maximum through PBAR to PSTAR, and
*     evaluates the function at PSTAR (to give HSTAR).
*
      DO 270 I=1,NOP                                                    
  270 PSTAR(I)=A*(PBAR(I)-G(IMAX,I))+PBAR(I)                            
      CALL GIVEF (PSTAR,HSTAR,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*                                                                       
*     The next 5 statements test whether a progress report is
*     required and, if so, provide one.  This procedure occurs
*     frequently in the program.
*
      NEVAL=NEVAL + 1                                                    
*                                                                           
*     If the number of function evaluations to date exceeds 750
*     (=MAX), the program prints out parameter values provided that
*     IPRINT and JPRINT have been set at 1.
*
      IF (NEVAL.LE.MAX) GO TO 279                                            
  279 IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 280
      GO TO 300
  280 J=NEVAL/IPRINT                                                    
      K=NEVAL-J*IPRINT                                                  
      IF (K) 290,290,300
  290 WRITE(2,160) NEVAL,HSTAR,(PSTAR(J),J=1,NOP)                          
  300 IF (HSTAR.LT.HMIN) GO TO 310                                         
      GO TO 380                                                            
*                                                                       
*     If HSTAR is less than HMIN, PBAR is reflected through PSTAR
*     (to give PSTST) and the function is evaluated there (to
*     give HSTST).
*
  310 DO 320 I=1,NOP                                                    
  320 PSTST(I)=C*(PSTAR(I)-PBAR(I))+PSTAR(I)                            
      CALL GIVEF (PSTST,HSTST,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
*     If IPRINT=1 the program prints out the progress of the
*     iteration.  This is not normally required.
*
      NEVAL=NEVAL+1                                                     
      IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 330
      GO TO 350
  330 J=NEVAL/IPRINT                                                    
      K=NEVAL-J*IPRINT                                                  
      IF (K) 340,340,350
  340 WRITE(2,160) NEVAL,HSTST,(PSTST(J),J=1,NOP)                          
  350 IF (HSTST.LT.HMIN) GO TO 360                                           
  351 GO TO 560                                                           
*                                                                       
*     If HSTST is less than HMIN, the maximum point of the current
*     simplex is replaced by PSTST and HMAX is replaced by HSTAR,
*     then a test is performed.
*
  360 DO 370 I=1,NOP                                                   
  370 G(IMAX,I)=PSTST(I)                                                
      H(IMAX)=HSTST                                                     
      GO TO 580                                                         
*                                                                       
*     If HSTAR is not less than HMIN, the program tests is HSTAR
*     is greater than the function value at all vertices other
*     than the maximum one.
*
  380 DO 400 I=1,NP1                                                    
      IF (I.EQ.IMAX)  GO TO 400                                               
  390 IF (HSTAR.LT.H(I)) GO TO 560                                         
  400 CONTINUE                                                          
*                                                                       
*     If it is less than at least one of these vertices, the
*     maximum point of the current simplex is replaced by PSTAR and
*     HMAX by HSTAR.  A test is then performed.
*
*     If HSTAR is greater than all the function values excluding
*     the maximum, the program tests is HSTAR is greater than HMAX.
*     If it is not, the maximum point of whichever simplex is
*     now in store (depending on whether HSTAR is greater or less
*     than HMAX) is replaced by PSTAR and HMAX by HSTAR, and the
*     contracted point PSTST and the function value there (HSTST)
*     are calculated.
*
      IF (HSTAR.GT.HMAX)  GO TO 430                                           
  410 DO 420 I=1,NOP                                                    
  420 G(IMAX,I)=PSTAR(I)                                                
      HMAX=HSTAR                                                            
      H(IMAX)=HSTAR                                                         
  430 DO 440 I=1,NOP                                                    
  440 PSTST(I)=B*G(IMAX,I)+(1.0-B)*PBAR(I)                              
      CALL GIVEF (PSTST,HSTST,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
      NEVAL =NEVAL+1                                                    
      IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 450
      GO TO 470
  450 J=NEVAL/IPRINT                                                    
      K=NEVAL-J*IPRINT                                                  
      IF (K) 460,460,470
  460 WRITE(2,160) NEVAL,HSTST,(PSTST(J),J=1,NOP)                          
  470 IF (HSTST.GT.HMAX) GO TO 500                                          
*                                                                       
*     If HSTST is less than HMAX, the maximum point is replaced by
*     PSTST and HMAX by HSTST.  A test is then applied.
*
  480 DO 490 I=1,NOP                                                    
  490 G(IMAX,I)=PSTST(I)                                                
      H(IMAX)=HSTST                                                     
      GO TO 580                                                         
*                                                                       
*     If HSTST is not less than HMAX, each point in the current 
*     simplex is replaced by a point midway between its current
*     position and the position of the minimum point of the
*     current simplex.  The function is evaluated at each new
*     vertex and the test performed.
*
  500 DO 510 I=1,NP1                                                    
      DO 510 J=1,NOP                                                    
  510 G(I,J)=(G(I,J)+G(IMIN,J))/2.0                                     
      DO 550 I=1,NP1                                                    
      DO 520 J=1,NOP                                                    
  520 F(J)=G(I,J)                                                       
      CALL GIVEF (F,H(I),VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
      NEVAL=NEVAL +1                                                    
      IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 530
      GO TO 550
  530 J=NEVAL/IPRINT                                                    
      K=NEVAL-J*IPRINT                                                  
      IF (K) 540,540,550
  540 WRITE(2,160) NEVAL, H(I),(F(J),J=1,NOP)                             
  550 CONTINUE                                                          
      GO TO 580                                                         
  560 DO 570 I=1,NOP                                                    
  570 G(IMAX,I)=PSTAR(I)                                                
      H(IMAX)=HSTAR                                                     
*                                                                       
*     If LOOP=NLOOP, tests for convergence begin.  Otherwise
*     computation goes back to the beginning of the basic loop.
*
  580 IF (LOOP.EQ.NLOOP) GO TO 590                                          
  581 GO TO 180                                                             
*
*                                                                       
*   Tests for Convergence -
*
*     The mean and standard deviation of the function values of the
*     current simplex are now calculated.
*
  590 HSTD=0.0                                                             
      HMEAN=0.0                                                             
      DO 600 I=1,NP1                                                    
      HSTD=HSTD+H(I)*H(I)                                               
  600 HMEAN=HMEAN+H(I)                                                  
      HMEAN=HMEAN/FLOAT(NP1)                                                 
      HSTD=(HSTD-FLOAT(NP1)*HMEAN*HMEAN)/FLOAT(NP1)                          
      HSTD=SQRT(ABS(HSTD))                                                   
*                                                                       
*     The parameter values (F) at the centroid of the current
*     simplex and the function value there (FUNC) are now
*     calculated.
*
      DO 620 I=1,NOP                                                    
      F(I)=0.0                                                          
      DO 610 J=1,NP1                                                    
  610 F(I)=F(I)+G(J,I)                                                  
      F(I)=F(I)/FLOAT(NP1)                                                   
  620 CONTINUE                                                         
      CALL GIVEF(F,FUNC,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
      NEVAL=NEVAL+1
*
*     If KWT=1, the program calculates a biweighted least squares
*     value.  This process begins once the number of iterations
*	  exceeds a predetermined value (=40).
*
      IF (NEVAL.GT.40) KWT=1
*
*     If the number of evaluations has exceeded the value of MAX
*     set (=750), the convergence process is judged not to have
*     succeeded in this case.  If so, this particular run of Loop
*     855 is assumed to have yielded a 'no result'.  A counting
*     variable (MFAIL) is given a value of 1, to be subtracted from
*     the value of MAXJB later, and the next run through the loop
*     begins.
*
      IF (NEVAL.GT.MAX) GO TO 630
      GO TO 700 
*
*     MFAIL is a parameter which counts the number of times a series
*     of iterations failed to converge on a minimum within the set
*     maximum number of iterations.
*	                                                            
  630 MFAIL = MFAIL + 1
      IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 640
      GO TO 845
  640 WRITE(2,650) MAX                                                     
  650 FORMAT(40H NUMBER OF FUNCTION EVALUATIONS EXCEEDS ,I4)            
      WRITE(2,660) HSTD                                                    
  660 FORMAT(51H STANDARD ERROR OF FUNCTION VALUES OF LAST SIMPLEX ,E13.  
     &6)
      WRITE(2,670) (F(I),I=1,NOP)                                          
  670 FORMAT(28H  CENTROID OF LAST SIMPLEX  ,8E13.5,     (/28X,8E13.5)) 
      WRITE(2,680) FUNC                                                    
  680 FORMAT(31H  FUNCTION VALUE AT CENTROID   ,E13.6)                  
      GO TO 845  
*
*
  700 IF (HSTD.LT.STOPC)  GO TO 720
*
*     If the standard deviation calculated above is not less than
*     the criterion set (STOPC), IFLAG and LOOP are set to zero
*     and the basic loop begins again.
*
      IFLAG=0                                                           
      LOOP=0                                                            
      GO TO 180                                                         
  720 KWT=1
      IF ((IPRINT.EQ.1) .AND. (JPRINT.EQ.1)) GO TO 730
      GO TO 750
  730 WRITE(2,740)                                                         
  740 FORMAT(2H */33H  INITIAL EVIDENCE OF CONVERGENCE)                 
      WRITE(2,745) (F(I),I=1,NOP)  
  745 FORMAT(28H  CENTROID OF LAST SIMPLEX  ,8E13.5,    (/28X,8E13.5))                                        
      WRITE(2,746) FUNC 
  746 FORMAT(31H  FUNCTION VALUE AT CENTROID   ,E13.6)                                                   
*                                                                       
*     If the standard deviation is less than the stopping 
*     criterion, IFLAG is set =0 if there was no evidence of
*     convergence on the last test and =1 if there was evidence
*     of convergence.
*
  750 IF (IFLAG) 760,760,770                                                 
*                                                                       
*     If IFLAG=0, reset IFLAG=1,  the mean of the function
*     values of the current simplex are saved (as SAVEMN), and  
*     computation goes back to the beginning of the basic loop.
*
  760 IFLAG=1                                                           
      SAVEMN=HMEAN                                                           
      LOOP=0                                                            
      GO TO 180                                                         
*                                                                       
*     If IFLAG=1, the program tests if the change in the mean is
*     less than the stopping criterion.  If it is, the process is
*     said to have converged.  If not, IFLAG and LOOP are both set
*     at zero and computation reverts to the start of the
*     basic loop.
*
  770 IF (HMEAN.EQ.0) GO TO 790
      TEST=SAVEMN/HMEAN
      IF (TEST.GT.0.99995 .AND. TEST.LT.1.00005)  GO TO 790                                           
  780 IFLAG=0                                                           
      LOOP=0                                                            
      GO TO 180                                                         
*
*     If a calculated value of F(4) [=dmax] is greater than the
*     upper limit (DLIM) set by the program, F(4) is set at
*     DLIM, STEP(4) is set at zero, and the search goes
*     back to the beginning.
*
  790 IF (DLIM.GE.F(4)) GO TO 797
      F(4)=DLIM
      STEP(4)=0.
      GO TO 180
*
*     If JPRINT=1 the program prints out the results of each successful
*     convergence on a minimum.
*
  797 IF (JPRINT) 850,850,800
  800 WRITE(2,810) NEVAL                                                   
  810 FORMAT(5(/),36H PROCESS CONVERGES ON MINIMUM AFTER ,I4,21H FUNCTIO 
     &N EVALUATIONS///)                                                 
      WRITE(2,820) (F(I),I=1,NOP)                                          
  820 FORMAT(14H MINIMUM AT   ,4(1X,E13.6))                              
      WRITE(2,830) FUNC                                                    
  830 FORMAT(//26H MINIMUM FUNCTION VALUE   ,E13.6)                     
      WRITE(2,840)                                                         
  840 FORMAT(///16H END  OF  SEARCH/1X,15(1H*)) 
*
*
*	  
  850 CONTINUE   
*                                                         
*
*     MTEST is set at 1 to flag that convergence has occurred.  This
*     is carried into Subroutine GIVEF to trigger computation of
*     MSFAIL where computation of S cannot occur.
*
      MTEST=1
*
*     Program execution returns to the subroutine to yield final pass values 
*     of F(1), F(2) and F(3).
*
      CALL GIVEF(F,FUNC,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
      LPRINT=1
      KPRINT=0                                                            
*                                                                       
*     
*     The estimated 'best fit' values of the parameters from the current
*     pass through Loop 855 are now computed, tabulated and printed out.
*
*     A density estimate (DEN) is calculated from D2L=F(3) by
*     correcting units to no./ha, and dividing by 2L in the case of
*     line transect data, and by 2Vt in the case of fixed-point data.
*
 1428 IF (KM) 1436,1436,1435
 1435 DEN(JB)=(1.E6*F(3))/(2.*DIST)
      GO TO 1437
 1436 IF (IFX) 1440,1440,1441
 1440 DEN(JB)=(1.E4*F(3))/(2.*DIST)
      GO TO 1437
 1441 DEN(JB)=(1.E4*F(3))/(2.*IV*IT)
*
*     Values of the other parameter estimates are obtained by redefining
*     the estimates of F(1) and F(2), thus:
*
 1437 COEFF1(JB)=F(1)
      COEFF2(JB)=F(2)
      COEFF3(JB)=S
*
*
*     Running totals of DEN, COEFF1 and COEFF2 are now made to
*     enable calculation of mean values after the loop ends.
*
      TDEN = TDEN + DEN(JB)
      TCOEFF1 = TCOEFF1 + COEFF1(JB)
      TCOEFF2 = TCOEFF2 + COEFF2(JB)
      TCOEFF3 = TCOEFF3 + COEFF3(JB)
* 
*                                  
*     Loop 855 now ends, returning calculations to Line 90 until the maximum
*     preset number of bootstraps (MAXJB) value is reached.  LOOP is reset 
*     to zero to enable a new series of iterations in the basic loop to
*     begin again, as are G(I,J) values.
*
*
  845 LOOP=0
      IFLAG=0
      KPRINT=0
*
      DO 853 I=1,NP1
      DO 853 J=1,NOP
  853 G(I,J) = 0.0
*
*     F and STEP values are reset to their original values also.
*
      DO 847 IE=1,NOP
      F(IE)=FT(IE)
      STEP(IE)=STEPT(IE)
  847 CONTINUE
*
*     	    	     
  855 CONTINUE  
*
*
*     Following completion of the runs through Loop 855, the means and
*     standard errors of each of the parameters are now calculated, based
*     on the values of the three arrays DEN(JB), COEFF1(JB) and COEFF2(JB).
*
*     The overall means ESTDEN, COEFFNT1 and COEFFNT2 are calculated first.
*     NUMEST is the number of parameter estimations made.
*
      NUMEST=MAXJB-MFAIL
      IF (NUMEST.GE.1) GO TO 1446
      NUMEST=1
 1446 ESTDEN=TDEN/NUMEST
      COEFFNT1=TCOEFF1/NUMEST
      COEFFNT2=TCOEFF2/NUMEST
*
*
*     Should MSFAIL be identical to NUMEST, 1 is added to NUMEST to prevent
*     division by zero at the next step.
*
      IF (NUMEST.EQ.MSFAIL) NUMEST=MSFAIL + 1
      COEFFNT3=TCOEFF3/(NUMEST-MSFAIL)
      IF (NUMEST.EQ.(MSFAIL+1)) NUMEST=NUMEST - 1
*
*
*     The next step is to calculate the standard errors of each 
*     parameter, provided that the number of analyses exceeds 1.
*     Each is the standard deviation of the parameter estimates.
*     If NUMEST is 0 or 1, standard error calculation is bypassed.
*
      IF (NUMEST.LE.1.) GO TO 1484
      DSUM=0.0
      DO 1450 JB=1,MAXJB	  
      IF ((DEN(JB)).EQ.0.) GO TO 1450
      DENDIF=(DEN(JB)-ESTDEN)**2.
      DSUM=DSUM+DENDIF
 1450 CONTINUE
      SDEN=SQRT(DSUM/(NUMEST-1))
*
      CF1SUM=0.0
      DO 1451 JB=1,MAXJB	  
      IF ((COEFF1(JB)).EQ.0.) GO TO 1451
      CF1DIF=(COEFF1(JB)-COEFFNT1)**2.
      CF1SUM=CF1SUM+CF1DIF
 1451 CONTINUE
      SCF1=SQRT(CF1SUM/(NUMEST-1))
*
      CF2SUM=0.0
      DO 1452 JB=1,MAXJB	  
      IF ((COEFF2(JB)).EQ.0.) GO TO 1452
      CF2DIF=(COEFF2(JB)-COEFFNT2)**2.
      CF2SUM=CF2SUM+CF2DIF
 1452 CONTINUE
      SCF2=SQRT(CF2SUM/(NUMEST-1))
      CF3SUM=0.0
      DO 1453 JB=1,MAXJB
      IF (COEFF3(JB).EQ.0) GO TO 1453
      CF3DIF=(COEFF3(JB)-COEFFNT3)**2.
      CF3SUM=CF3SUM+CF3DIF
 1453 CONTINUE
      IF ((NUMEST-MSFAIL-1).LE.0) GO TO 1484
      SCF3=SQRT(CF3SUM/(NUMEST-MSFAIL-1))
*
*
*     The products of the program are now output.
*
*
*     The estimated topographical cover (TCOV) is first, followed
*     by the number of actual parameter estimations made (NUMEST).
*
 1484 WRITE(2,1471) TCOV
 1471 FORMAT(/33H Estimated Topographical Cover = ,F7.6) 
*
      WRITE(2,1472) NUMEST
 1472 FORMAT(/34H Number of Parameter Estimations =,I4)
* 
*
*     Now follow general headings for the results table.
*
*
      WRITE(2,835)                                                             
  835 FORMAT(//24H Calculated values were:/)                                
      KPRINT=1   
*
      WRITE(2,1485)
 1485 FORMAT(//28H ESTIMATED PARAMETER VALUES://)
      WRITE(2,1486)
 1486 FORMAT(1H ,77(1Hx)/2H x,4(18X,1Hx))
      WRITE(2,1487)
 1487 FORMAT(15H x    PARAMETER,5X,1Hx,6X,5HVALUE,7X,1Hx,
     &18H  STANDARD ERROR  ,1Hx,7X,4HUNIT,7X,1Hx)
      WRITE(2,1488)
 1488 FORMAT(2H x,4(18X,1Hx)/1H ,77(1Hx)/2H x,4(18X,1Hx)/12H x ESTIMATED
     &,8X,1Hx,3(18X,1Hx))
*
*     Density estimates are printed.
*
 1479 IF (KM) 1489,1489,1491
 1489 IF (MAXJB-1) 2000,2000,2002
 2000 WRITE(2,2001) ESTDEN
 2001 FORMAT(15H x DENSITY  (D),5X,1Hx,4X,F10.3,4X,
     &40Hx (indeterminate)  x  indivs./hectare  x)
      GO TO 1493
 2002 WRITE(2,1490) ESTDEN,SDEN
 1490 FORMAT(15H x DENSITY  (D),5X,1Hx,4X,F10.3,4X,1Hx,4X,F10.3,4X,
     &20Hx  indivs/hectare  x)
      GO TO 1493
*
 1491 IF (MAXJB-1) 2003,2003,2004
 2003 WRITE(2,1492) ESTDEN
 1492 FORMAT(15H x DENSITY  (D),5X,1Hx,4X,F10.2,4X,
     &39Hx (indeterminate)  x   indivs./sq.km. x)
      GO TO 1493
 2004 WRITE(2,2005) ESTDEN,SDEN
 2005 FORMAT(15H x DENSITY  (D),5X,1Hx,4X,F10.3,4X,1Hx,4X,F10.3,4X,
     &20Hx   indivs./sq.km. x)
*
*     The conspicuousness coefficient is next printed.
*
 1493 WRITE(2,1494)
 1494 FORMAT(2H x,4(18X,1Hx)/1H ,77(1Hx)/2H x,4(18X,1Hx)/
     &21H x CONSPICUOUSNESS  x,3(18X,1Hx))
      IF (MAXJB-1) 2006,2006,2007
 2006 WRITE(2,1495) COEFFNT1
 1495 FORMAT(21H x COEFFICIENT  (a) x,4X,F10.4,4X,
     &39Hx (indeterminate)  x      metres      x)
      GO TO 1501
 2007 WRITE(2,2008) COEFFNT1,SCF1
 2008 FORMAT(21H x COEFFICIENT  (a) x,4X,F10.3,4X,1Hx,4X,F10.3,4X,1Hx,
     &6X,6Hmetres,6X,1Hx)
*
*     The second coefficient is either a cover proportion or a sound
*     attenuation coefficient, decided by the values of KDT.
*
 1501 IF (KDT) 1496,1496,1500
 1496 WRITE(2,1497)
 1497 FORMAT(2H x,4(18X,1Hx)/1H ,77(1Hx)/2H x,4(18X,1Hx)/14H x COVER   
     &   ,6X,1Hx,3(18X,1Hx))
      IF (MAXJB-1) 1509,1509,1512
 1509 WRITE(2,1511) COEFFNT2
 1511 FORMAT(21H x PROPORTION  (c)  x,4X,F10.4,4X,
     &39Hx (indeterminate)  x                  x)
      GO TO 2013
 1512 WRITE(2,1513) COEFFNT2,SCF2
 1513 FORMAT(21H x PROPORTION  (c)  x,4X,F10.4,4X,1Hx,4X,F10.4,4X,1Hx,
     &4X,9H         ,5X,1Hx)
      GO TO 2013
*
 1500 WRITE(2,1502)
 1502 FORMAT(2H x,4(18X,1Hx)/1H ,77(1Hx)/2H x,4(18X,1Hx)/
     &14H x ATTENUATION,6X,1Hx,3(18X,1Hx))
      IF (MAXJB-1) 2014,2014,2016
 2014 WRITE(2,2015) COEFFNT2
 2015 FORMAT(21H x COEFFICIENT  (b) x,4X,F10.4,4X,
     &39Hx (indeterminate)  x    per  metre    x)
      GO TO 2013
 2016 WRITE(2,1504) COEFFNT2,SCF2
 1504 FORMAT(21H x COEFFICIENT  (b) x,4X,F10.4,4X,1Hx,4X,F10.4,4X,1Hx,
     &4X,9Hper metre,5X,1Hx)
*
*     The table is now ruled off.
*
 2013 WRITE(2,1499)
 1499 FORMAT(2H x,4(18X,1Hx)/1H ,77(1Hx)//)
*
      WRITE(2,2017) COEFFNT3,SCF3
 2017 FORMAT(/X,31HDetectability Coefficient (S) =,F8.2,6H, SE =,F6.2)
*
*
*     Key model estimates are now output if ISHOW was originally set
*     at 1.  The best estimates of F(1), F(2) and F(3) must be entered
*     in the subroutine GIVEF to do this.  The totals in the class
*     intervals are those from the initial data, derived from the
*     saved array variable VALT().  The output is then produced
*     within the subroutine.
*
*
      F(1) = COEFFNT1
      F(2) = COEFFNT2
      DO 1530 JV=1,80
 1530 VAL(JV)=VALT(JV)
 1531 IF(KM) 1517,1517,1516
 1516 F(3)=(ESTDEN*2.*DIST*PD*(SNS/2))/1.E6
      GO TO 1520
 1517 IF (IFX) 1518,1518,1519
 1518 F(3)=(ESTDEN*2.*DIST*PD*(SNS/2))/1.E4
      GO TO 1520
 1519 F(3)=ESTDEN*2.*IV*IT*PD/1.E4	       	  	  
 1520 CALL GIVEF(F,FUNC,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
*	  
 1525 CONTINUE
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
*     STOP
      RETURN
*
*
  283 WRITE(6,285)INFILE,IOS
  285 FORMAT(' Error opening ',A40,' - IOS = ',I6)
*     STOP
      RETURN
  284 WRITE(6,285)OUTFILE,IOS
*     STOP
      END SUBROUTINE CALCULATE_DENSITY
*
*     *************************************************************
*
*
*     This version of Subroutine GIVEF will handle ground survey
*     data [and also aerial survey data for which there is complete
*     visibility ahead of the aircraft].  The subroutine
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
*
      SUBROUTINE GIVEF(F,FUNC,VAL,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH,S,SNS,
     & IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT,LTMAX,LTMIN,
     & NUMA,NUMO,NVALS,NS,MSFAIL,MTEST)
*
*     Double precision for all real numbers is set, together with
*     common values, dimensions and symbols for key variables.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IMV,IRY,ISHOW
      INTEGER KDT,KWT,LTMIN,LTMAX,MAX,MAXJB,NS,NUMA,NUMEST,NUMO,NVALS 
      DIMENSION R(5000),BSTR(5000),Y(5000),BSTY(5000),
     & DEN(5000),COEFF1(5000),COEFF2(5000)
      DIMENSION G(5,4),STEP(4),STEPT(4),F(4),VAL(80)
      DIMENSION ROUT(80),CALCNR(80),OBSDNR(80) 
      DIMENSION YOUT(80),CALCNY(80),OBSDNY(80)           
*
*     COMMON / common1 /VAL(80),H(4),PBAR(4),PSTAR(4),PSTST(4)
*     COMMON / common2 /HSTAR,HSTST,CLINT,PD,PS,R3S,STT,TCOV,THH,VGH
*    &,S,SNS
*     COMMON / common3 /KM,IFX,IMV,IRY,ISHOW,IT,IV,KDT,KPRINT,KWT,LPRINT
*     COMMON / common4 /LTMAX,LTMIN,MAX,NUMA,NUMO,NVALS,NS,MSFAIL,MTEST
*
*     TOT is the progressive value of the sum of squares of the
*     differences between observed and calculated values.  It is
*     initially set at zero.
*
      TOT=0.                             
      WTOT=0.                              
*
*     The program now sets upper and lower limits to the
*     parameter values supplied, usually by effectively giving
*     FUNC (through HTOT) a very high value if a parameter falls
*     outside predetermined limits.
*
      P=F(1)                                                            
*
*     HTOT is a variable used to deflect the search away from
*     highly improbable values of the parameters.
*
      HTOT=0.
      DMAX=F(4)
*
*     HTOT is set high if DMAX < 1
*
      IF (DMAX.GE.1.) GO TO 5002
*
*     Initially, HTOT is set at 1 000 000
*
      HTOT=1.E+6
 5002 Q=F(2)
*
*     If Q is negative, the program sets IMV=1 and so uses the
*     median value of d in an interval as the basis of
*     computations, in order to avoid logarithms of negative
*     values appearing in Approximation 1 below.
*
      IF (Q) 5003,5005,5005
*
 5003 IMV=1
*
*
*     F(3) is renamed D2L for its run through Subroutine GIVEF.
*
 5005 D2L=F(3)
*
*     HTOT is set to a higher value if D2L<0.0001
*
      IF (D2L-.0001) 2020,5400,5400
 2020 HTOT=1.E+6 + HTOT
*
*     HTOT is set high if a>400
*
 5400 IF (P-400.) 5800,5800,5500
 5500 HTOT=1.E+6 + HTOT
*
*     HTOT is set high if a<0.0099
*
 5800 IF (P-.0099) 5900,5900,5950
*
 5900 HTOT=ABS((P-2.)*1.E06)+HTOT
*
*     RMAX is set at zero if DMAX is equal to or less than THH
*
 5950 IF (DMAX.GT.THH) GO TO 6000
      RMAX=0.0
      GO TO 6001
*
*     RMAX - the maximum horizontal recognition distance - is
*     computed from the direct line distance and height difference.
*
 6000 RMAX=DSQRT(DMAX*DMAX-THH*THH)
 6001 IF (Q) 6003,6002,6002                                             
*
*     An upper limit of 0.4 is set for the attenuation coefficient
*
 6002 IF (Q.LT.0.4) GO TO 6006                                            
      HTOT=1.E+6 + HTOT
      GO TO 6006                                                        
*
*     A lower limit of -2./DMAX is set to the attenuation coeffnt.
*
 6003 QMIN=-2./DMAX                                                     
      IF (Q.GT.QMIN) GO TO 6006                                           
      Q=QMIN                                                            
*
*     The theoretical probability of detecting an individual at
*     the maximum recognition distance [PRMAX=Pr(rmax)] is now
*     calculated.
*
*     PA is the square of the conspicuousness coefficient.
*
*
 6006 PA=P*P
*
*     Visibility and audibility will be affected by topographical
*     features in habitats where the ground is not level. The
*     total probability of visibility at DMAX (VISMAX) will be the product of
*     the probability VEGMAX that the animal is unobscured by vegetation
*     at distance DMAX, and the probability TOPMAX that it is unobscured
*     by topography.  VEGMAX will be a function of dMAX (VEGMAX=(1-Q)**DMAX),
*     while TOPMAX is approximated by the function  
*     TOPMAX=(1-T)**(LTMAX-DMAX)=(EXP(ln(0.001)/(LTMAX-LTMIN))**(LTMAX-DMAX),
*     and VISMAX=VEGMAX*TOPMAX.  A first step is to calculate TOPMAX, 
*     provided that either LTMIN is not very large or RMAX is currently
*     less than LTMIN.
*
      IF ((LTMIN.LT.999.).AND.(RMAX.GT.LTMIN))  GO TO 6010
      TOPMAX=1.
      GO TO 6011
 6010 TOPMAX=(EXP(-6.9078/(LTMAX-LTMIN)))**(DMAX-LTMIN)
*
*     For observing situations (e.g. aerial survey) where there is 'ground'
*     cover for only the first part of the direct-line distance
*     d between animal and observer (indicated by the vegetation
*     height (VGH) exceeding zero), this cover will obscure
*     some animals.  The proportion (VISMAX) visible at DMAX will
*     be a function of d (VISMAX=(1-c)**DMAX). The distance obscured
*     (DVG) will have a value equal to (cover height) x (distance d)/
*     (observer-animal height difference).
*
 6011 IF (KDT) 5959,5959,5961
* 
*     If the calculated cover proportion is less than zero, the 
*     program adds an arbitrary 100 to the accumulating total 
*     to deflect the search away from such unreal values.
*
 5959 IF (Q.GE.0.) GO TO 5960
      HTOT=1.E+2 + HTOT
*
*     If visual data are supplied, and there is vegetation close
*     to the ground only, and the animals are within that
*     vegetation, the cover proportion correction applies only
*     to the proportion of the direct-line distance potentially
*     obscured by that ground vegetation.
*
 5960 IF (KDT) 5958,5958,5961
 5958 IF (VGH) 5963,5963,5964
 5964 DVG=VGH*DMAX/THH
      IF (DMAX.LE.DVG) GO TO 5963
      VEGMAX=(1.-Q)**DVG
      VISMAX=VEGMAX*TOPMAX
      GO TO 5965
 5963 VEGMAX=(1.-Q)**DMAX
      VISMAX=VEGMAX*TOPMAX
 5965 DDSM=DMAX*DMAX
      PRMAX=PA*VISMAX/DDSM
*
*     For visual data collected where there is cover, IMV is put
*     =1 to avoid the possibility of negative values being taken
*     to logarithms later in the program.
*
      IMV=1
      GO TO 5962
*
 5961 QDMAX=Q*DMAX                                                      
*
*     To prevent overflow during computations, and upper limit of
*     70 and a lower limit of -68 are set to DMAX.
*
      IF (QDMAX.LT.70.) GO TO 6007                                       
      QDMAX=70.                                                         
      GO TO 6008                                                        
 6007 IF ((-QDMAX).LT.68.) GO TO 6008                                    
      QDMAX=-68.                                                        
 6008 SSMAX=DEXP(QDMAX)                                                  
      DDSM=DMAX*DMAX                                                    
      PAM=PA/DDSM                                                       
      PRMAX=PAM/SSMAX                                                   
*
*     The control variable L10 determines the number of classes
*     for which an expected value is calculated.  It is set equal
*     to RMAX/CLINT for all situations except those where a strip
*     below an aircraft is hidden from the observer.
*     0.5 is added to remove errors due to 'chopping'.
*
 5962 L10=(RMAX/CLINT) + 0.5
      TOTE=0.
*
*     TR is the highest r value in the current frequency class;
*     because the first class computed is that furthest
*     from the observer or from the transect line - TR is initially
*     set at the class width (CLINT) times the number of classes.
*     This is adjusted upwards by STT if this starting value is 
*     greater than zero.
*
      TR=CLINT*FLOAT(L10) + STT
*
*
*     To compute 'expected' values within each class, each class
*     is subdivided into subclasses, each of width 800/L10, so that
*     CINT is small and L10 x L20=800.  Each subclass
*     corresponds to an observing arc of width CINT ('delta-r') 
*     that sweeps forwards ahead of the observer.  L20 is the
*     number of classes in the inner loop.
*     0.5 is added to remove errors due to 'chopping'.
*
      L20=(800/L10) + 0.5
*
*     CINT ('delta-r') is the width of each subclass; it is set 
*     at the class width divided by the number of subclasses used.
*
      CINT=CLINT/FLOAT(L20)
      HCINT=CINT/2.
*     
*     With perpendicular distance data, RMAX is altered slightly
*     to become an exact multiple of CINT.
*
      IF (IRY) 7007,7006,7007
 7006 IERMAX=(F(4))/CINT
      ERMAX=IERMAX*CINT
*
*     In the event that PRMAX > 1., HTOT is set to a high value.
*
 7007 IF (PRMAX.LT.1.) GO TO 7009                                        
      HTOT=PRMAX*1.E06+HTOT                                                  
*
*     QR is the proportion of the number of animals originally in
*     a subclass still present after animals detected in arcs
*     further out have been disregarded.  It is initially unity.
*
 7009 QR=1.
*
*     JRLOW - a control variable used to control the printout
*     of RLOW - is set initially at zero.
*
      JRLOW=0
*
*     EXPD is the cumulative expected number within a class; it is
*     set at zero before computations begin; so is EXPDV.
*
      EXPD=0.                                                          
      EXPDV=0.
      S=0.
*
*     Loop 10, which calculates expected values and compares them
*     with observed values, now begins .....
*
*
      DO 10 JJ=1,L10
*
*     MD is the hth class in the series from MD=80 to MD=1 (or 2,
*     in situations where the most central belt is obscured).
*     0.5 is added to enable rounding to the nearest integer
*     despite chopping by the program.)  To calculate it, TR is
*     reduced by the starting value (STT) if this exceeds zero.
*
      MD=INT((TR-STT)/CLINT + 0.5)
*
*     OBSD is the observed value for the hth arc, as supplied
*     in the input to the program.
*
      OBSD=VAL(MD)
*
*     The observed frequency value in the hth class is corrected to
*     allow for animals that overtake the observer from behind.
*
      OBSD=((NUMA+NUMO)*OBSD)/NUMA
*
*     Where the data are radial distance or fixed point data, EXPD 
*     accumulates the expected frequency values within a class; it
*     is set at zero before computations for the class begin.
*
*     Where the data are perpendicular distance data (IRY=1), TEXPD
*     accumulates the expected frequency values within each class;
*     it is set at zero before computations for the class begin,
*     while EXPD is allowed to accumulate.
*
      IF (IRY) 7021,7022,7021                                          
 7021 EXPD=0.                                                          
      GO TO 7014
 7022 TEXPD=0.
 7014 IF (TOTE.LT.0.) GO TO 7017
      GO TO 7023
 7017 TOTE=0.
*
*     WL is the lower boundary of the current class.
*
 7023 WL=TR-CLINT
*
*     If the entire class lies above both the highest calculated
*     and highest observed recognition distances, the expected 
*     value (EXPDV) is set at zero and the calculations in Loop
*     20 are bypassed.
*
      IF (WL.GT.RMAX .AND. WL.GT.ERMAX) GO TO 7015
      GO TO 7020
 7015 EXPDV=0.
      GO TO 8100                                                        
*
*     DNR is the central r value in an observing arc.
*     It is initially set half of CINT below TR
*
 7020 DNR=TR-HCINT
*
*     DNRH is the highest r value in an observing arc.  It is
*     initially equal to TR.
*
      DNRH=DNR+HCINT
*
*     DH is the direct-line distance to the outer edge of the 
*     observing arc.
*
      DH=DSQRT(THH*THH+DNRH*DNRH)
*
*     Two alternative ways of calculating the probability of
*     detection are provided.  If the method is based on the area
*     under the probability curve (IMV=0), the area under
*     the curve from y=DH to y=(infinity) is calculated and
*     expressed as YYH.  If IMV=1, this calculation is bypassed.
*
      IF (IMV) 7024,7024,8000
*
*     This method of computation uses one of two approximations to
*     the exponential integral (APREXI).  If bd < 1, Approximation
*     2 is used; if bd is equal to or greater than 1, Approximation
*     1 is used.
*
 7024 ZH=Q*DH
*
*     To avoid overflow during computations, ZH is set at 68 if 
*     ZH > 68.
*
      IF (ZH.LT.68.) GO TO 7025
      ZH=68.
      GO TO 7026
*
*     ZH is set at -70 if ZH < 70.
*
 7025 IF ((-ZH).LT.70.) GO TO 7026
      ZH=-70.
*
*     If ZH is zero, VHOWH (=VH/WH) is set at 0, and a few lines of
*     calculations are bypassed.
*
 7026 IF (ZH) 7030,7027,7030
 7027 VHOWH=0.
      GO TO 7040
*
*     The choice of approximation depends on the value of bd.
*
 7030 IF (ZH.GE.1.0) GO TO 7035
*    
*     Approximation 2 is used if bd<1.
*
      APREXI=-.577216+.999992*ZH-.249910*ZH**2+.055200*ZH**3
     $-.009760*ZH**4+.0010792*ZH**5-LOG(ZH)
      GO TO 7040
*
*     Approximation 1 is used if bd=1 or bd>1.
*
 7035 VH=ZH*ZH+2.334733*ZH+0.250621
      WH=ZH*ZH+3.330457*ZH+1.681534
      VHOWH=VH/WH
 7040 SSH=DEXP(ZH)
      YH=PA/(DH*SSH)
*
*     YYH is the area under the detectability curve from
*     the current DH value to infinity.
*
      IF (ZH.GT.0. .AND. ZH.LT.1.0) GO TO 7050
      YYH=YH*(1.-VHOWH)
      GO TO 8000
 7050 YYH=YH-PA*Q*APREXI
*
*     Loop 20, which calculates the expected number in each 
*     arc within the class, and adds them together to produce
*     a progressive total (EXPD and TEXPD), now begins .....
*
*
 8000 DO 20 JL=1,L20
*
*     DNRL is the lowest r value in each arc; it is calculated
*     simply by subtracting the arc width from the highest r
*     value in the arc.
*
      DNRL=DNRH-CINT
*   
*     If, in the last class evaluated, DNRL comes to a negative 
*     value, it is then set at zero.
*
      IF (DNRL) 8001,8002,8002
 8001 DNRL=0.
*
*     DL is the direct-line distance from the observer to the
*     inner edge of the current observing arc.
*
 8002 DL=DSQRT(THH*THH+DNRL*DNRL)
*
*     DINT is the difference between DL and DH, and thus the width
*     of the area under the probability density curve between DH
*     and DL.
*
      DINT=DH-DL
*
*     If IMV has been set at zero, calculation of the probability
*     of detection now follows, using the same approximation to the
*     exponential integral as previously.  If IMV has been set at 
*     1, computation moves to address 8020.
*
      IF (IMV) 8004,8004,8020
 8004 ZL=Q*DL
      IF (ZL.LT.68.) GO TO 8008
      ZL=68.
      GO TO 8009
 8008 IF ((-ZL).LT.70.) GO TO 8009
      ZL=-70.
 8009 IF (ZL) 8011,8010,8011
 8010 VLOWM=0.
      GO TO 8013
*
*     The choice of approximation depends on the value of bd.
*
 8011 IF (ZL.GE.1.0) GO TO 8012
*    
*     Approximation 2 is used if bd<1.
*
      APREXI=-.577216+.999992*ZL-.249910*ZL**2+.055200*ZL**3
     &-.009760*ZL**4+.0010792*ZL**5-LOG(ZL)
      GO TO 8013
*
*     Approximation 1 is used if bd=1 or bd>1.
*
 8012 VL=ZL*ZL+2.334733*ZL+0.250621
      WM=ZL*ZL+3.330457*ZL+1.681534
      VLOWM=VL/WM
 8013 SSL=DEXP(ZL)
      YL=PA/(DL*SSL)
*
*     YYL is the area under the detectability curve from
*     the current DL value to infinity.
*
      IF (ZL.GT.0. .AND. ZL.LT.1.0) GO TO 8014
      YYL=YL*(1.-VLOWM)
      GO TO 8016
 8014 YYL=YL-PA*Q*APREXI
*
*     AUC is the area under the detectability curve between
*     d=DL and d=DH.
*
 8016 AUC=ABS(YYH-YYL)
*
*     The mean corrected probability [PR=P(r)] of detecting an
*     individual animal present between d=DL and d=DH is given by
*     dividing the area under the curve by the arc width, then
*     subtracting Pr(rmax).  It is the mean height of the curve.
*
      PR=AUC/DINT
      GO TO 7083
*
*     Where IMV was set at 1, computation of the probability 
*     is based on its median value in the arc rather than
*     the mean height of its density curve.
*
*     DD is the direct-line distance to the centre of the arc.
*
 8020 DD=DSQRT(THH*THH+DNR*DNR)
*
*     Visibility and audibility will be affected by topographical
*     features in habitats where the ground is not level. The
*     total probability of visibility at DD (VISDD) will be the product of
*     the probability VEGDD that the animal is unobscured by vegetation
*     at distance DD, and the probability TOPDD that it is unobscured
*     by topography.  VEGDD will be a function of d (VISDD=(1-Q)**DD),
*     while TOPDD is approximated by the function  
*     TOPDD=(1-TCOV)**(DD-LTMIN)=(EXP(ln(0.001)/(LTMAX-LTMIN)))**(DD-LTMIN).
*     and VISDD=VEGDD*TOPDD.  A first step is to calculate TOPDD, 
*     provided that LTMIN is not very large value or d is 
*     currently less than LTMIN.  [If LTMIN is set at 999 or higher, 
*     topography is assumed not to affect detectability, so TOPDD is set
*     at 1 and TCOV at zero.  If DD is less than LTMIN, topography is also
*     assumed not to affect detectability, so TOPDD is again set at 1.]
*
      IF (LTMIN.GE.999.) GO TO 8003
      IF (DD.GE.LTMIN) GO TO 8019
      TOPDD=1.
      GO TO 8021
 8003 TOPDD=1.
      TCOV=0.
      GO TO 8021
 8019 TOPDD=(EXP(-6.9078/(LTMAX-LTMIN)))**(DD-LTMIN)
      TCOV=1-(EXP(-6.9078/(LTMAX-LTMIN)))
*
*     For observing situations where there is 'ground'
*     cover for only the first part of the direct-line distance
*     d between animal and observer (indicated by the vegetation
*     height (VGH) exceeding zero), this cover will obscure
*     some animals.  The proportion visible at d (VISDD) will
*     be a function of d (VISDD=(1-Q)**DD). The distance obscured
*     (DVG) will have a value equal to (cover height) x (distance d)/
*     (observer-animal height difference).
*
 8021 IF (KDT) 8017,8017,8018
 8017 IF (VGH) 8026,8026,8025
 8025 DVG=VGH*DD/THH
      IF (DD.LE.DVG) GO TO 8026
      VEGDD=(1.-Q)**DVG
      VISDD=VEGDD*TOPDD
      GO TO 8027
 8026 VEGDD=(1.-Q)**DD
      VISDD=VEGDD*TOPDD
 8027 DDS=DD*DD
      PR=PA*VISDD/DDS
      GO TO 7083
*
 8018 QDD=Q*DD
      IF (QDD.LT.68.) GO TO 8022
      QDD=68.
      GO TO 8023
 8022 IF ((-QDD).LT.70.) GO TO 8023
      QDD=-70.
 8023 SS=DEXP(QDD)
      DDS=DD*DD
      PAD=PA/DDS
      PR=PAD/SS
*
*     PRC [=P(r)] is the height of the detectability curve at DD,
*     corrected by subtracting PRMAX.
*
 7083 PRC=PR-PRMAX
*
*     Because 0<PRC<1, an upper limit of 1 and a lower limit of
*     0 are set to the probability function.
*
      IF (PRC.LE.1.) GO TO 7084
      PRC=1.
*
 7084 IF (PRC.GE.0) GO TO 8033
      PRC=0
*
*     The probability (PRR=g(r)) of detection within an arc of width 
*     CINT is different from that in an arc of unit width, and is
*     equal to the product of the corrected probability and the
*     arc width.  
*
 8033 PRR=PRC*CINT
*
*     If the median r value in an arc is above both the expected
*     and observed highest r values, PRR is set at zero and QR at 1
*
      IF (DNR.GE.ERMAX .AND. DNR.GE.RMAX) GO TO 8034
*
*     If the median r value is above the maximum expected r value
*     but less than the maximum observed r value, PRR is allowed to
*     have a negative value.
*
      IF (DNR.LE.RMAX .OR. DNR.LE.ERMAX) GO TO 8036
 8034 PRR=0.
*
*     E [= g(r).P(r)] is a probability density function which
*     describes the probability that an animal is both present and
*     detected in an observing arc.  It is calculated by multiplying
*     together the probability of detection in the arc [PRC=g(r)]
*     and the probability [QR=P(r)] that an individual is still
*     available for detection in that arc. TOTE is the accumulating
*     total of the probability density function.
*
 8036 E=PRR*QR
      TOTE=TOTE+E
*
*     The subroutine now calculates the number of detections (ED)
*     expected for radial, fixed-point and perpendicular distance data.
*     For radial data, computation is based on the assumption 
*     that the arc of
*     radius r and width delta-r sweeps over a series of plots of
*     unit area.  The total number of detections expected in such
*     a plot will be [plot area] x [apparent density] x [probability
*     of detection].
*
*     With perpendicular distance data, ED is the expected number
*     detected at a distance r from the observer in the strips
*     CINT units wide at a perpendicular distance y=r from the
*     transect line.  Because there are one or two such strips, each
*     L times the area of the individual plot at distance r,
*     the expected no. = D x NS x LJ x Pd x CINT x g(r) x P(r) x Q(r).
*     NS must be converted to a floating-point number (SNS) first.
*
      SNS=FLOAT(NS)
      IF (IRY) 8043,8042,8043
 8042 IF (PD.GT.0.) GO TO 8047
      ED=D2L*CINT*E*(SNS/2.)
      GO TO 8050
 8047 ED=D2L*CINT*E*PD*(SNS/2.)
      GO TO 8050
*
*     With radial distance data, each observing arc sweeps out an
*     area L units long and 2r units wide, within which there are
*     2r/CINT x L individual plots at radial distance r.  Thus the
*     expected total number detected  =  [number expected in one
*     plot] x [number of plots] = D x NS x LJ x Pd x r x g(r) x P(r).
*
*     With fixed-point data, the expected total number detected
*     = [2Vt(=D2L)] x Pd x Ps x r x g(r) x P(r).
*
 8043 IF (IFX) 8044,8044,8046
 8044 IF (PD.GT.0.) GO TO 8048
      ED=D2L*E*DNRH*(SNS/2.)
      GO TO 8050
 8048 ED=D2L*E*DNRH*PD*(SNS/2.)
      GO TO 8050
 8046 ED=D2L*PD*PS*E*DNRH
*
*     With radial distance data, the total number expected in a
*     class (EXPD) is obtained by progressively adding together
*     the ED values from each arc with radii between DNRH and DNRL.
*
 8050 EXPD=ED+EXPD
*
*     In the case of perpendicular distance data, all the EXPD
*     values within a class are added together to give a TEXPD
*     value for the class (while EXPD continues to accumulate
*     from class to class).
*
      IF (IRY) 7180,7175,7180
 7175 TEXPD=EXPD+TEXPD
*
*     In preparation for the next arc in the range, DNRL is
*     now renamed DNRH (because the lowest value in one arc
*     becomes the highest in the next, going inward), while DL
*     similarly becomes DH.
*
 7180 DNRH=DNRL
      DH=DL
*
*     Similarly, in the case of radial distance data, YYL 
*     becomes YYH.
*
      IF (IMV) 8051,8051,8052
 8051 YYH=YYL
*
*     The DNR value for the next arc is DNR minus the
*     arc width (CINT).
*
 8052 DNR=DNR-CINT
*
*     Where PRR is negative, QR is kept at 1, because no animals
*     are in fact removed from the sample at such distances from
*     the observer.
*
      IF (PRR) 8055,8056,8056
 8055 QR=1.
*
*     The proportion still present in the next arc will be
*     the proportion present in the present arc less the
*     proportion expected to have been detected already.
*
*     For fixed-point data, QR remains at 1.
*
 8056 QR=(1.-PRR)*QR
*
*     If QR has fallen to a value of less than 0.001 - one 
*     animal in a thousand - the program is set to print out the
*     approximate value of rmin where this happens.
*
      IF (QR .GT. .001) GO TO 20
      IF (JRLOW) 8065,8065,20
*
*     RLOW is the inner boundary of the arc by which 99.9% of the
*     detections are expected to have been made.
*
 8065 RLOW=DNRL-CINT
      JRLOW=1
*
*     Loop 20 now ends.
*
   20 CONTINUE                                                          
*
*
*     Loop 10 is completed by determining the expected value
*     (EXPDV) in the class; this is arrived at differently for
*     radial or fixed-point and perpendicular distance data. In
*     both cases, negative EXPD values (EXPDN) are then added in.
*
      IF (IRY) 8070,8075,8070
 8070 EXPDV=EXPD
      GO TO 8100
 8075 EXPDV=TEXPD
*
*     Calculated values are printed out once the program has
*     converged on a minimum, and KPRINT has been set at 1.
*
 8100 IF (WL.GT.RMAX) GO TO 8140
*
*     If D2L has a trial value of zero, calculation of S at this stage 
*     is bypassed and a record of this non-computation retained
*     as the temporary variable MSFAIL.  This is done only when the
*     program has converged on a minimum (MTEST=1) and in the final
*     pass through Loop 10 (JJ=L10).
*
      IF (D2L) 8103,8103,8102
 8102 IF (SNS.EQ.0.) SNS=2
      S=EXPDV/(D2L*PD*(SNS/2)) + S
      GO TO 8104
 8103 IF ((MTEST.EQ.1) .AND. (JJ.EQ.L10)) MSFAIL=MSFAIL + 1
 8104 IF (KPRINT) 8140,8140,8110
*
*     For output purposes only, EXPDV is redefined as EXPDR in
*     the case of radial distance data, and as EXPDY in the case
*     of perpendicular distance data.  Negative values of EXPDV
*     are printed as '0.0' in the output because the 'observations'
*     are unreal.
*
 8110 IF (IRY) 8112,8114,8112                                          
 8112 RR=TR-(CLINT/2.)
      EXPDR=EXPDV
      IF (EXPDR) 811,812,812
  811 EXPDR=0.
  812 IF (ISHOW.LE.0) GO TO 816
      WRITE(2,8113) RR,EXPDR,OBSD
 8113 FORMAT(4H  r=,F8.1,5X,10HCalc.N(r)=,F9.2,5X,10HObsd.N(r)=,F9.1)
  816 ROUT(JJ)=RR
      CALCNR(JJ)=EXPDR 
      OBSDNR(JJ)=OBSD            
      GO TO 8116                                                        
 8114 YY=TR-(CLINT/2.)
      EXPDY=EXPDV
      IF (EXPDY) 813,814,814
  813 EXPDY=0.
  814 IF (ISHOW.LE.0) GO TO 817
      WRITE(2,8115) YY,EXPDY,OBSD       
 8115 FORMAT(4H  y=,F8.1,5X,10HCalc.N(y)=,F9.2,5X,10HObsd.N(y)=,F9.1)
  817 YOUT(JJ)=YY
      CALCNY(JJ)=EXPDY 
      OBSDNY(JJ)=OBSD   
*
*     The difference between each observed (OBSD) and expected 
*     value (EXPDV) for the class is now calculated.
*
 8116 IF (ISHOW) 8140,8140,8117
 8117 WRITE(2,8118) PRC,PRR,QR,E,TOTE,ED,EXPD,EXPDV,S
 8118 FORMAT (1X,9(F12.6,3X))
 8140 DIF=OBSD-EXPDV
*
*     Each difference between the expected (EXPDV) and observed
*     (OBSD) value is squared and a running sum of squares total
*     (TOT) computed.  The final overall sum (FUNC) - which
*     includes HTOT values where parameter values are inappropriate
*     - is then returned to the main program.
*
*     When curve-fitting is in progress, data are perpendicular
*     distance data, and KWT=1, the robust weighting procedure of
*     Wonnacott and Wonnacott is employed.  The difference between
*     observed and expected values is weighted according to the
*     magnitude of the difference between observed and calculated
*     values.  Weighting commences once initial convergence has
*     occurred (when KWT becomes 1).  A weighting term (DIF/R3S)
*     has R3S equal to 3 x Interquartile Range (Q2-Q1).  If no
*     value of R3S is supplied to the program, this term is
*     approximated by putting R3S=100.
*
*     Once the program has converged on a minimum (and LPRINT has
*     been put =1), weighting is removed to calculate the sum of
*     squares (for subsequent variance computations).  Hence the
*     procedure over the next few lines is governed by the values
*     of KWT and LPRINT.
*
      IF (IRY) 8133,8141,8133
 8141 IF (KWT) 8133,8133,8130
 8130 IF (R3S) 8134,8134,8132
 8134 R3S=100.0
      GO TO 8132
 8133 W=1.
      GO TO 8122
 8132 Z=DIF/R3S
      IF (Z.LE.1.) GO TO 8120
      W=0.
      GO TO 8122
 8120 W=(1-Z*Z)**2.
 8122 WDIFSQ=W*DIF*DIF
      WTOT=WTOT+WDIFSQ
*
      IF (LPRINT) 8146,8146,8145
 8145 DIFSQ=DIF*DIF                                                     
      TOT=TOT+DIFSQ                                                     
*
*     TR is reduced by CLINT before the subroutine goes to the
*     next class inward.
*
 8146 TR=TR-CLINT
*
*     Loop 10 ends.
*
   10 CONTINUE                                                          
*
*
      IF (KPRINT) 8148,8148,8147
 8147 WRITE(2,8150) RLOW
 8150 FORMAT(//,23H 99.9% r value (rmin) =,F7.2,3H m /)
      OPEN(UNIT=3,FILE='Outfile2.txt',STATUS='NEW',IOSTAT=IOS,ERR=286)
      IF (IRY) 8151,8155,8151
 8151 WRITE (3,8152)
 8152 FORMAT (3X,38H   Midpt.   Calculated     Observed   )
      DO 8153 JJ=1,L10
      JV=L10-JJ+1
      LIMIT=DMAX/CLINT
      IF (JJ.GT.LIMIT) GO TO 8148
      WRITE (3,8154) ROUT(JV),CALCNR(JV),OBSDNR(JV)
 8154 FORMAT (5X,F6.1,6X,F7.2,7X,F6.1,4X)
 8153 CONTINUE
      GO TO 8148
 8155 WRITE (3,8156)
 8156 FORMAT (3X,38H   Midpt.   Calculated     Observed   )
      DO 8157 JJ=1,L10
      JV=L10-JJ+1
      LIMIT=DMAX/CLINT
      IF (JJ.GT.LIMIT) GO TO 8148
      WRITE (3,8158) YOUT(JV),CALCNY(JV),OBSDNY(JV)
 8158 FORMAT (5X,F6.1,6X,F7.2,7X,F6.1,4X)
 8157 CONTINUE                                 
      CLOSE(UNIT=3)
*
 8148 IF (LPRINT) 8160,8160,8161
 8160 FUNC=WTOT+HTOT
      GO TO 8180
 8161 FUNC=TOT+HTOT
      GO TO 8180
  286 WRITE(6,287)OUTFILE,IOS
  287 FORMAT(' Error opening ',A40,' - IOS = ',I6)
*
*     If F(2) has been given a high negative value, FUNC is set to
*     a high value before returning to the main program.
*
 8180 IF (F(2).GE.QMIN)GO TO 8190                 
       FUNC=(ABS(F(2)-QMIN-1.))*FUNC                                     
 8190 IF (KPRINT) 8200,8200,8191
 8191 WRITE(2,8192) FUNC
 8192 FORMAT(X,28HOLS Difference at Minimum = ,F15.6/)
 8200 RETURN                                                            
      END SUBROUTINE GIVEF
*
*
       SUBROUTINE SRANDNAG(ISEED)
*
*  This subroutine sets the integer seed to be used with the
*  companion RANDNAG function to the value of ISEED.  A flag is
*  set to indicate that the sequence of pseudo-random numbers
*  for the specified seed should start from the beginning.
*
      DOUBLE PRECISION JSEED,IFRST
      COMMON / SEED /JSEED,IFRST
*
      JSEED=ISEED
      IFRST=0
*
      END SUBROUTINE SRANDNAG

      REAL FUNCTION RANDNAG()
*
*
*  This function returns a pseudo-random number for each invocation.
*  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
*  standard number generator whose Pascal code appears in the
*  article:
*     Park, Steven K. and Miller, Keith W., "Random Number
*     Generators: Good Ones are Hard to Find", Communications of
*     the ACM, October, 1988.
*  It is used in the form recommended by the Numerical
*  Algorithms Group (NAG) for UNIX systems.
*
*
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)
*
      COMMON / SEED /JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN
*
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST=1
      ENDIF
*
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
	  TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RANDNAG = REAL(NEXTN)/REAL(MODLUS)
*
      RETURN
      END FUNCTION RANDNAG

      BLOCKDATA RANDBD
      COMMON / SEED /JSEED,IFRST
*
      DATA JSEED,IFRST/123456789,0/
*
      END BLOCKDATA RANDBD
