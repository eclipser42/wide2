/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@class DensityController;

@interface Observation : NSObject
{
    double distance;
    int groupSize;
    double angle;
    double elevation;
}
- (double)distance;
- (void)setDistance:(double)newDistance;
- (int)groupSize;
- (void)setGroupSize:(int)newGroupSize;
- (double)angle;
- (void)setAngle:(double)newAngle;
- (double)elevation;
- (void)setElevation:(double)newElevation;
@end

@interface Density : NSDocument
#if __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 1060
                                < NSWindowDelegate,NSTableViewDelegate >
#endif
{
    IBOutlet DensityController *controller;
	NSString *header;
	int currentResultsIndex;
    BOOL elevationsAreSupplied;
    double elevations[MAX_OBSERVATIONS];
	NSMutableArray *observations;
    IBOutlet NSArrayController *observationsController;
    calc_params params;
    BOOL complete;
    NSString *completeMsg;
}

- (void)setHeader:(NSString*)newHeader;
- (int)ifx;
- (void)setIfx:(int)ifx;
- (int)iry;
- (void)setIry:(int)iry;
- (int)kdt;
- (void)setKdt:(int)kdt;
- (double)dist;
- (void)setDist:(double)dist;
- (int)km;
- (void)setKm:(int)km;
- (double)ltmin;
- (void)setLtmin:(double)ltmin;
- (int)ns;
- (void)setNs:(int)ns;
- (double)pd;
- (void)setPd:(double)pd;
- (double)vgh;
- (void)setVgh:(double)vgh;
- (double)durn;
- (void)setDurn:(double)durn;
- (double)movementRate;
- (void)setMovementRate:(double)rate;
- (double)ps;
- (void)setPs:(double)ps;
- (void)setElevationsAreSupplied:(BOOL)newElevationsAreSupplied;
- (double)thh;
- (void)setThh:(double)thh;

- (double)stt;
- (void)setStt:(double)stt;

- (void)insertObject:(Observation *)o inObservationsAtIndex:(int)index;
- (void)removeObjectFromObservationsAtIndex:(int)index;
- (void)startObservingObservation:(Observation *)o;
- (void)stopObservingObservation:(Observation *)o;
- (int)maxjb;
- (void)setMaxjb:(int)maxjb;
- (double)clint;
- (void)setClint:(double)clint;
- (double)f0;
- (void)setF0:(double)f0;
- (double)f1;
- (void)setF1:(double)f1;
- (double)f2;
- (void)setF2:(double)f2;
- (double)f3;
- (void)setF3:(double)f3;
- (double)step0;
- (void)setStep0:(double)step0;
- (double)step1;
- (void)setStep1:(double)step1;
- (double)step2;
- (void)setStep2:(double)step2;
- (int)iprint;
- (void)setIprint:(int)iprint;
- (int)ishow;
- (void)setIshow:(int)ishow;
- (int)jprint;
- (void)setJprint:(int)jprint;
- (int)nvals;
- (int)currentIteration;
- (int)maxIteration;

/**
 * Parse input with the observation data in rows, separated by either commas or spaces
 **/
- (BOOL)parseInputRows:(NSString *)input withCommas:(BOOL)requireCommas;

/**
 * Parse input with the observation data in columns
 **/
- (BOOL)parseInputColumns:(NSString *)input;
- (BOOL)parseInputOldColumns:(NSString *)input;

- (void)endEditingSavingCurrentResponder;
- (NSString *)resultsFileName;
- (NSString *)graphDataFileName;
- (IBAction)viewCurrentResults:(id)sender;
- (IBAction)graphCurrentResults:(id)sender;

- (void)launchCalculation;

- (void)calculationThread:(id)ignored;

- (void)calculationWork;
@end
