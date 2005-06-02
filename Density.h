/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@class DensityController;

@interface Density : NSDocument
{
    IBOutlet DensityController *controller;
	NSString *header;
    NSString *outFile;
    NSString *graphFile;
    int iprint, jprint, ishow, maxjb;
    double elevations[MAX_OBSERVATIONS];
    calc_params params;
    NSTimer *calculationTimer;
    BOOL complete;
    NSString *completeMsg;
}

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

- (void)calculate;

- (void)calculationThread:(id)ignored;

- (void)calculationWork;
@end
