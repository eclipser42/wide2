/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@class DensityController;

@interface Density : NSDocument
{
    //IBOutlet DensityController *controller;
	NSString *header;
    NSString *outFile;
    NSString *graphFile;
    double elevations[MAX_OBSERVATIONS];
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
- (double)thh;
- (void)setThh:(double)thh;

- (double)stt;
- (void)setStt:(double)stt;
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
- (double)step3;
- (void)setStep3:(double)step3;
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

- (void)launchCalculation;

- (void)calculationThread:(id)ignored;

- (void)calculationWork;
@end
