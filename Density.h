/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@interface Density : NSDocument
{
	NSString *header;
    NSString *outFile;
    NSString *graphFile;
    int iprint,jprint,ishow,maxjb;
    calc_params params;
    BOOL complete;
    NSString *completeMsg;
}

- (int)nvals;

- (void)calculate;

/**
 * Parse input with the observation data in rows, separated by either commas or spaces
 **/
- (BOOL)parseInputRows:(NSString *)input withCommas:(BOOL)requireCommas;

/**
 * Parse input with the observation data in columns
 **/
- (BOOL)parseInputColumns:(NSString *)input;
@end
