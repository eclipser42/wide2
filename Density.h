/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@interface Density : NSDocument
{
	NSString *header;
    NSString *outFile;
    int maxjb;
    calc_params params;
}

- (int)nvals;

- (void)calculate;

/**
 * Parse input with the observation data separated by spaces
 **/
- (BOOL)parseInputSpaces:(NSString *)input;

/**
 * Parse input with the observation data separated by commas
 **/
- (BOOL)parseInputCommas:(NSString *)input;
@end
