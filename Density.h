/* Density */

#import "mnps2.h"
#import <Cocoa/Cocoa.h>

@interface Density : NSDocument
{
	NSString *header;
    NSString *outFile;
    calc_params params;
}

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
