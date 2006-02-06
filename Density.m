//
//  Density.m
//  CocoaTest
//
//  Created by David Morgan on Mon Jul 26 2004.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

#import "Density.h"
#import "DensityController.h"

#import "Foundation/NSFileManager.h"
#import <unistd.h>

@implementation Density

- (id)init
{
    self = [super init];
    if (self) {

        // Add your subclass-specific initialization here.
		header = [[NSString allocWithZone:[self zone]] init];
		outFile = [NSString allocWithZone:[self zone]];
		graphFile = [NSString allocWithZone:[self zone]];
        completeMsg = [[NSString allocWithZone:[self zone]] init];

        // Set defaults
        params.stt = 0;
        params.clint = 20;
        maxjb = 500;
        params.f[0] = 10;
        params.f[1] = 0.1;
        params.f[2] = 10;
        params.f[3] = 10;
        params.step[0] = 0;
        params.step[1] = 0;
        params.step[2] = 0;
        params.step[3] = 0;
        params.iprint = 0;
        params.jprint = 0;
        params.ishow = 0;
        params.imv = 0;
        params.r3s = 100;

        // If an error occurs here, send a [self release] message and return nil.

    }
    return self;
}

- (void)dealloc
{
    [header release];
    [outFile release];
    [graphFile release];
    [completeMsg release];
    [super dealloc];
}

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers,
    // you should remove this method and override -makeWindowControllers instead.
    return @"Density";
}

- (NSData *)dataRepresentationOfType:(NSString *)aType
{
    NSAssert([aType isEqualToString:@"DensityCensus"], @"Unknown type");
    NSMutableString *contents = [NSMutableString stringWithCapacity:2048];
    [contents appendFormat:@"'%@'%c", header, 10];
    [contents appendFormat:@"%d, %d, %d, %g, %d, %g, %d, %g, %g%c", params.ifx, params.iry, params.kdt, params.dist, params.km, params.ltmin, params.ns, params.pd, params.vgh, 10];
    [contents appendFormat:@"%g, %g, %g, %g%c%c", params.durn, params.rate, params.ps, params.thh, 10, 10];
    for (int i = 0; i < params.nvals; i++) {
        if (params.iry == 1) {
            [contents appendFormat:@"%g, %d, %g", params.r[i], params.nsize[i], params.angle[i]];
        } else {
            [contents appendFormat:@"%g, %d", params.r[i], params.nsize[i]];
        }
        if (elevations[i] == -1) {
            [contents appendFormat:@"%c", 10];
        } else {
            [contents appendFormat:@", %g%c", elevations[i], 10];
        }
    }
    [contents appendFormat:@"%c%g, %d, %g, %g, %g, %g, %g%c", 10, params.stt, maxjb, params.clint, params.f[0], params.f[1], params.f[2], params.f[3], 10];
    [contents appendFormat:@"%g, %g, %g%c", params.step[0], params.step[1], params.step[2], 10];
    [contents appendFormat:@"%d, %d, %d%c", iprint, jprint, ishow, 10];
    [contents appendFormat:@"%d, %g, %g%c", params.imv, params.r3s, params.step[3], 10];
    return [contents dataUsingEncoding:NSUTF8StringEncoding];
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
    // Insert code here to read your document from the given data.
    // You can also choose to override -loadFileWrapperRepresentation:ofType: or -readFromFile:ofType: instead.
	NSString *contents = [NSString allocWithZone:[self zone]];
    contents = [contents initWithData:data encoding:NSUTF8StringEncoding];

    if ([self parseInputColumns:contents]) {
        NSLog(@"Parsed column-arranged data");
    } else if ([self parseInputOldColumns:contents]) {
        NSLog(@"Parsed old column-arranged data");
    } else if ([self parseInputRows:contents withCommas:YES]) {
        NSLog(@"Parsed comma-separated data");
    } else if ([self parseInputRows:contents withCommas:NO]) {
        NSLog(@"Parsed space-separated data");
    } else {
        NSLog(@"Parsing failed");
        return NO;
    }
    iprint = params.iprint;
    jprint = params.jprint;
    ishow = params.ishow;
    maxjb = params.maxjb;

    outFile = [[[self fileName] stringByDeletingPathExtension] stringByAppendingPathExtension:@"results"];
    if ([outFile isEqualToString:[self fileName]]) {
        outFile = [outFile stringByAppendingPathExtension:@"results"];
    }
    [outFile retain];

    graphFile = [[[self fileName] stringByDeletingPathExtension] stringByAppendingPathExtension:@"graphData"];
    if ([graphFile isEqualToString:[self fileName]]) {
        graphFile = [graphFile stringByAppendingPathExtension:@"graphData"];
    }
    [graphFile retain];

    return YES;
}

- (int)nvals
{
    return params.nvals;
}

- (int)currentIteration
{
    return params.bootstrap;
}

- (int)maxIteration
{
    return maxjb;
}

- (BOOL)parseInputColumns:(NSString *)input
{
	NSScanner* scanner = [NSScanner scannerWithString:input];

	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the opening quote */
	if (![scanner scanUpToString:@"'" intoString:&header]) return NO;
    [header retain];
	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the closing quote */

	if (![scanner scanInt:&params.ifx]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iry]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.kdt]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.dist]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
    if (![scanner scanInt:&params.km]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ltmin]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ns]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.pd]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.vgh]) return NO;

	if (![scanner scanDouble:&params.durn]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.rate]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ps]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.thh]) return NO;

    int i = 0;
    unsigned bottomParams;
    while (i < MAX_OBSERVATIONS) {
        unsigned potentialLine = [scanner scanLocation];
        if ([scanner scanDouble:(params.r + i)]) {
            bottomParams = potentialLine; /* potentialLine really was a new line */
        } else {
            i--; /* potentialLine was not a new line of observations, but the middle of the bottom parameters */
            break; /* finish and back up */
        }
        if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
        if (![scanner scanInt:(params.nsize + i)]) break;
        if (params.iry == 1) {
            if (![scanner scanString:@"," intoString:nil]) break; /* Skip the comma separator */
            if (![scanner scanDouble:(params.angle + i)]) return NO;
        }
        // and maybe there's an elevation parameter here too
        if ([scanner scanString:@"," intoString:nil]) {
            if (![scanner scanDouble:(elevations + i)]) return NO;
        } else {
            elevations[i] = -1;
        }
        i++;
    }
    if (i >= 1) {
        params.nvals = i;
        [scanner setScanLocation:bottomParams];
    } else {
        return NO;
    }

	if (![scanner scanDouble:&params.stt]) return NO;
    if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.maxjb]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.clint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
    for (int i = 0; i < NUM_SHAPE_PARAMS - 1; i++) {
		if (![scanner scanDouble:&params.f[i]]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
    if (![scanner scanDouble:&params.f[NUM_SHAPE_PARAMS - 1]]) return NO;

    if (![scanner scanDouble:&params.step[0]]) return NO;
    if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
    if (![scanner scanDouble:&params.step[1]]) return NO;
    if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
    if (![scanner scanDouble:&params.step[2]]) return NO;

	if (![scanner scanInt:&params.iprint]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
	if (![scanner scanInt:&params.jprint]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
	if (![scanner scanInt:&params.ishow]) return YES;

    if (![scanner scanInt:&params.imv]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
	if (![scanner scanDouble:&params.r3s]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
    if (![scanner scanDouble:&params.step[NUM_SHAPE_PARAMS - 1]]) return YES;

    return YES;
}

- (BOOL)parseInputOldColumns:(NSString *)input
{
	NSScanner* scanner = [NSScanner scannerWithString:input];

	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the opening quote */
	if (![scanner scanUpToString:@"'" intoString:&header]) return NO;
    [header retain];
	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the closing quote */

	if (![scanner scanDouble:&params.clint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.stt]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.numa]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.numo]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.dist]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.thh]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ltmin]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ltmax]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ifx]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iry]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ns]) return NO;

    if (![scanner scanInt:&params.km]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.imv]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.kdt]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iprint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.jprint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ishow]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.maxjb]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.r3s]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.vgh]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.durn]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.rate]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.pd]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ps]) return NO;

    for (int i = 0; i < NUM_SHAPE_PARAMS - 1; i++) {
		if (![scanner scanDouble:&params.f[i]]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
    if (![scanner scanDouble:&params.f[NUM_SHAPE_PARAMS - 1]]) return NO;

	for (int i = 0; i < NUM_SHAPE_PARAMS - 1; i++) {
		if (![scanner scanDouble:&params.step[i]]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
    if (![scanner scanDouble:&params.step[NUM_SHAPE_PARAMS - 1]]) return NO;

    int i = 0;
    while (YES) {
		if (![scanner scanDouble:(params.r + i)]) break;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
		if (![scanner scanInt:(params.nsize + i)]) return NO;
        if (params.iry == 1) {
            if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
			if (![scanner scanDouble:(params.angle + i)]) return NO;
		}
        elevations[i] = -1;
        i++;
	}
    params.nvals = i;

    return (i >= 1);
}

- (BOOL)parseInputRows:(NSString *)input withCommas:(BOOL)requireCommas
{
	NSScanner* scanner = [NSScanner scannerWithString:input];

	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the opening quote */
	if (![scanner scanUpToString:@"'" intoString:&header]) return NO;
    [header retain];
	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the closing quote */
    if (requireCommas) {
        if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
    }

	if (![scanner scanInt:&params.nvals]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.clint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.stt]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.numa]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.numo]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.dist]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.thh]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ltmin]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ltmax]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ifx]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iry]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ns]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.km]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.imv]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.kdt]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iprint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.jprint]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ishow]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.maxjb]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.r3s]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.vgh]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.durn]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.rate]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.pd]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanDouble:&params.ps]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	for (int i = 0; i < NUM_SHAPE_PARAMS; i++) {
		if (![scanner scanDouble:&params.f[i]]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
	for (int i = 0; i < NUM_SHAPE_PARAMS; i++) {
		if (![scanner scanDouble:&params.step[i]]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}

	for (int i = 0; i < params.nvals; i++) {
        elevations[i] = -1;
		if (![scanner scanDouble:(params.r + i)]) return NO;
        if (requireCommas) {
            if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
        }
	}
	for (int i = 0; i < params.nvals; i++) {
		if (![scanner scanInt:(params.nsize + i)]) return NO;
        if (requireCommas) {
            if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
        }
    }
    if (params.iry == 1) {
        for (int i = 0; i < params.nvals; i++) {
            if (![scanner scanDouble:(params.angle + i)]) return NO;
            if (requireCommas) {
                if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
            }
        }
    }

	return YES;
}

- (void)calculate
{
    [self setValue:@"" forKey:@"completeMsg"];
    //[self calculationWork];

    [NSThread detachNewThreadSelector:@selector(calculationThread:) toTarget:self withObject:nil];
}

- (void)calculationThread:(id)ignored
{
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    [self calculationWork];
    [pool release];
}

const double PI = 3.14159254;
double deg2rad(double deg) {
    return deg / 180 * PI;
}

- (void)calculationWork
{
    NSLog(@"Calculating to %@ and %@", outFile, graphFile);
    params.iprint = iprint;
    params.jprint = jprint;
    params.ishow = ishow;
    params.maxjb = maxjb;
    params.complete = 0;

    const char * hdr = [header UTF8String];
    const char * out = [outFile UTF8String];
    const char * graph = [graphFile UTF8String];
    unlink([[NSFileManager defaultManager] fileSystemRepresentationWithPath:outFile]);
    unlink([[NSFileManager defaultManager] fileSystemRepresentationWithPath:graphFile]);

    // Compute THH and LTMAX
    int elevationCount = 0;
    double sumOfSquaredElevations = 0;
    params.ltmax = 0;
    for (int i = 0; i < params.nvals; i++) {
        if (elevations[i] != -1) {
            double horizontalDistance = params.r[i];
            if (params.km == 2) horizontalDistance *= 1000;
            elevationCount++;
            double elevation = horizontalDistance * tan(deg2rad(elevations[i]));
            sumOfSquaredElevations += elevation * elevation;
        }
        if (params.r[i] > params.ltmax) {
            params.ltmax = params.r[i];
        }
    }
    if (elevationCount) {
        double thh = sqrt(sumOfSquaredElevations / elevationCount);
        NSLog(@"Overriding THH %g with %g", params.thh, thh);
        params.thh = thh;
    }

    // Compute NUMA and NUMO
    params.numa = params.numo = 0;
    if (params.ifx) {
        for (int i = 0; i < params.nvals; i++) {
            params.numa += params.nsize[i];
        }
    } else {
        for (int i = 0; i < params.nvals; i++) {
            if (params.r[i] > 0) {
                params.numa += params.nsize[i];
            } else {
                params.numo += params.nsize[i];
            }
        }
    }

    calculate_density_(&params, hdr, out, graph, strlen(hdr), strlen(out), strlen(graph));

    NSString *message;
    if (params.complete) {
        message = [NSString stringWithFormat:@"Calculation complete %c%cDensity estimate: %.3g  Standard error: %.2g%c%cDetailed results in %@ and %@",
            10, 10, params.estden, params.sden, 10, 10, outFile, [graphFile lastPathComponent]];
    } else {
        message = [NSString stringWithFormat:@"Calculation failed: Detailed results in %@", outFile];
    }
    [self setValue:message forKey:@"completeMsg"];
    [controller finishCalculation];
}

@end
