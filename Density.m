//
//  Density.m
//  CocoaTest
//
//  Created by David Morgan on Mon Jul 26 2004.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

#import "Density.h"

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
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"Density";
}

- (NSData *)dataRepresentationOfType:(NSString *)aType
{
    NSAssert([aType isEqualToString:@"DensityCensus"], @"Unknown type");
    NSMutableString *contents = [NSMutableString stringWithCapacity:2048];
    [contents appendFormat:@"'%@',%c", header, 10];
    [contents appendFormat:@"%d, %f, %f, %d, %d, %f, ", params.nvals, params.clint, params.stt, params.numa, params.numo, params.dist];
    [contents appendFormat:@"%f, %d, %d, %d, %d, %d,%c", params.thh, params.ltmin, params.ltmax, params.ifx, params.iry, params.ns, 10];
    [contents appendFormat:@"%d, %d, %d, %d, %d, %d, ", params.km, params.imv, params.kdt, iprint, jprint, ishow];
    [contents appendFormat:@"%d, %f, %f, %d, %d, %f, %f,%c", maxjb, params.r3s, params.vgh, params.it, params.iv, params.pd, params.ps, 10];
    for (int i = 0; i < NUM_SHAPE_PARAMS; i++) {
        [contents appendFormat:@"%f, ", params.f[i]];
    }
    for (int i = 0; i < NUM_SHAPE_PARAMS; i++) {
        [contents appendFormat:@"%f, ", params.step[i]];
    }
    [contents appendFormat:@"%c", 10];
    for (int i = 0; i < params.nvals; i += 10) {
        for (int j = i; j < i + 10 && j < params.nvals; j++) {
            [contents appendFormat:@"%f, ", params.r[j]];
        }
        [contents appendFormat:@"%c", 10];
    }
	for (int i = 0; i < params.nvals; i += 10) {
        for (int j = i; j < i + 10 && j < params.nvals; j++) {
            [contents appendFormat:@"%d, ", params.nsize[j]];
        }
        [contents appendFormat:@"%c", 10];
	}
	if (params.iry == 1) {
		for (int i = 0; i < params.nvals; i += 10) {
			for (int j = i; j < i + 10 && j < params.nvals; j++) {
				[contents appendFormat:@"%d, ", params.angle[j]];
			}
			[contents appendFormat:@"%c", 10];
		}
	}
    return [contents dataUsingEncoding:NSUTF8StringEncoding];
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
    // Insert code here to read your document from the given data.  You can also choose to override -loadFileWrapperRepresentation:ofType: or -readFromFile:ofType: instead.
	NSString *contents = [NSString allocWithZone:[self zone]];
    contents = [contents initWithData:data encoding:NSUTF8StringEncoding];

    if ([self parseInputCommas:contents]) {
        NSLog(@"Parsed comma-separated data");
    } else if ([self parseInputSpaces:contents]) {
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

- (void)calculate
{
    NSLog(@"Calculating to %@ and %@", outFile, graphFile);
    const char * hdr = [header UTF8String];
    const char * out = [outFile UTF8String];
    const char * graph = [graphFile UTF8String];
    unlink([[NSFileManager defaultManager] fileSystemRepresentationWithPath:outFile]);
    unlink([[NSFileManager defaultManager] fileSystemRepresentationWithPath:graphFile]);

    params.iprint = iprint;
    params.jprint = jprint;
    params.ishow = ishow;
    params.maxjb = maxjb;
    params.complete = 0;
    calculate_density(&params, hdr, out, graph, strlen(hdr), strlen(out), strlen(graph));
    NSString *message;
    if (params.complete) {
        message = [NSString stringWithFormat:@"Calculation complete %c%cDensity estimate: %3g  Standard error: %3g%c%cDetailed results in %@ and %CgraphData",
            10, 10, params.estden, params.sden, 10, 10, outFile, 0x2026];
    } else {
        message = [NSString stringWithFormat:@"Calculation failed: Detailed results in %@", outFile];
    }
    [self setValue:message forKey:@"completeMsg"];
}

- (BOOL)parseInputSpaces:(NSString *)input
{
	NSScanner* scanner = [NSScanner scannerWithString:input];

	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the opening quote */
	if (![scanner scanUpToString:@"'" intoString:&header]) return NO;
    [header retain];
	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the closing quote */

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
	if (![scanner scanInt:&params.ltmin]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ltmax]) return NO;
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
	if (![scanner scanInt:&params.it]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iv]) return NO;
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
		if (![scanner scanDouble:(params.r + i)]) return NO;
	}
	for (int i = 0; i < params.nvals; i++) {
		if (![scanner scanInt:(params.nsize + i)]) return NO;
	}
	if (params.iry == 1) {
		for (int i = 0; i < params.nvals; i++) {
			if (![scanner scanDouble:(params.angle + i)]) 
			return NO;
		}
	}
	
	return YES;
}

- (BOOL)parseInputCommas:(NSString *)input
{
	NSScanner* scanner = [NSScanner scannerWithString:input];

	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the opening quote */
	if (![scanner scanUpToString:@"'" intoString:&header]) return NO;
    [header retain];
	if (![scanner scanString:@"'" intoString:nil]) return NO; /* Skip the closing quote */
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */

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
	if (![scanner scanInt:&params.ltmin]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.ltmax]) return NO;
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
	if (![scanner scanInt:&params.it]) return NO;
	if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	if (![scanner scanInt:&params.iv]) return NO;
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
		if (![scanner scanDouble:(params.r + i)]) return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
	for (int i = 0; i < params.nvals; i++) {
		if (![scanner scanInt:(params.nsize + i)]) {
            return NO;
        }
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
	}
	if (params.iry == 1) {
	for (int i = 0; i < params.nvals; i++) {
		if (![scanner scanDouble:(params.angle + i)]) 
			return NO;
		if (![scanner scanString:@"," intoString:nil]) return NO; /* Skip the comma separator */
		}
	}

	return YES;
}

@end
