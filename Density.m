//
//  Density.m
//  CocoaTest
//
//  Created by James Clough on 9 Jan 2007.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

#import "Density.h"
#import "DensityController.h"
#import "GraphData.h"

#import "Foundation/NSFileManager.h"
#import <unistd.h>

@implementation Observation
- (double)distance                        { return distance; }
- (void)setDistance:(double)newDistance   { distance = newDistance; }
- (int)groupSize                          { return groupSize; }
- (void)setGroupSize:(int)newGroupSize;   { groupSize = newGroupSize; }
- (double)angle                           { return angle; }
- (void)setAngle:(double)newAngle         { angle = newAngle; }
- (double)elevation                       { return elevation; }
- (void)setElevation:(double)newElevation { elevation = newElevation; }
@end

@implementation Density

- (id)init
{
    self = [super init];
    if (self) {

        // Add your subclass-specific initialization here.
		header = @"[Data Set Description]";
        currentResultsIndex = 0;
		observations = [[NSMutableArray alloc] init];
        completeMsg = @"";

        // Set defaults
        params.ns = 2;
        params.thh = 0;
        params.stt = 0;
        params.clint = 20;
        params.maxjb = 500;
        params.pd = 1;
        params.f[0] = 15;
        params.f[1] = 0.1;
        params.f[2] = 5;
        params.f[3] = 250;
        params.step[0] = 6;
        params.step[1] = 0.1;
        params.step[2] = 5;
        params.step[3] = 0;
        params.iprint = 0;
        params.jprint = 0;
        params.ishow = 0;

        // If an error occurs here, send a [self release] message and return nil.

    }
    return self;
}

- (void)dealloc
{
    [header release];
	[observations release];
    [completeMsg release];

    [super dealloc];
}

#pragma mark NSDocument overrides

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers,
    // you should remove this method and override -makeWindowControllers instead.
    return @"Density";
}

- (NSData *)dataRepresentationOfType:(NSString *)aType
{
    NSAssert([aType isEqualToString:@"WildlifeDensity Dataset"], @"Save requested to file of unknown type");
    NSMutableString *contents = [NSMutableString stringWithCapacity:2048];
    [contents appendFormat:@"'%@'%c", header, 10];
    [contents appendFormat:@"%d, %d, %d, %g, %d, %g, %d, %g, %g%c", params.ifx, params.iry, params.kdt, params.dist, params.km, params.ltmin, params.ns, params.pd, params.vgh, 10];
    [contents appendFormat:@"%g, %g, %g, %g%c%c", params.durn, params.rate, params.ps, params.thh, 10, 10];
    for (int i = 0; i < [observations count]; ++i) {
        Observation *observation = [observations objectAtIndex:i];
        if (params.iry == 1) {
            [contents appendFormat:@"%g, %d, %g",
                [observation distance],
                [observation groupSize],
                [observation angle]];
        } else {
            [contents appendFormat:@"%g, %d",
                [observation distance],
                [observation groupSize]];
        }
        double elevationAngle = [observation elevation];
        if (elevationAngle == -1) {
            [contents appendFormat:@"%c", 10];
        } else {
            [contents appendFormat:@", %g%c", elevationAngle, 10];
        }
    }
    if ([observations count] == 0) {
        // write at least one observation line so the file format is valid
        if (params.iry == 1) {
            [contents appendFormat:@"0, 0, 0%c", 10];
        } else {
            [contents appendFormat:@"0, 0%c", 10];
        }
    }
    [contents appendFormat:@"%c%g, %d, %g, %g, %g, %g, %g%c", 10, params.stt, params.maxjb, params.clint, params.f[0], params.f[1], params.f[2], params.f[3], 10];
    [contents appendFormat:@"%g, %g, %g%c", params.step[0], params.step[1], params.step[2], 10];
    [contents appendFormat:@"%d, %d, %d%c", params.iprint, params.jprint, params.ishow, 10];
    [contents appendFormat:@"%d, %g, %g%c", 0 /* imv */, 0.0 /* r3s */, params.step[3], 10];
    return [contents dataUsingEncoding:NSUTF8StringEncoding];
}

- (BOOL)canConcurrentlyReadDocumentsOfType:(NSString *)aType
{
    // Our document-reading code can be safely executed concurrently, in non-main threads.
    return YES;
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
    // Insert code here to read your document from the given data.
    // You can also choose to override -loadFileWrapperRepresentation:ofType: or -readFromFile:ofType: instead.
	NSString *contents = [NSString allocWithZone:[self zone]];
    contents = [contents initWithData:data encoding:NSUTF8StringEncoding];

    elevationsAreSupplied = NO;
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

    [observations removeAllObjects];
    for (int i = 0; i < params.nvals; ++i) {
        Observation *observation = [[Observation alloc] init];
        [observation setDistance:params.r[i]];
        [observation setGroupSize:params.nsize[i]];
        [observation setAngle:params.angle[i]];
        [observation setElevation:elevations[i]];
        [self startObservingObservation:observation];
        [observations addObject:observation];
    }

    NSFileManager *fileMgr = [NSFileManager defaultManager];
    do {
        ++currentResultsIndex;
    } while ([fileMgr fileExistsAtPath:[self resultsFileName]]
             || [fileMgr fileExistsAtPath:[self graphDataFileName]]);
    --currentResultsIndex;

    return YES;
}

- (void)windowControllerDidLoadNib:(NSWindowController *)windowController
{
    [controller updateTableColumns];
    [[windowController window] setDelegate:self];
}

#define NL_STRING ([NSString stringWithFormat:@"%c", 10])
#define CR_STRING ([NSString stringWithFormat:@"%c", 13])
- (void)saveDocumentWithDelegate:(id)delegate didSaveSelector:(SEL)didSaveSelector contextInfo:(void *)contextInfo
{
    [self endEditingSavingCurrentResponder];

    // Trim the header if necessary
    //NSArray* headerLines = [header componentsSeparatedByCharactersInSet:[NSCharacterSet newlineCharacterSet]];
    NSArray* tigerLines = [header componentsSeparatedByString:NL_STRING];
    NSArray* headerLines = [[tigerLines objectAtIndex:0] componentsSeparatedByString:CR_STRING];
    if ([headerLines count] > 1 || [tigerLines count] > 1) {
        NSLog(@"Trimming header to just first line (discarding %d later lines)", [headerLines count] - 1);
        [self setValue:[headerLines objectAtIndex:0] forKey:@"header"];
    }
    
    [super saveDocumentWithDelegate:delegate didSaveSelector:didSaveSelector contextInfo:contextInfo];
}

- (BOOL)prepareSavePanel:(NSSavePanel *)savePanel
{
    [savePanel setAllowsOtherFileTypes:YES];
    [savePanel setCanSelectHiddenExtension:YES];
    [savePanel setRequiredFileType:@"WDdata"];
    return YES;
}

- (void) document:(NSDocument *)doc didSave:(BOOL)didSave contextInfo:(void *)contextInfo
{
    NSLog(@"Saved %d", didSave);
}

#pragma mark NSWindow delegate

- (NSUndoManager *)windowWillReturnUndoManager:(NSWindow *)window
{
    return [self undoManager];
}

#pragma mark NSTableView delegate

- (void)paste:(id)sender
{
    NSPasteboard *pb = [NSPasteboard generalPasteboard];
    NSArray *pasteTypes = [NSArray arrayWithObjects:
        NSTabularTextPboardType, NSStringPboardType, nil];
    NSString *bestType = [pb availableTypeFromArray:pasteTypes];
    if (bestType != nil) {
        NSString *contents = [pb stringForType:bestType];
        [self willChangeValueForKey:@"nvals"];
        
        // Parse pasteboard into cells
        NSScanner *scanner = [NSScanner scannerWithString:contents];
        while (true) {
            double distance, groupSize, angle = 0, elevation = 0;
            bool rowIsValid = YES;
            Observation *newRow = [[Observation alloc] init];
            
            if ([scanner scanDouble:&distance])
                [newRow setDistance:distance];
            else
                rowIsValid = NO;
            if (rowIsValid && [scanner scanDouble:&groupSize])
                [newRow setGroupSize:groupSize];
            else
                rowIsValid = NO;
            if (rowIsValid && params.iry == 1) {
                if ([scanner scanDouble:&angle])
                    [newRow setAngle:angle];
                else
                    rowIsValid = NO;
            }
            if (rowIsValid && elevationsAreSupplied) {
                if ([scanner scanDouble:&elevation])
                    [newRow setElevation:elevation];
                else
                    rowIsValid = NO;
            }
            if (rowIsValid)
            {
                // TODO: insert all rows or put up an alert
                NSIndexSet *selection = [observationsController selectionIndexes];
                unsigned insertionPoint;
                if ([selection count] == 0) {
                    insertionPoint = [observations count];
                    NSLog(@"inserting at end (%d)", insertionPoint);
                } else {
                    insertionPoint = [selection lastIndex] + 1;
                    NSLog(@"inserting at %d", insertionPoint);
                }
                [self insertObject:newRow inObservationsAtIndex:insertionPoint];
                [[self undoManager] setActionName:@"Paste"];
                if ([scanner isAtEnd])
                    break;
            } else {
                NSLog(@"Failed parsing %@ after %d", contents, [scanner scanLocation]);
                break;
            }
        }
        [self didChangeValueForKey:@"nvals"];
    }
}

#pragma mark Key-Value Observing

// Wrap [NSObject setValue:forKey:] in a name that can be captured by [NSUndoManager forwardInvocation:]
- (void)restoreValue:(id)value forKey:(NSString *)key
{
    [self setValue:value forKey:key];
}

- (void)setHeader:(NSString*)newHeader
{
    [[self undoManager] registerUndoWithTarget:self selector:@selector(setHeader:) object:header];
	[header release];
	header = newHeader;
	[header retain];
}

#pragma mark - Method Tab

- (int)ifx
{
    return params.ifx;
}

- (void)setIfx:(int)ifx
{
    // Note this is a model change that doesn't register to undo itself: all callers must regsiter
    // undo & redo of calls to setIfx:
    params.ifx = ifx;
}

- (int)iry
{
    return params.iry;
}

- (void)setIry:(int)iry;
{
    // Note this is a model change that doesn't register to undo itself: all callers must regsiter
    // undo & redo of calls to setIry:
    params.iry = iry;
    [controller updateTableColumns];
}

- (int)kdt
{
    return params.kdt;
}

- (void)setKdt:(int)kdt;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setKdt:params.kdt];
    params.kdt = kdt;
}

- (double)dist
{
    return params.dist;
}

- (void)setDist:(double)dist;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setDist:params.dist];
    params.dist = dist;
}

- (int)km
{
    return params.km;
}

- (void)setKm:(int)km;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setKm:params.km];
    params.km = km;
}

- (double)ltmin
{
    return params.ltmin;
}

- (void)setLtmin:(double)ltmin;
{
    params.ltmin = ltmin;
}

- (int)ns
{
    return params.ns;
}

- (void)setNs:(int)ns;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setNs:params.ns];
    params.ns = ns;
}

#pragma mark - Sample Details Tab

- (double)pd
{
    return params.pd;
}

- (void)setPd:(double)pd;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setPd:params.pd];
    params.pd = pd;
}

- (double)vgh
{
    return params.vgh;
}

- (void)setVgh:(double)vgh;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setVgh:params.vgh];
    params.vgh = vgh;
}

- (double)durn
{
    return params.durn;
}

- (void)setDurn:(double)durn;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setDurn:params.durn];
    params.durn = durn;
}

- (double)movementRate
{
    return params.rate;
}

- (void)setMovementRate:(double)rate;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setMovementRate:params.rate];
    params.rate = rate;
}

- (double)ps
{
    return params.ps;
}

- (void)setPs:(double)ps;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setPs:params.ps];
    params.ps = ps;
}

- (void)setElevationsAreSupplied:(BOOL)newElevationsAreSupplied
{
    [[[self undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithBool:elevationsAreSupplied]
                                                                 forKey:@"elevationsAreSupplied"];
    [[self undoManager] setActionName:@"Set Elevation Method"];
    elevationsAreSupplied = newElevationsAreSupplied;
    [controller updateTableColumns];
}

- (double)thh
{
    return params.thh;
}

- (void)setThh:(double)thh;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setThh:params.thh];
    params.thh = thh;
}

- (double)stt
{
    return params.stt;
}

- (void)setStt:(double)stt;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setStt:params.stt];
    params.stt = stt;
}

#pragma mark - Observations Tab

- (void)insertObject:(Observation *)o inObservationsAtIndex:(int)index
{
    //NSLog(@"adding %@ to %@", o, observations);
    // Add the inverse of this operation to the undo stack
    NSUndoManager *undo = [self undoManager];
    [[undo prepareWithInvocationTarget:self] removeObjectFromObservationsAtIndex:index];
    if (![undo isUndoing]) {
        [undo setActionName:@"Insert Observation"];
    }
    
    // Add the Person to the array
    [self startObservingObservation:o];
    [observations insertObject:o atIndex:index];
}

- (void)removeObjectFromObservationsAtIndex:(int)index
{
    Observation *o = [observations objectAtIndex:index];
    //NSLog(@"removing %@ from %@", o, observations);
    // Add the inverse of this operation to the undo stack
    NSUndoManager *undo = [self undoManager];
    [[undo prepareWithInvocationTarget:self] insertObject:o
                                    inObservationsAtIndex:index];
    if (![undo isUndoing]) {
        [undo setActionName:@"Delete Observation"];
    }
    
    [observations removeObjectAtIndex:index];
    [self stopObservingObservation:o];
}

- (void)startObservingObservation:(Observation *)o
{
    [o addObserver:self
        forKeyPath:@"distance"
           options:NSKeyValueObservingOptionOld
           context:NULL];

    [o addObserver:self
        forKeyPath:@"groupSize"
           options:NSKeyValueObservingOptionOld
           context:NULL];

    [o addObserver:self
        forKeyPath:@"angle"
           options:NSKeyValueObservingOptionOld
           context:NULL];

    [o addObserver:self
        forKeyPath:@"elevation"
           options:NSKeyValueObservingOptionOld
           context:NULL];
}

- (void)stopObservingObservation:(Observation *)o
{
    [o removeObserver:self forKeyPath:@"distance"];
    [o removeObserver:self forKeyPath:@"groupSize"];
    [o removeObserver:self forKeyPath:@"angle"];
    [o removeObserver:self forKeyPath:@"elevation"];
}

- (void)changeKeyPath:(NSString *)keyPath
             ofObject:(id)obj
              toValue:(id)newValue
{
    // setValue:forKeyPath: will cause the key-value observing method
    // to be called, which takes care of the undo stuff
    [obj setValue:newValue forKeyPath:keyPath];
}

- (void)observeValueForKeyPath:(NSString *)keyPath
                      ofObject:(id)object
                        change:(NSDictionary *)change
                       context:(void *)context
{
    NSUndoManager *undo = [self undoManager];
    id oldValue = [change objectForKey:NSKeyValueChangeOldKey];
    
    // NSNull objects are used to represent nil in a dictionary
    if (oldValue == [NSNull null]) {
        oldValue = nil;
    }
    //NSLog(@"oldValue = %@", oldValue);
    [[undo prepareWithInvocationTarget:self] changeKeyPath:keyPath
                                                  ofObject:object
                                                   toValue:oldValue];
    [undo setActionName:@"Edit Observation"];
}

#pragma mark - Options Tab

- (int)maxjb
{
    return params.maxjb;
}

- (void)setMaxjb:(int)maxjb
{
    [[[self undoManager] prepareWithInvocationTarget:self] setMaxjb:params.maxjb];
    params.maxjb = maxjb;
}

- (double)clint
{
    return params.clint;
}

- (void)setClint:(double)clint;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setClint:params.clint];
    params.clint = clint;
}

- (double)f0
{
    return params.f[0];
}

- (void)setF0:(double)f0;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setF0:params.f[0]];
    params.f[0] = f0;
}

- (double)f1
{
    return params.f[1];
}

- (void)setF1:(double)f1;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setF1:params.f[1]];
    params.f[1] = f1;
}

- (double)f2
{
    return params.f[2];
}

- (void)setF2:(double)f2;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setF2:params.f[2]];
    params.f[2] = f2;
}

- (double)f3
{
    return params.f[3];
}

- (void)setF3:(double)f3;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setF3:params.f[3]];
    params.f[3] = f3;
}

- (double)step0
{
    return params.step[0];
}

- (void)setStep0:(double)step0;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setStep0:params.step[0]];
    params.step[0] = step0;
}

- (double)step1
{
    return params.step[1];
}

- (void)setStep1:(double)step1;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setStep1:params.step[1]];
    params.step[1] = step1;
}

- (double)step2
{
    return params.step[2];
}

- (void)setStep2:(double)step2;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setStep2:params.step[2]];
    params.step[2] = step2;
}

- (int)iprint
{
    return params.iprint;
}

- (void)setIprint:(int)iprint
{
    [[[self undoManager] prepareWithInvocationTarget:self] setIprint:params.iprint];
    params.iprint = iprint;
}

- (int)ishow
{
    return params.ishow;
}

- (void)setIshow:(int)ishow;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setIshow:params.ishow];
    params.ishow = ishow;
}

- (int)jprint
{
    return params.jprint;
}

- (void)setJprint:(int)jprint;
{
    [[[self undoManager] prepareWithInvocationTarget:self] setJprint:params.jprint];
    params.jprint = jprint;
}

#pragma mark - Estimate Tab

- (int)nvals
{
    return [observations count];
}

- (int)currentIteration
{
    if (params.complete) {
        return -1;
    } else {
        return params.bootstrap;
    }
}

- (int)maxIteration
{
    return params.maxjb;
}

#pragma mark Implementation

- (BOOL)parseInputColumns:(NSString *)input
{
    // Approximate newlineCharacterSet on Tiger
    NSCharacterSet *newlineCharacterSet;
    SInt32 systemVersion = 0;
    OSStatus err = Gestalt(gestaltSystemVersion, &systemVersion);
    if (systemVersion >= 0x1050) {
        newlineCharacterSet = [NSCharacterSet newlineCharacterSet];
    } else {
        newlineCharacterSet = [NSCharacterSet characterSetWithCharactersInString:[NSString stringWithFormat:@"%c%c", 10, 13]];
    }

	NSScanner* scanner = [NSScanner scannerWithString:input];

    // Read the first line as the header
    NSString* firstLine;
    if (![scanner scanUpToCharactersFromSet:newlineCharacterSet intoString:&firstLine]) return NO;
    // Trim any opening or closing quote
    header = [firstLine stringByTrimmingCharactersInSet:[NSCharacterSet characterSetWithCharactersInString:@"'"]];
    [header retain];

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
    if ([scanner scanString:@"," intoString:nil]) {
        if (![scanner scanDouble:&params.thh]) return NO;
    }

    int i = 0;
    unsigned bottomParams;
    while (i < MAX_OBSERVATIONS) {
        unsigned potentialLine = [scanner scanLocation];
        BOOL potentialElevations = NO;
        if ([scanner scanDouble:(params.r + i)]) {
            bottomParams = potentialLine; /* potentialLine really was a new line */
            if (potentialElevations) {
                elevationsAreSupplied = YES;
            }
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
            potentialElevations = YES;
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

	if (![scanner scanDouble:&params.step[NUM_SHAPE_PARAMS - 1]]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
    // if there's another value then that last was a leftover imv value
	if (![scanner scanDouble:&params.step[NUM_SHAPE_PARAMS - 1]]) return YES;
	if (![scanner scanString:@"," intoString:nil]) return YES; /* Skip the comma separator */
    // and if there's a final value then that last was a leftover r3s value
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
	if (![scanner scanInt:nil]) return NO; // was imv
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
	if (![scanner scanDouble:nil]) return NO; // was r3s
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
	if (![scanner scanInt:nil]) return NO; // was imv
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
	if (![scanner scanDouble:nil]) return NO; // was r3s
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

// From http://www.red-sweater.com/blog/229/
- (void)endEditingSavingCurrentResponder
{
    // Save the current first responder, respecting the fact
    // that it might conceptually be the delegate of the
    // field editor that is \u201cfirst responder.\u201d
    NSWindow *oMainDocumentWindow = [self windowForSheet];
    //oMainDocumentWindow = [(NSWindowController*)[[self windowControllers] objectAtIndex:0] window];
    id oldFirstResponder = [oMainDocumentWindow firstResponder];
    if ((oldFirstResponder != nil) &&
        [oldFirstResponder isKindOfClass:[NSTextView class]] &&
        [(NSTextView*)oldFirstResponder isFieldEditor])
    {
        // A field editor\u2019s delegate is the view we\u2019re editing
        oldFirstResponder = [oldFirstResponder delegate];
        if ([oldFirstResponder isKindOfClass:[NSResponder class]] == NO)
        {
            // Eh \u2026 we\u2019d better back off if
            // this thing isn\u2019t a responder at all
            oldFirstResponder = nil;
        }
    }

    // Gracefully end all editing in our window (from Erik Buck).
    // This will cause the user\u2019s changes to be committed.
    if([oMainDocumentWindow makeFirstResponder:oMainDocumentWindow])
    {
        // All editing is now ended and delegate messages sent etc.
    }
    else
    {
        // For some reason the text object being edited will
        // not resign first responder status so force an
        /// end to editing anyway
        [oMainDocumentWindow endEditingFor:nil];
    }

    // If we had a first responder before, restore it
    if (oldFirstResponder != nil)
    {
        [oMainDocumentWindow makeFirstResponder:oldFirstResponder];
    }
}

- (NSString *)resultsFileName
{
    NSString *datasetBaseName = [[self fileName] stringByDeletingPathExtension];
    return [NSString stringWithFormat:@"%@ %d.results", datasetBaseName, currentResultsIndex];
}
- (NSString *)graphDataFileName
{
    NSString *datasetBaseName = [[self fileName] stringByDeletingPathExtension];
    return [NSString stringWithFormat:@"%@ %d.graphData", datasetBaseName, currentResultsIndex];
}

- (IBAction)viewCurrentResults:(id)sender
{
    NSString *resultsFile = [self resultsFileName];
    if ([[NSFileManager defaultManager] fileExistsAtPath:resultsFile]) {
        LSLaunchURLSpec viewerSpec;
        OSStatus err = LSFindApplicationForInfo(kLSUnknownCreator, CFSTR("com.macromates.textmate"), NULL, NULL, &viewerSpec.appURL);
        if (err != noErr)
            err = LSFindApplicationForInfo(kLSUnknownCreator, CFSTR("com.barebones.textwrangler"), NULL, NULL, &viewerSpec.appURL);
        if (err != noErr)
            viewerSpec.appURL = (CFURLRef)[NSURL fileURLWithPath:@"/Applications/TextEdit.app"];
        viewerSpec.itemURLs = (CFArrayRef)[NSArray arrayWithObject:[NSURL fileURLWithPath:resultsFile]];
        viewerSpec.passThruParams = NULL;
        viewerSpec.launchFlags = kLSLaunchDefaults | kLSLaunchDontAddToRecents;
        viewerSpec.asyncRefCon = NULL;
        
        err = LSOpenFromURLSpec(&viewerSpec, NULL);
    }
}

- (IBAction)graphCurrentResults:(id)sender
{
    NSString *graphFile = [self graphDataFileName];
    if ([[NSFileManager defaultManager] fileExistsAtPath:graphFile]) {
        CFURLRef plotURL;
        OSStatus err = LSFindApplicationForInfo(kLSUnknownCreator, CFSTR("de.micw.plot"), NULL, NULL, &plotURL);
        if (err == noErr) {
            CFStringRef plotPath = CFURLCopyFileSystemPath(plotURL, kCFURLPOSIXPathStyle);
            char plotPathC[PATH_MAX + 1];
            if (CFStringGetCString(plotPath, plotPathC, PATH_MAX, kCFStringEncodingUTF8)) {
                strncat(plotPathC, "/Contents/MacOS/Plot", PATH_MAX);
                plotPathC[PATH_MAX] = 0;

                NSString *graphMacro = [[NSBundle mainBundle] pathForResource:@"Graph Results Plot macro" ofType:nil];
                if (fork() == 0) {
                    /* Child process */
                    execl(plotPathC,
                          plotPathC,
                          "-m",
                          [graphMacro UTF8String],
                          "-i",
                          [graphFile UTF8String],
                          NULL);
                    perror("exec Plot.app failed");
                }
            }
            CFRelease(plotURL);
            CFRelease(plotPath);
        }
    }
}

const double PI = 3.14159254;
double deg2rad(double deg) {
    return deg / 180 * PI;
}

- (void)launchCalculation
{
    [self endEditingSavingCurrentResponder];
    params.complete = 0;
    [self setValue:@"" forKey:@"completeMsg"];

    // Trim any empty observations
    [self willChangeValueForKey:@"nvals"];
    for (int i = 0; i < [observations count]; ) {
        int nsize = [[observations objectAtIndex:i] groupSize];
        if (nsize == 0) {
            [self removeObjectFromObservationsAtIndex:i];
        } else {
            ++i;
        }
    }
    [[self undoManager] setActionName:@"Trim empty observations"];
    [self didChangeValueForKey:@"nvals"];

    // Compute THH and LTMAX
    int elevationCount = 0;
    double sumOfSquaredElevations = 0;
    params.ltmax = 0;
    params.nvals = [observations count];
    for (int i = 0; i < params.nvals; i++) {
        Observation *observation = [observations objectAtIndex:i];
        params.r[i] = [observation distance];
        params.angle[i] = [observation angle];
        params.nsize[i] = [observation groupSize];
        
        double elevationAngle = [observation elevation];
        if (elevationAngle != -1) {
            double horizontalDistance = params.r[i];
            if (params.km == 2) horizontalDistance *= 1000;
            elevationCount++;
            double elevation = horizontalDistance * tan(deg2rad(elevationAngle));
            sumOfSquaredElevations += elevation * elevation;
        }
        if (params.r[i] > params.ltmax) {
            params.ltmax = params.r[i];
        }
    }
    if (elevationsAreSupplied && elevationCount) {
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
    
    NSFileManager *fileMgr = [NSFileManager defaultManager];
    do {
        ++currentResultsIndex;
    } while ([fileMgr fileExistsAtPath:[self resultsFileName]]
             || [fileMgr fileExistsAtPath:[self graphDataFileName]]);

    SInt32 systemVersion = 0;
    OSStatus err = Gestalt(gestaltSystemVersion, &systemVersion);
    if (systemVersion >= 0x1040) {
        [NSThread detachNewThreadSelector:@selector(calculationThread:) toTarget:self withObject:nil];
    } else {
        [self calculationWork];
    }
}

- (void)calculationThread:(id)ignored
{
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    [self calculationWork];
    [pool release];
}

- (void)calculationWork
{
    NSString *outFile = [self resultsFileName];
    NSString *graphFile = [self graphDataFileName];
    NSLog(@"Calculating to %@ and %@", outFile, graphFile);
    const char * hdr = [header UTF8String];
    const char * out = [outFile UTF8String];
    const char * graph = [graphFile UTF8String];

    calculate_density(&params, hdr, out, graph, strlen(hdr), strlen(out), strlen(graph));

    [completeMsg release];
    if (params.complete) {
        completeMsg = [NSString stringWithFormat:@"Calculation complete %c%cDensity estimate: %.3g  Standard error: %.2g%c%cDetailed results in %@ and %@",
            10, 10, params.results.estden, params.results.sden, 10, 10, outFile, [graphFile lastPathComponent]];
    } else {
        completeMsg = [NSString stringWithFormat:@"Calculation failed: Detailed results in %@", outFile];
        params.complete = 1;
    }
    [completeMsg retain];

    GraphData *graphData = [[GraphData alloc] initForURL:[NSURL fileURLWithPath:graphFile] withContents:&params.results];
    [[NSDocumentController sharedDocumentController] addDocument:graphData];
    [graphData makeWindowControllers];
    //[graphData showWindows];
}

- (calc_results *)internalResults
{
    return &params.results;
}
@end
