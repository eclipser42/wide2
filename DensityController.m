#import "GraphData.h"
#import "DensityController.h"

@implementation DensityController

+ (void)initialize {
    //[self setKeys:[NSArray arrayWithObjects:@"censusType", nil]
    //    triggerChangeNotificationsForDependentKey:@"dataSuppliedEnabled"];
    [self setKeys:[NSArray arrayWithObjects:@"detectionUnit", nil]
        triggerChangeNotificationsForDependentKey:@"transectUnit"];
    [self setKeys:[NSArray arrayWithObjects:@"censusType", @"transectUnit", nil]
        triggerChangeNotificationsForDependentKey:@"detectionUnitEnabled"];
    [self setKeys:[NSArray arrayWithObjects:@"transectUnit", nil]
        triggerChangeNotificationsForDependentKey:@"detectionUnit"];
    [self setKeys:[NSArray arrayWithObjects:@"censusType", @"detectionUnit", nil]
        triggerChangeNotificationsForDependentKey:@"detectionUnitString"];
    [self setKeys:[NSArray arrayWithObjects:@"censusType", @"xyz", nil]
        triggerChangeNotificationsForDependentKey:@"f3Enabled"];
}

- (id)init
{
    self = [super init];
    if (self) {
        maximumClassDistance = 500;
        detectionUnitString = @"m";
    }
    return self;
}

- (void)awakeFromNib
{
    [anglesColumn retain];
    [elevationsColumn retain];
    [document addObserver:self
               forKeyPath:@"kdt"
                  options:0
                  context:NULL];
    if ([document kdt] > 1) {
        [self willChangeValueForKey:@"maximumClassDistance"];
        maximumClassDistance = [document kdt];
        [self didChangeValueForKey:@"maximumClassDistance"];
    }
}

- (void)dealloc
{
    [document removeObserver:self forKeyPath:@"kdt"];
    [anglesColumn release];
    [elevationsColumn release];

    [super dealloc];
}

// Wrap [NSObject setValue:forKey:] in a name that can be captured by [NSUndoManager forwardInvocation:]
- (void)restoreValue:(id)value forKey:(NSString *)key
{
    [self setValue:value forKey:key];
}

- (void)observeValueForKeyPath:(NSString *)keyPath
                      ofObject:(id)object
                        change:(NSDictionary *)change
                       context:(void *)context
{
    // if kdt has changed in the document, then some of our keys' values may have too
    [self willChangeValueForKey:@"observationType"];
    [self willChangeValueForKey:@"classDistancesEnabled"];
    [self willChangeValueForKey:@"maximumClassDistance"];
    [self didChangeValueForKey:@"observationType"];
    [self didChangeValueForKey:@"classDistancesEnabled"];
    [self didChangeValueForKey:@"maximumClassDistance"];
}

#pragma mark Key-Value Coding

- (void)setNilValueForKey:(NSString *)key
{
    if ([key isEqualToString:@"maximumClassDistance"]) {
        [self setValue:[NSNumber numberWithInt:2] forKey:key];
    } else {
        [self setValue:[NSNumber numberWithInt:0] forKey:key];
    }
}

- (int)censusType
{
    if ([document ifx]) {
        return 2;
    } else if ([document iry]) {
        return 1;
    } else {
        return 0;
    }
}

- (void)setCensusType:(int)censusType
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithInt:[self censusType]] forKey:@"censusType"];
    switch (censusType)
    {
        default:
            [document setValue:[NSNumber numberWithInt:0] forKey:@"ifx"];
            [document setValue:[NSNumber numberWithInt:0] forKey:@"iry"];
            break;
        case 1:
            [document setValue:[NSNumber numberWithInt:0] forKey:@"ifx"];
            [document setValue:[NSNumber numberWithInt:1] forKey:@"iry"];
            break;
        case 2:
            [document setValue:[NSNumber numberWithInt:1] forKey:@"ifx"];
            [document setValue:[NSNumber numberWithInt:0] forKey:@"iry"];
            break;
    }
}

- (BOOL)distancesSuppliedEnabled
{
    return ([document iry] > 0);
}

- (int)observationType
{
    if ([document kdt] == 0) {
        return 0;
    } else if ([document kdt] == 1) {
        return 1;
    } else {
        // TODO: [doc kdt] should instead be in (int)maximumClassDistance
        maximumClassDistance = [document kdt];
        return 2;
    }
}

- (void)setObservationType:(int)observationType
{
    // rely on [document kdt] to manage undo and redo
    // also rely on KVObservation to update classDistancesEnabled
    if (observationType < 2)
        [document setValue:[NSNumber numberWithInt:observationType] forKey:@"kdt"];
    else
        [document setValue:[NSNumber numberWithInt:maximumClassDistance] forKey:@"kdt"];
}

- (BOOL)classDistancesEnabled
{
    return ([document kdt] > 1);
}

//- (int)minimumClassDistance
//- (void)setClassMinimum:(int)minimumClassDistance;

- (void)setMaximumClassDistance:(int)maximumDistance
{
    // rely on [document kdt] to manage undo and redo
    [document setValue:[NSNumber numberWithInt:maximumDistance] forKey:@"kdt"];
    maximumClassDistance = maximumDistance;
}

- (int)transectUnit
{
    return ([document km] > 0);
}

- (void)setTransectUnit:(int)transectUnit
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithInt:[self transectUnit]] forKey:@"transectUnit"];
    [self setDetectionUnit:transectUnit];
    //[document setValue:[NSNumber numberWithInt:transectUnit] forKey:@"km"];
}

- (BOOL)detectionUnitEnabled
{
    if ([document ifx]) {
        return YES;
    } else {
        return ([document km] > 0);
    }
}

- (int)detectionUnit
{
    return ([document km] == 2) ? 2 : 1;
}

- (void)setDetectionUnit:(int)detectionUnit
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithInt:[self detectionUnit]] forKey:@"detectionUnit"];
    [document setKm:detectionUnit];
    switch (detectionUnit)
    {
        default:
            detectionUnitString = @"m";
            break;
        case 2:
            detectionUnitString = @"km";
            break;
    }
}

- (BOOL)topographyDoesObscure
{
    return ([document ltmin] < 999);
}

- (void)setTopographyDoesObscure:(BOOL)topographyDoesObscure
{
    BOOL previousValue = [self topographyDoesObscure];
    if (topographyDoesObscure) {
        [document setLtmin:0];
    } else {
        // this has the unfortunate effect that an original change from ltmin=N (0<N<999) to doesNotObscure (ltmin=999) leaves N showing
        // in the greyed-out text field, but redoing that change leaves 0 showing. But undoing the re-done change does restore N.
        [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithDouble:[self minimumObscuringDistance]] forKey:@"minimumObscuringDistance"];
        [document setLtmin:999];
    }
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithBool:previousValue] forKey:@"topographyDoesObscure"];
}

- (double)minimumObscuringDistance
{
    if ([document ltmin] < 999) {
        return [document ltmin];
    } else {
        return 0;
    }
}

- (void)setMinimumObscuringDistance:(double)minimumObscuringDistance
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithDouble:[self minimumObscuringDistance]] forKey:@"minimumObscuringDistance"];
    if (minimumObscuringDistance < 999) {
        [document setLtmin:minimumObscuringDistance];
    }
}

-(void)updateTableColumns
{
    [anglesColumn setHidden:[document iry] != 1];
    [elevationsColumn setHidden:![[document valueForKey:@"elevationsAreSupplied"] boolValue]];
}

- (IBAction)calculate:(id)sender
{
    [calculateButton setHidden:YES];
    [progressBar setMaxValue:[document maxIteration]];
    [progressBar setDoubleValue:0];
    [progressBar setHidden:NO];
	[document launchCalculation];
    [NSTimer scheduledTimerWithTimeInterval:0.3 target:self selector:@selector(updateProgress:) userInfo:nil repeats:YES];
}

- (void)updateProgress:(NSTimer*)timer
{
    double iteration = [document currentIteration];
    if (iteration >= 0)
    {
        [progressBar setDoubleValue:iteration];
    } else {
        [timer invalidate];
        [progressBar setHidden:YES];
        [document setValue:[document valueForKey:@"completeMsg"] forKey:@"completeMsg"];
        [calculateButton setHidden:NO];

        GraphData *graphData = [[GraphData alloc] initForURL:[NSURL fileURLWithPath:[document graphDataFileName]] withContents:[document internalResults]];
        [[NSDocumentController sharedDocumentController] addDocument:graphData];
        [graphData makeWindowControllers];
        [graphData showWindows];
    }
}

@end
