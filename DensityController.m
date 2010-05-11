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
    [observationsTable setDelegate:document];
    if ([document kdt] > 1) {
        [self willChangeValueForKey:@"maximumClassDistance"];
        maximumClassDistance = [document kdt];
        [self didChangeValueForKey:@"maximumClassDistance"];
    }
}

- (void)dealloc
{
    [anglesColumn release];
    [elevationsColumn release];

    [super dealloc];
}

// Wrap [NSObject setValue:forKey:] in a name that can be captured by [NSUndoManager forwardInvocation:]
- (void)restoreValue:(id)value forKey:(NSString *)key
{
    [self setValue:value forKey:key];
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
    [[[document undoManager] prepareWithInvocationTarget:document] restoreValue:[NSNumber numberWithInt:[document iry]] forKey:@"iry"];
    [[[document undoManager] prepareWithInvocationTarget:document] restoreValue:[NSNumber numberWithInt:[document ifx]] forKey:@"ifx"];
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
        maximumClassDistance = [document kdt];
        return 2;
    }
}

- (void)setObservationType:(int)observationType
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithInt:[self observationType]] forKey:@"observationType"];
    [self willChangeValueForKey:@"classDistancesEnabled"];
    if (observationType < 2)
        [document setValue:[NSNumber numberWithInt:observationType] forKey:@"kdt"];
    else
        [document setValue:[NSNumber numberWithInt:maximumClassDistance] forKey:@"kdt"];
    [self didChangeValueForKey:@"classDistancesEnabled"];
}

- (BOOL)classDistancesEnabled
{
    return ([document kdt] > 1);
}

//- (int)minimumClassDistance
//- (void)setClassMinimum:(int)minimumClassDistance;

- (void)setMaximumClassDistance:(int)maximumDistance
{
    [[[document undoManager] prepareWithInvocationTarget:self] restoreValue:[NSNumber numberWithInt:maximumDistance] forKey:@"maximumDistance"];
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

// perpenType is 1 for radial+angle perpendicular data; otherwise zero
// setting perpenType is only valid when ifx==0 and iry>0
/*- (int)perpenType
{
    return ([document iry] == 2);
}

- (void)setPerpenType:(int)perpenType
{
    NSAssert([document iry] > 0, @"wrong mode to set perpenType");
    [document setValue:[NSNumber numberWithInt:(perpenType + 1)] forKey:@"iry"];
}*/

-(void)updateTableColumns
{
    if ([document iry] == 1) {
        if ([anglesColumn tableView] == nil) {
            [observationsTable addTableColumn:anglesColumn];
        }
    } else {
        [observationsTable removeTableColumn:anglesColumn];
    }
    if ([[document valueForKey:@"elevationsAreSupplied"] boolValue]) {
        if ([elevationsColumn tableView] == nil) {
            [observationsTable addTableColumn:elevationsColumn];
        }
    } else {
        [observationsTable removeTableColumn:elevationsColumn];
    }
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
    if (iteration > 0)
    {
        [progressBar setDoubleValue:iteration];
    } else {
        [timer invalidate];
        [progressBar setHidden:YES];
        [document setValue:[document valueForKey:@"completeMsg"] forKey:@"completeMsg"];
        [calculateButton setHidden:NO];
    }
}

@end
