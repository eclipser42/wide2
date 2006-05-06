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
    [self setKeys:[NSArray arrayWithObjects:@"censusType", @"detectionUnit", nil]
        triggerChangeNotificationsForDependentKey:@"estdenUnitString"];
}

- (id)init
{
    self = [super init];
    if (self) {
        detectionUnitString = @"m";
    }
    return self;
}

- (int)censusType
{
    if ([document ifx])
    {
        return 2;
    }
    else if ([document iry])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

- (void)setCensusType:(int)censusType
{
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

- (int)transectUnit
{
    return ([document km] > 0);
}

- (void)setTransectUnit:(int)transectUnit
{
    [self setDetectionUnit:transectUnit];
    //[document setValue:[NSNumber numberWithInt:transectUnit] forKey:@"km"];
}

- (BOOL)detectionUnitEnabled
{
    if ([document ifx])
    {
        return YES;
    }
    else
    {
        return ([document km] > 0);
    }
}

- (int)detectionUnit
{
    return ([document km] == 2) ? 2 : 1;
}

- (void)setDetectionUnit:(int)detectionUnit
{
    NSLog(@"setting detectionUnit(KM) to %d", detectionUnit);
    [document setKm:detectionUnit];
    switch (detectionUnit)
    {
        default:
            detectionUnitString = @"m";
            estdenUnitString = @"ind/ha";
            break;
        case 2:
            detectionUnitString = @"km";
            break;
    }
    if (detectionUnit > 0)
    {
        estdenUnitString = @"ind/km2";
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
        NSLog(@"updating to %g", iteration);
        [progressBar setDoubleValue:iteration];
    } else {
        [timer invalidate];
        [progressBar setHidden:YES];
        [document setValue:[document valueForKey:@"completeMsg"] forKey:@"completeMsg"];
        [calculateButton setHidden:NO];
        NSLog(@"done updating");
    }
}

@end
