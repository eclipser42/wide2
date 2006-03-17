#import "DensityController.h"

@implementation DensityController

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
