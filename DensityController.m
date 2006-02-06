#import "DensityController.h"

@implementation DensityController

- (IBAction)calculate:(id)sender
{
    [calculateButton setHidden:YES];
    [progressBar setMaxValue:[document maxIteration]];
    [progressBar setDoubleValue:0];
    [progressBar setHidden:NO];
    [NSTimer scheduledTimerWithTimeInterval:0.3 target:self selector:@selector(updateProgress:) userInfo:nil repeats:YES];
	[document calculate];
}

- (void)updateProgress:(NSTimer*)timer
{
    double iteration = [document currentIteration];
    NSLog(@"updating to %g", iteration);
    [progressBar setDoubleValue:iteration];
    if (iteration >= [document maxIteration])
    {
        [timer invalidate];
        NSLog(@"done updating");
    }
}

- (void)finishCalculation
{
    [progressBar setHidden:YES];
    [calculateButton setHidden:NO];
}
@end
