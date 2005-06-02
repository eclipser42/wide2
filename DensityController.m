#import "DensityController.h"

@implementation DensityController

- (IBAction)calculate:(id)sender
{
    [calculateButton setEnabled:NO];
    [progressBar setMaxValue:[document maxIteration]];
    progressTimer = [NSTimer scheduledTimerWithTimeInterval:0.3 target:self selector:@selector(updateProgress:) userInfo:nil repeats:YES];
	[document calculate];
}

- (void)updateProgress:(id)ignored
{
    [progressBar setDoubleValue:[document currentIteration]];
}

- (void)finishCalculation
{
    [progressTimer fire];
    [progressTimer invalidate];
    progressTimer = nil;
    [calculateButton setEnabled:YES];
}
@end
