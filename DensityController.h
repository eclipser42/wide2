/* DensityController */

#import <Cocoa/Cocoa.h>
#import "Density.h"

@interface DensityController : NSObject
{
    IBOutlet Density *document;
    IBOutlet NSButton *calculateButton;
    IBOutlet NSProgressIndicator *progressBar;
    NSTimer *progressTimer;
}
- (IBAction)calculate:(id)sender;
- (void)updateProgress:(id)ignored;
- (void)finishCalculation;
@end
