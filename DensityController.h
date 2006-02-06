/* DensityController */

#import <Cocoa/Cocoa.h>
#import "Density.h"

@interface DensityController : NSObject
{
    IBOutlet Density *document;
    IBOutlet NSButton *calculateButton;
    IBOutlet NSProgressIndicator *progressBar;
}
- (IBAction)calculate:(id)sender;
- (void)updateProgress:(NSTimer*)timer;
- (void)finishCalculation;
@end
