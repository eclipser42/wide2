/* DensityController */

#import <Cocoa/Cocoa.h>
#import "Density.h"

@interface DensityController : NSObject
{
    IBOutlet Density *document;
    IBOutlet NSButton *calculateButton;
    IBOutlet NSProgressIndicator *progressBar;
    NSString *detectionUnitString;
    NSString *estdenUnitString;
}
- (int)censusType;
- (void)setCensusType:(int)censusType;
- (BOOL)distancesSuppliedEnabled;
- (int)transectUnit;
- (void)setTransectUnit:(int)transectUnit;
- (BOOL)detectionUnitEnabled;
- (int)detectionUnit;
- (void)setDetectionUnit:(int)detectionUnit;
- (IBAction)calculate:(id)sender;
- (void)updateProgress:(NSTimer*)timer;
@end
