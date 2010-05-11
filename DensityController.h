/* DensityController */

#import <Cocoa/Cocoa.h>
#import "Density.h"

@interface DensityController : NSObject
{
    IBOutlet Density *document;
    IBOutlet NSTableView *observationsTable;
    IBOutlet NSTableColumn *anglesColumn;
    IBOutlet NSTableColumn *elevationsColumn;
    IBOutlet NSButton *calculateButton;
    IBOutlet NSProgressIndicator *progressBar;
    NSString *detectionUnitString;
    int maximumClassDistance;
}
- (int)censusType;
- (void)setCensusType:(int)censusType;
- (BOOL)distancesSuppliedEnabled;
- (int)observationType;                         // KDT
- (void)setObservationType:(int)observationType;
- (BOOL)classDistancesEnabled;
//- (int)minimumClassDistance;                    // STT
//- (void)setClassMinimum:(int)minimumClassDistance;
- (void)setMaximumClassDistance:(int)maximumClassDistance;
- (int)transectUnit;
- (void)setTransectUnit:(int)transectUnit;
- (BOOL)detectionUnitEnabled;
- (int)detectionUnit;
- (void)setDetectionUnit:(int)detectionUnit;
- (BOOL)topographyDoesObscure;
- (void)setTopographyDoesObscure:(BOOL)topographyDoesObscure;
- (double)minimumObscuringDistance;
- (void)setMinimumObscuringDistance:(double)minimumObscuringDistance;
- (void)updateTableColumns;
- (IBAction)calculate:(id)sender;
- (void)updateProgress:(NSTimer*)timer;
@end
