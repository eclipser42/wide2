/* DensityController */

#import <Cocoa/Cocoa.h>
#import "Density.h"

@interface DensityController : NSObject
{
    IBOutlet Density *document;
}
- (IBAction)calculate:(id)sender;
@end
