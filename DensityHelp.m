#import "DensityHelp.h"

@implementation DensityHelp
- (IBAction)openManual:(id)sender {
    NSString *helpFile = [[NSBundle mainBundle] pathForResource:@"WildlifeDensity Manual" ofType:@"pdf"];
    NSURL *helpURL = [NSURL fileURLWithPath:helpFile];
    LSLaunchURLSpec launchSpec = {0};
    launchSpec.appURL = (CFURLRef)[NSURL fileURLWithPath:@"/Applications/Preview.app"];
    launchSpec.itemURLs = (CFArrayRef)[NSArray arrayWithObject:helpURL];
    launchSpec.launchFlags = kLSLaunchDontAddToRecents;
    OSStatus result = LSOpenFromURLSpec(&launchSpec, NULL);
    //OSStatus result = LSOpenCFURLRef((CFURLRef) helpURL, NULL);
    if (result) {
        NSLog(@"Opening WD help file %@ in Preview failed with error %d, retrying with default app", helpFile, result);
        launchSpec.appURL = NULL;
        launchSpec.launchFlags = kLSLaunchAndDisplayErrors | kLSLaunchDontAddToRecents;
        result = LSOpenFromURLSpec(&launchSpec, NULL);
        if (result) {
            NSLog(@"Opening WD help file %@ failed with error %d", helpFile, result);
        }
    }
}
@end
