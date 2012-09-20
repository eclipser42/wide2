//
//  Result.m
//  WildlifeDensity C
//
//  Created by Guest User on 19/09/12.
//
//

#import "Result.h"

@implementation Result
@synthesize testField;
- (id)init
{
    self = [super init];
    if (self) {
        ;
    }
    return self;
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
    return YES;
}

- (NSString *)windowNibName
{
    return @"Result";
}
- (void)windowControllerDidLoadNib:(NSWindowController *)windowController
{
    [testField setStringValue:[[self fileURL] path]];
}
@end
