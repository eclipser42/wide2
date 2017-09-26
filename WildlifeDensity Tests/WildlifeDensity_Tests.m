//
//  WildlifeDensity_Tests.m
//  WildlifeDensity Tests
//
//  Created by James Clough on 8-5-16.
//
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>
#import "Density.h"

@interface WildlifeDensity_Tests : XCTestCase

@end

@implementation WildlifeDensity_Tests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)oneDataset:(NSString *)dataset expectingResults:(NSString *)results andGraphData:(NSString *)graphData {
}
- (void)oneDataset:(NSString *)dataset {
    NSBundle *bundle = [NSBundle bundleForClass:[self class]];
    NSURL *datasetURL = [bundle URLForResource:dataset withExtension:@"WDdata"];
    Density *test = [[Density alloc] initWithContentsOfURL:datasetURL ofType:@"WildlifeDensity Dataset" error:nil];
    XCTAssertNotNil(test, @"Dataset %@ didn't load", dataset);

    [test launchCalculation];
    do {
        struct timespec rqtp = {0};
        rqtp.tv_nsec = 50 * 1000 * 1000;
        nanosleep(&rqtp, NULL);
    } while ([test currentIteration] >= 0);

    NSString *results = [bundle pathForResource:dataset ofType:@"results"];
    NSString *resultsComparison = [NSString stringWithFormat:@"diff -u \"%@\" \"%@\"", results, [test resultsFileName]];
    int diffExit = system([resultsComparison UTF8String]);
    XCTAssertEqual(0, diffExit);

    NSString *graphData = [bundle pathForResource:dataset ofType:@"graphData"];
    NSString *graphComparison = [NSString stringWithFormat:@"diff -u \"%@\" \"%@\"", graphData, [test graphDataFileName]];
    diffExit = system([graphComparison UTF8String]);
    XCTAssertEqual(0, diffExit);
}

- (void)testRadial {
    [self oneDataset:@"Run04096"];
}

- (void)testPerpendicular {
    [self oneDataset:@"Run04038y"];
}

- (void)testMean {
    [self oneDataset:@"Run05074"];
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
