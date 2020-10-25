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

@interface Density (TestAccess)
- (calc_params *)internalParameters;
@end
@implementation Density (TestAccess)
- (calc_params *)internalParameters {
    return &params;
}
@end

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

- (Density *)loadDataset:(NSString *)dataset {
    NSBundle *bundle = [NSBundle bundleForClass:[self class]];
    NSURL *datasetURL = [bundle URLForResource:dataset withExtension:@"WDdata"];
    Density *testData = [[Density alloc] initWithContentsOfURL:datasetURL ofType:@"WildlifeDensity Dataset" error:nil];
    XCTAssertNotNil(testData, @"Dataset %@ didn't load", dataset);
    return testData;
}

- (void)oneDataset:(Density *)dataset expectingResults:(NSString *)results andGraphData:(NSString *)graphData {
    [dataset launchCalculation];
    do {
        struct timespec rqtp = {0};
        rqtp.tv_nsec = 50 * 1000 * 1000;
        nanosleep(&rqtp, NULL);
    } while ([dataset currentIteration] >= 0);

    XCTAssertNotNil(results, @"Expected results didn't load");
    NSString *resultsComparison = [NSString stringWithFormat:@"diff -u \"%@\" \"%@\"", results, [dataset resultsFileName]];
    int diffExit = system([resultsComparison UTF8String]);
    XCTAssertEqual(0, diffExit);

    XCTAssertNotNil(graphData, @"Expected graphData didn't load");
    NSString *graphComparison = [NSString stringWithFormat:@"diff -u \"%@\" \"%@\"", graphData, [dataset graphDataFileName]];
    diffExit = system([graphComparison UTF8String]);
    XCTAssertEqual(0, diffExit);
}

- (void)oneDataset:(NSString *)dataset {
    Density *test = [self loadDataset:dataset];
    NSBundle *bundle = [NSBundle bundleForClass:[self class]];
    NSString *results = [bundle pathForResource:dataset ofType:@"results"];
    NSString *graphData = [bundle pathForResource:dataset ofType:@"graphData"];
    [self oneDataset:test expectingResults:results andGraphData:graphData];
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

- (void)testFixedPoint {
    NSString *filename = @"Run16070";
    Density *dataset = [self loadDataset:filename];
    calc_params *params = [dataset internalParameters];
    // This dataset needs 1000 evaluations for full precision, but 100 is adequate here
    params->maxjb = 100;

    NSBundle *bundle = [NSBundle bundleForClass:[self class]];
    NSString *manualResults = [bundle pathForResource:[filename stringByAppendingString:@"-manual a+c"] ofType:@"results"];
    NSString *autoResults = [bundle pathForResource:[filename stringByAppendingString:@"-auto a+c"] ofType:@"results"];
    NSString *manualGraphData = [bundle pathForResource:[filename stringByAppendingString:@"-manual a+c"] ofType:@"graphData"];
    NSString *autoGraphData = [bundle pathForResource:[filename stringByAppendingString:@"-auto a+c"] ofType:@"graphData"];

    for (double step0 = 0 ; step0 <= 1; ++step0) {
        for (double step1 = 0 ; step1 <= 1; ++step1) {
            params->enteredStep[0] = step0;
            params->enteredStep[1] = step1;
            NSString * expectedResults = (step0 && step1) ? manualResults : autoResults;
            NSString * expectedGraphData = (step0 && step1) ? manualGraphData : autoGraphData;
            [self oneDataset:dataset expectingResults:expectedResults andGraphData:expectedGraphData];
        }
    }
}

@end
