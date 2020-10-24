//
//  Aggregates_Tests.m
//  WildlifeDensity Tests
//
//  Created by James Clough on 25-10-20.
//

#import "aggregates.h"
#import <XCTest/XCTest.h>

@interface Aggregates_Tests : XCTestCase

@end

@implementation Aggregates_Tests

- (void)testExample {
    aggregation agg = {0};
    add(&agg, 3);
    add(&agg, 4);
    add(&agg, 5);
    XCTAssertEqual(count(&agg), 3);
    XCTAssertEqual(mean(&agg), 4);
}

- (void)testEmpty {
    aggregation empty = {0};
    XCTAssertTrue(isnan(mean(&empty)));
    XCTAssertTrue(isnan(variance(&empty)));
    XCTAssertTrue(isnan(population_sd(&empty)));
    XCTAssertTrue(isnan(sample_sd(&empty)));
    add(&empty, 0);
    XCTAssertFalse(isnan(mean(&empty)));
    XCTAssertFalse(isnan(variance(&empty)));
    XCTAssertFalse(isnan(population_sd(&empty)));
    XCTAssertTrue(isnan(sample_sd(&empty)));
    add(&empty, 0);
    XCTAssertFalse(isnan(sample_sd(&empty)));
}

- (void)testFemaleFulmars {
    aggregation agg = {0};
    add(&agg, 727.7);
    add(&agg, 1086.5);
    add(&agg, 1091);
    add(&agg, 1361.3);
    add(&agg, 1490.5);
    add(&agg, 1956.1);
    XCTAssertEqualWithAccuracy(count(&agg), 6, 0);
    XCTAssertEqualWithAccuracy(mean(&agg), 1285.5, 0.02);
    XCTAssertEqualWithAccuracy(variance(&agg), 886047.09/6, 0.5);
    XCTAssertEqualWithAccuracy(sample_sd(&agg), 420.96, 0.005);
}

- (void)testMaleFulmars {
    aggregation agg = {0};
    add(&agg, 525.8);
    add(&agg, 605.7);
    add(&agg, 843.3);
    add(&agg, 1195.5);
    add(&agg, 1945.6);
    add(&agg, 2135.6);
    add(&agg, 2308.7);
    add(&agg, 2950);
    XCTAssertEqualWithAccuracy(count(&agg), 8, 0);
    XCTAssertEqualWithAccuracy(mean(&agg), 1563.77, 0.02);
    XCTAssertEqualWithAccuracy(variance(&agg), 5599317.68/8, 0.02);
    XCTAssertEqualWithAccuracy(sample_sd(&agg), 894.37, 0.005);
}

@end
