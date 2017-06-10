//
//  PRNGTest.m
//  
//
//  Created by James Clough on 8-5-16.
//
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>
#import "mnps2.h"


@interface PRNGTest : XCTestCase
{
    int64_t seed;
}

@end

extern void srandnag (int iseed);
extern void randgen (double *result);

@implementation PRNGTest

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (double)parkAndMiller {
    const int64_t multiplier = 16807;
    const int64_t modulus = 2147483647;

    XCTAssertEqual(modulus, INT32_MAX);
    seed *= multiplier;
    seed %= modulus;

    return (double) seed / modulus;
}

- (void)testSeed0 {
    double oneResult;

    seed = 1;
    oneResult = [self parkAndMiller];
    XCTAssertEqualWithAccuracy(16807, oneResult * INT32_MAX, 1e-7);

    for (int i = 1; i < 9999; ++i) {
        oneResult = [self parkAndMiller];
    }
    [self parkAndMiller];
    XCTAssertEqual(1043618065, seed);

    int64_t tenThousandthResult = oneResult * INT32_MAX;
    int64_t nextResult = (tenThousandthResult * 16807) % INT32_MAX;
    XCTAssertEqual(1043618065, nextResult);
}

- (void)testSeed {
    double oneResult;

    srandnag(1);
    randgen(&oneResult);
    XCTAssertEqualWithAccuracy(16807, oneResult * INT32_MAX, 1e-7);

    for (int i = 1; i < 9999; ++i) {
        randgen(&oneResult);
    }
    int64_t tenThousandthResult = oneResult * INT32_MAX;
    int64_t nextResult = (tenThousandthResult * 16807) % INT32_MAX;
    XCTAssertEqual(1043618065, nextResult);
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
        [self parkAndMiller];
    }];
}

@end
