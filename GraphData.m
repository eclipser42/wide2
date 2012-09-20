//
//  GraphData.m
//
//  Created by James Clough on 25/11/2011.
//

#import "GraphData.h"
#import "GRLineDataSet.h"

@interface DataSetContent : NSObject
{
    double distance;
    int groupSize;
    double angle;
    double elevation;
}
@end

@implementation DataSetContent
@end

@implementation GraphData
- (void) dealloc
{
    [modelData release];
    [observedData release];
    [super dealloc];
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)type
{
    NSString *string = [[NSString alloc] initWithData:data encoding:NSASCIIStringEncoding];
    NSUInteger length = [string length];
    NSUInteger paraStart = 0, paraEnd = 0, contentsEnd = 0;
    // discard the first line
    [string getParagraphStart:&paraStart end:&paraEnd
                  contentsEnd:&contentsEnd forRange:NSMakeRange(paraEnd, 0)];
    // and parse remaining lines
    distance = [[NSMutableArray alloc] init];
    model = [[NSMutableArray alloc] init];
    observed = [[NSMutableArray alloc] init];
    while (paraEnd < length) {
        [string getParagraphStart:&paraStart end:&paraEnd
                      contentsEnd:&contentsEnd forRange:NSMakeRange(paraEnd, 0)];
        NSString *currentLine = [string substringWithRange:NSMakeRange(paraStart, contentsEnd - paraStart)];
        NSScanner *sc = [NSScanner scannerWithString:currentLine];
        double d, o, m;
        [sc scanDouble:&d];
        [sc scanDouble:&m];
        [sc scanDouble:&o];
        [distance addObject:[NSNumber numberWithDouble:d]];
        [model addObject:[NSNumber numberWithDouble:m]];
        [observed addObject:[NSNumber numberWithDouble:o]];
    }
    //NSLog(@"calc %@", model);
    //NSLog(@"obsd %@", observed);
    [string release];
    return TRUE;
}
- (NSString *)windowNibName
{
    return @"GraphData";
}
- (void) awakeFromNib
{
    NSLog(@"Hello");
	[chartView setProperty: [NSNumber numberWithInt: 0] forKey: GRChartDrawBackground];

	// Force the Y-axis to display from zero
	GRAxes * axes = [chartView axes];
	[axes setProperty: [NSNumber numberWithInt: 0] forKey: @"GRAxesYPlotMin"];
	[axes setProperty: [NSNumber numberWithInt: 1] forKey: @"GRAxesFixedYPlotMin"];

	//modelData = [[GRLineDataSet alloc] initWithOwnerChart: chartView];
	[chartView addDataSet: modelData loadData: YES];
	//observedData = [[GRLineDataSet alloc] initWithOwnerChart: chartView];
    [observedData setProperty:[NSArray arrayWithObjects:[NSNumber numberWithInt:6],
                               [NSNumber numberWithInt:4],
                               nil] forKey:GRDataSetPlotLineDashPattern];
    NSNumber *dp = [observedData propertyForKey:GRDataSetPlotLineWidth];//2
    [observedData setProperty:[NSNumber numberWithDouble:0.8] forKey:GRDataSetPlotLineWidth];
	[chartView addDataSet: observedData loadData: YES];
}
// Delegate methods for GRChartView
// Or maybe these are the dataSource methods

- (NSInteger) chart: (GRChartView *) chartView numberOfElementsForDataSet: (GRDataSet *) dataSet
{
	return [model count];
}

- (double) chart: (GRChartView *) chartView yValueForDataSet: (GRDataSet *) dataSet element: (NSInteger) element
{
    if (dataSet == modelData)
        return [[model objectAtIndex: element] doubleValue];
    else
        return [[observed objectAtIndex: element] doubleValue];
}

- (NSColor *) chart: (GRChartView *) chartView colorForDataSet: (GRDataSet *) dataSet element: (NSInteger) element
{
	return [NSColor blueColor];
}
@end
