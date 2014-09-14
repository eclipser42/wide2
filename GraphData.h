//
//  GraphData.h
//
//  Created by James Clough on 25/11/2011.
//

#import "mnps2.h"
#import <Cocoa/Cocoa.h>
#import "GRChartView.h"

@interface GraphData : NSDocument {
    bool internalResults;
    double densityEstimate;
    double standardError;
    NSString * densityUnitString;
    NSString * distanceAxisLabel;
	NSMutableArray * distance;
	NSMutableArray * model;
	NSMutableArray * observed;
	IBOutlet	GRChartView * chartView;
    GRDataSet * modelData;
    GRDataSet * observedData;
}
- (id)initForURL:(NSURL *)url withContents:(calc_results *)results;

@end
