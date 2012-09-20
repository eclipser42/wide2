//
//  GraphData.h
//
//  Created by James Clough on 25/11/2011.
//

#import <Cocoa/Cocoa.h>
#import "GRChartView.h"

@interface GraphData : NSDocument {
    bool internalResults;
    double densityEstimate;
    double standardError;
	NSArray * distance;
	NSArray * model;
	NSArray * observed;
	IBOutlet	GRChartView * chartView;
    GRDataSet * modelData;
    GRDataSet * observedData;
}

@end
