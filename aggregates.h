//
//  aggregates.h
//  WildlifeDensity
//
//  Created by James Clough on 25-10-20.
//

#ifndef aggregates_h
#define aggregates_h

typedef struct {
    // This format allows rapid accumulation and retains ~7 digits of precision as long as
    // the sum is no more than ~7 digits bigger than the standard deviation.
    int _n;
    double _sum, _sum2;
} aggregation;

void add(aggregation *agg, double sample);
int count(aggregation *agg);
// These aggregates are NaN unless count>0
double mean(aggregation *agg);
double variance(aggregation *agg);
double population_sd(aggregation *agg);
// This aggregate is NaN unless count>1
double sample_sd(aggregation *agg);

#endif /* aggregates_h */
