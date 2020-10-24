//
//  aggregates.c
//  WildlifeDensity
//
//  Created by James Clough on 25-10-20.
//

#include "aggregates.h"
#include <math.h>
void add(aggregation *agg, double sample) {
    agg->_n += 1;
    agg->_sum += sample;
    agg->_sum2 += sample * sample;
}

int count(aggregation *agg) {
    return agg->_n;
}

double mean(aggregation *agg) {
    return agg->_sum / agg->_n;
}

double variance(aggregation *agg) {
    double mean = agg->_sum / agg->_n;
    return (agg->_sum2 / agg->_n) - mean * mean;
}

double population_sd(aggregation *agg) {
    return sqrt(variance(agg));
}

double sample_sd(aggregation *agg) {
    double ssd = agg->_n * agg->_sum2 - agg->_sum * agg->_sum;
    return sqrt(ssd / agg->_n / (agg->_n - 1));
}
