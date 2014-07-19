/*
 *  mnps2.h
 *  WildlifeDensity
 *
 *  Created by David Morgan on Mon Aug 02 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#define MAX_INTERVALS 80

typedef struct {
    double midpoints[MAX_INTERVALS], calcn[MAX_INTERVALS], obsdn[MAX_INTERVALS];
    int num_intervals;
    double estden, sden;
} calc_results;

#define NUM_SHAPE_PARAMS 4
#define MAX_OBSERVATIONS 10000

typedef struct {
    int nvals,ifx,iry,ns,km,kdt;
    int iprint,jprint,ishow,maxjb;
    double durn,rate,clint,stt,dist,thh,ltmin,vgh,pd,ps,f[NUM_SHAPE_PARAMS];
    double step[NUM_SHAPE_PARAMS], r[MAX_OBSERVATIONS];
    int nsize[MAX_OBSERVATIONS];
    double angle[MAX_OBSERVATIONS];
    int complete, bootstrap;
    calc_results results;
} calc_params;

extern void calculate_density (calc_params *params,
                              const char *header, const char *outfile, const char *graphfile);
