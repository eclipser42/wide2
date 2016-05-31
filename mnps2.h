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
    int num_intervals, km;
    double estden, sden;
} calc_results;

#define NUM_SHAPE_PARAMS 4
#define MAX_OBSERVATIONS 10000

typedef struct {
    int nvals,ifx,iry,ns,km,kdt;
    int iprint,jprint,ishow,maxjb;
    double durn,rate,clint,stt,dist,thh,ltmin,vgh,pd,ps;
    double enteredValue[NUM_SHAPE_PARAMS], enteredStep[NUM_SHAPE_PARAMS];
    int nsize[MAX_OBSERVATIONS];
    double r[MAX_OBSERVATIONS], angle[MAX_OBSERVATIONS];
    int complete, bootstrap;
    calc_results results;
} calc_params;

typedef struct {
    int numa, numo, ngroups, nclass;
    float distance[MAX_OBSERVATIONS]; // Corrected to metres if entered in km
    double clint, stt, estj, val[MAX_INTERVALS], f[NUM_SHAPE_PARAMS], step[NUM_SHAPE_PARAMS];
} search_params;

extern void calculate_density (calc_params *params,
                               const char *header, const char *outfile, const char *graphfile);
