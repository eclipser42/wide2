/*
 *  mnps2.h
 *  WildlifeDensity
 *
 *  Created by David Morgan on Mon Aug 02 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#define NUM_SHAPE_PARAMS 4
#define MAX_OBSERVATIONS 10000

// This structure definition MUST be kept in sync with the definition in mnps2.f

typedef struct {
    int nvals,numa,numo,ltmin,ltmax,ifx,iry,ns,km,imv,kdt;
    int iprint,jprint,ishow,maxjb,it,iv;
    double clint,stt,dist,thh,r3s,vgh,pd,ps,f[NUM_SHAPE_PARAMS];
    double step[NUM_SHAPE_PARAMS], r[MAX_OBSERVATIONS];
    int nsize[MAX_OBSERVATIONS];
    double angle[MAX_OBSERVATIONS];
    int complete, bootstrap;
    double estden, sden;
} calc_params;

//extern void calculate_density(calc_params *params, const char *header,
//                              int headerlen);
extern void calculate_density(calc_params *params,
                              const char *header, const char *outfile, const char *graphfile,
                              int headerlen, int outfilelen, int graphfilelen);
