/*
  compile with
  cc -o blur blur.c -Wall -O2 -I$TEEM_INCLUDE -L$TEEM_LIB -ltpz

  example use:
  ./blur -i mean.nrrd -b 0 80 -o blury.nrrd
 */
#include <teem/nrrd.h>
#include <teem/hest.h>

int
doblur(Nrrd *nout, Nrrd *nin, const double blur[2], double cutoff,
       airArray *mop) {

    NrrdResampleContext *rsmc = nrrdResampleContextNew();
    /* your code for freeing up rsmc when done may be different;
       one way or the other you need to call nrrdResampleContextNix */
    airMopAdd(mop, rsmc, (airMopper)nrrdResampleContextNix, airMopAlways);
    double parm[NRRD_KERNEL_PARMS_NUM];

    int E = 0;
    if (!E) E |= nrrdResampleNrrdSet(rsmc, nin);
    if (!E) E |= nrrdResampleDefaultCenterSet(rsmc, nrrdCenterCell);
    if (!E) E |= nrrdResampleBoundarySet(rsmc, nrrdBoundaryBleed);
    if (!E) E |= nrrdResampleTypeOutSet(rsmc, nin->type);
    if (!E) E |= nrrdResampleRenormalizeSet(rsmc, AIR_TRUE);
    // no resampling on slowest axis (2)
    if (!E) E |= nrrdResampleKernelSet(rsmc, 2, NULL, NULL);
    // same cutoff regardless of stdv
    parm[1] = cutoff;
    for (unsigned int axi=0; axi<2; axi++) {
        if (blur[axi] > 0) {
            if (!E) E |= nrrdResampleSamplesSet(rsmc, axi, nin->axis[0].size);
            if (!E) E |= nrrdResampleRangeFullSet(rsmc, axi);
            parm[0] = blur[axi];
            if (!E) E |= nrrdResampleKernelSet(rsmc, axi, nrrdKernelGaussian, parm);
        } else {
            if (!E) E |= nrrdResampleKernelSet(rsmc, axi, NULL, NULL);
        }
    }
    if (!E) E |= nrrdResampleExecute(rsmc, nout);
    return E;
}

int
main(int argc, const char *argv[]) {
    const char *me = argv[0];
    hestParm *hparm = hestParmNew();
    hestOpt *hopt=NULL;
    airArray *mop = airMopNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);

    Nrrd *nin = NULL;
    hestOptAdd(&hopt, "i", "nin", airTypeOther, 1, 1, &nin, NULL,
               "input image, assumed to be 3D array with slowest "
               "axis size 2", NULL, NULL, nrrdHestNrrd);
    double blur[2];
    hestOptAdd(&hopt, "b", "Bx By", airTypeDouble, 2, 2, blur, "10 10",
               "how much to blur in faster (Bx) and slower (By) "
               "axes of input; use 0 for no blurring at all");
    double cutoff;
    hestOptAdd(&hopt, "c", "cutoff", airTypeDouble, 1, 1, &cutoff, "4",
               "how many stdv to cut off blurring at");
    char *oFile;
    hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &oFile, NULL,
               "output image");

    hestParseOrDie(hopt, argc-1, argv+1, hparm, me,
                   "demo of nrrd blurring",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
    airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
    airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

    if (!( 3 == nin->dim && 2 == nin->axis[2].size )) {
        fprintf(stderr, "%s: expected 3-D X-by-Y-by-2 image\n", me);
        airMopError(mop); return 1;
    }
    Nrrd *nout = nrrdNew();
    airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);

    char *err;
    if (doblur(nout, nin, blur, cutoff, mop)) {
        airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
        fprintf(stderr, "%s: problem making output:\n%s\n", me, err);
        airMopError(mop); return 1;
    }
    if (nrrdSave(oFile, nout, NULL)) {
        airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
        fprintf(stderr, "%s: problem saving output:\n%s\n", me, err);
        airMopError(mop); return 1;
    }

    airMopOkay(mop);
    return 0;
}
