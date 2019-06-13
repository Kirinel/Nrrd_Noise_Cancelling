#include <stdio.h>
#include <teem/nrrd.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <math.h>

#define BLURY  80
#define BLURX  80
#define PASS1  5.5
#define PASS2  1.9

void mean_sd_from_list(char *basePath) {
	char meanname[]="mean.nrrd", sdname[] = "sd.nrrd" ,*err;
    DIR *dir;
    struct dirent *ptr;
    char base[1024];
 
    int first = 1;
    Nrrd* nmean = nrrdNew();
    Nrrd* nstd = nrrdNew();
    int count = 0;
    unsigned long sx, sy, sc, sz;
    double* meanVal;
    double* sdVal;

    if ((dir=opendir(basePath)) == NULL) {
         perror("Open dir error...");
         exit(1);
    }
    while ((ptr=readdir(dir)) != NULL) {
    	// read the file and modify the file name string
    	if(*ptr->d_name == '.') continue;
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0) ///current dir OR parrent dir
            continue;
        printf("d_name:%s/%s\n",basePath,ptr->d_name);
        memset(base, '\0', sizeof(base));
        strcpy(base, basePath);
        strcat(base, "/");
        strcat(base, ptr->d_name);
        count++;
        char* filename = base;

        // load the file into an empty Nrrd
        Nrrd* nin = nrrdNew();
        if (nrrdLoad(nin, filename, NULL)) {
		    err = biffGetDone(NRRD);
		    fprintf(stderr, "trouble reading \"%s\":\n%s", filename, err);
		    free(err);
		    return;
		}

		// if it's the first file, initialize the output nrrd
		if(first) {
			sx = nin->axis[0].size, sy = nin->axis[1].size, sc = nin->axis[2].size, sz = nin->axis[3].size;
			const size_t sizes[3] = {sx, sy, sc};
			nrrdAlloc_nva(nmean, nrrdTypeDouble, 3, sizes);
			first = 0;
		}

		//calculate the sum value(not the mean)
		meanVal= (double*)nmean->data;
  		unsigned short* data = (unsigned short*)(nin->data);
  		for(int i = 0; i < sx; i++) {
		 	for(int j = 0; j < sy; j++) {
		  		for(int c = 0; c < sc; c++) {
		  			double sum = 0;
		  			for(int k = 0; k < sz; k++) {
		  				sum += (double)(data[i + sx * (j + sy * (c + sc * k))]);
		  			}
		  			if(first) meanVal[i + sx * (j + sy * c)] = sum;
		  			else meanVal[i + sx * (j + sy * c)] += sum;
		  		}
		  	}
		}
		nrrdNuke(nin);
    }
    closedir(dir);

    // divide by the size of the fourth dimension to get the mean value
    for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	  			meanVal[i + sx * (j + sy * c)] /= sz * count;
	  		}
	  	}
	}

	if (nrrdSave(meanname, nmean, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", meanname, err);
	    free(err);
	    return;
	}
	// finish calculating mean values

	// open the directory to read the files again to calculate the standard deviation
	if ((dir=opendir(basePath)) == NULL) {
         perror("Open dir error...");
         exit(1);
    }

    first = 1;
    while ((ptr=readdir(dir)) != NULL) {
    	// read the file and modify the name string
    	if(*ptr->d_name == '.') continue;
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir
            continue;
        printf("d_name:%s/%s\n",basePath,ptr->d_name);
        memset(base, '\0', sizeof(base));
        strcpy(base, basePath);
        strcat(base, "/");
        strcat(base, ptr->d_name);
        count++;
        // else if(ptr->d_type == 10)    ///link file
        //     printf("d_name:%s/%s\n",basePath,ptr->d_name);
        char* filename = base;

        // load the file into an empty Nrrd
        Nrrd* nin = nrrdNew();
        if (nrrdLoad(nin, filename, NULL)) {
		    err = biffGetDone(NRRD);
		    fprintf(stderr, "trouble reading \"%s\":\n%s", filename, err);
		    free(err);
		    return;
		}

		// if it is the first file, initialize the output nrrd
		if(first) {
			const size_t sizes[3] = {sx, sy, sc};
			nrrdAlloc_nva(nstd, nrrdTypeDouble, 3, sizes);
			first = 0;
		}

		// calculate the sum of the square of the difference
		sdVal = (double*)nstd->data;
  		unsigned short* data = (unsigned short*)(nin->data);
  		for(int i = 0; i < sx; i++) {
		 	for(int j = 0; j < sy; j++) {
		  		for(int c = 0; c < sc; c++) {
		  			double sum = 0;
		  			for(int k = 0; k < sz; k++) {
  						sum += ((double)(data[i + sx * (j + sy * (c + sc * k))]) - meanVal[i + sx * (j + sy * c)]) * ((double)(data[i + sx * (j + sy * (c + sc * k))]) - meanVal[i + sx * (j + sy * c)]);
		  			}
		  			if(first) sdVal[i + sx * (j + sy * c)] = sum;
		  			else sdVal[i + sx * (j + sy * c)] += sum;
		  		}
		  	}
		}
		nrrdNuke(nin);
    }

    closedir(dir);

    // divide by the size of the fourth dimension and take square root to get standard deviation
    for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	  			sdVal[i + sx * (j + sy * c)] = sqrt(sdVal[i + sx * (j + sy * c)] / (sz * count));
	  		}
	  	}
	}

	if (nrrdSave(sdname, nstd, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", sdname, err);
	    free(err);
	    return;
	}

    return;
}



int nrrd_scale(Nrrd *nout, Nrrd *nin, double scale) {
	unsigned long sx, sy, sc;
	sx = nin->axis[0].size, sy = nin->axis[1].size, sc = nin->axis[2].size;
	const size_t sizes[3] = {sx, sy, sc};
	nrrdAlloc_nva(nout, nrrdTypeDouble, 3, sizes);
	double *data1 = (double*)nin->data, *odata = (double*)nout->data;
	for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	 			int pos = i + sx * (j + sy * c);
	 			odata[pos] = data1[pos] * scale;
	  		}
	  	}
	}
	return 0;
}

int nrrd_1op(Nrrd *nout, Nrrd *nin, char nop) {
	char *err;
	int op = ('r' == nop ? nrrdUnaryOpReciprocal :
				('e' == nop ? nrrdUnaryOpExp :
					('l' == nop ? nrrdUnaryOpLog : 0)));

	if (!op) {
		fprintf(stderr, "wrong operator\n");
	 	return 1;
	}
	if(nrrdArithUnaryOp(nout, op, nin)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "error doing unary operation:\n%s", err);
	    free(err);
	    return 1;
	}
	// printf("computed 1op: %zu\n", nout->axis[0].size);
	return 0;
}

int nrrd_2op(Nrrd *nout, Nrrd *nin1, Nrrd *nin2, char nop) {
	char * err;
	int op = ('p' == nop ? nrrdBinaryOpAdd :
				('m' == nop ? nrrdBinaryOpSubtract :
					('l' == nop ? nrrdUnaryOpLog : 0)));

	if (!op) {
		fprintf(stderr, "wrong operator\n");
	 	return 1;
	}
	if(nrrdArithBinaryOp(nout, op, nin1, nin2)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "error doing binary operation:\n%s", err);
	    free(err);
	    return 1;
	}
	return 0;
}


int nrrd_blur(Nrrd *nout, Nrrd *nin, const double bx, const double by, double cutoff) {
	const double blur[2] = {bx, by};
	NrrdResampleContext *rsmc = nrrdResampleContextNew();
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



int main() {
    // char basePath[1000];

    // ///get the current absoulte path
    // memset(basePath,'\0',sizeof(basePath));
    // getcwd(basePath, 999);
    // printf("the current dir is : %s\n",basePath);

    // ///get the file list
    // memset(basePath,'\0',sizeof(basePath));
    // strcpy(basePath,"./data");
    // mean_sd_from_list(basePath);

    
    char *err;
    char *meanPath = "./mean_sd_nrrd/mean.nrrd";
    char *sdPath = "./mean_sd_nrrd/sd.nrrd";
    Nrrd *nmean = nrrdNew();
    if (nrrdLoad(nmean, meanPath, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", meanPath, err);
		nrrdNuke(nmean);
		free(err);
		return 1;
	}


	/*
	unu resample -b weight -i mean.nrrd  -s = x1 = -k gauss:$BLURY,4 -o blury.nrrd;
   unu resample -b weight -i blury.nrrd -s x1 = = -k gauss:$PASS1,4 | unu 2op - blury.nrrd - | unu 2op - mean.nrrd - -o mean-loX.nrrd
   unu resample -b weight -i mean.nrrd -s x1 = = -k gauss:$BLURX,4 -o blurx.nrrd;
   
   doint this: unu resample -b weight -i blurx.nrrd -s = x1 = -k gauss:$PASS1,4 | unu 2op - blurx.nrrd - | unu 2op - mean-loX.nrrd - -o mean-loXloY.nrrd
   unu resample -b weight -i mean-loXloY.nrrd -s x1 x1 = -k gauss:$PASS2,4 -o mean-lo.nrrd
   unu 2op - mean.nrrd mean-lo.nrrd -o mean-noise.nrrd
	*/


// unu resample -b weight -i mean.nrrd  -s = x1 = -k gauss:$BLURY,4 -o blury.nrrd;
   	Nrrd *meanBlury = nrrdNew();
   	char *meanBluryPath = "./results/meanBlury.nrrd";
	if(nrrd_blur(meanBlury, nmean, 0, BLURY, 4)
		| nrrdSave(meanBluryPath, meanBlury, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble doing bluring or saving %s :\n%s", meanBluryPath, err);
	    nrrdNuke(meanBlury);
	    free(err);
	    return 1;
	}

// unu resample -b weight -i blury.nrrd -s x1 = = -k gauss:$PASS1,4 | unu 2op - blury.nrrd - | unu 2op - mean.nrrd - -o mean-loX.nrrd
	Nrrd *temp1 = nrrdNew();
	Nrrd *temp2 = nrrdNew();
	Nrrd *mean_loX = nrrdNew();
	char *meanloxPath = "./results/mean-loX.nrrd";

	if(nrrd_blur(temp1, meanBlury, PASS1, 0, 4)
		| nrrd_2op(temp2, meanBlury, temp1, 'm')
		| nrrd_2op(mean_loX, nmean, temp2, 'm')
		| nrrdSave(meanloxPath, mean_loX, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating meanloX or saving %s :\n%s", meanloxPath, err);
	    nrrdNuke(temp1);
	    nrrdNuke(temp2);
	    nrrdNuke(mean_loX);
	    free(err);
	    return 1;
	}

// unu resample -b weight -i mean.nrrd -s x1 = = -k gauss:$BLURX,4 -o blurx.nrrd;
	Nrrd *meanBlurx = nrrdNew();
	char *meanBlurxPath = "./results/meanBlurx.nrrd";
	if(nrrd_blur(meanBlurx, nmean, BLURX, 0, 4)
		| nrrdSave(meanBlurxPath, meanBlurx, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble doing bluring or saving %s :\n%s", meanBlurxPath, err);
	    nrrdNuke(meanBlurx);
	    free(err);
	    return 1;
	}

// unu resample -b weight -i blurx.nrrd -s = x1 = -k gauss:$PASS1,4 | unu 2op - blurx.nrrd - | unu 2op - mean-loX.nrrd - -o mean-loXloY.nrrd
	Nrrd *mean_loXloY = nrrdNew();
	char *meanloxloyPath = "./results/mean-loXloY.nrrd";
	if(nrrd_blur(temp1, meanBlurx, 0, PASS1, 4)
		| nrrd_2op(temp2, meanBlurx, temp1, 'm')
		| nrrd_2op(mean_loXloY, mean_loX, temp2, 'm')
		| nrrdSave(meanloxloyPath, mean_loXloY, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating meanloXloY or saving %s :\n%s", meanloxloyPath, err);
	    nrrdNuke(temp1);
	    nrrdNuke(temp2);
	    nrrdNuke(mean_loXloY);
	    free(err);
	    return 1;
	}

// unu resample -b weight -i mean-loXloY.nrrd -s x1 x1 = -k gauss:$PASS2,4 -o mean-lo.nrrd
	Nrrd *mean_lo = nrrdNew();
	char *meanloPath = "./results/mean-lo.nrrd";
	if(nrrd_blur(mean_lo, mean_loXloY, PASS2, PASS2, 4)
		| nrrdSave(meanloPath, mean_lo, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", meanloPath, err);
	    nrrdNuke(mean_lo);
	    free(err);
	    return 1;
	}

// unu 2op - mean.nrrd mean-lo.nrrd -o mean-noise.nrrd
	Nrrd *mean_noise = nrrdNew();
	char *meannoisePath = "./results/mean-noise.nrrd";
	if(nrrd_2op(mean_noise, nmean, mean_lo, 'm')
		| nrrdSave(meannoisePath, mean_noise, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", meannoisePath, err);
	    nrrdNuke(mean_noise);
	    free(err);
	    return 1;
	}



    Nrrd *nstd = nrrdNew();
	if (nrrdLoad(nstd, sdPath, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", sdPath, err);
		nrrdNuke(nstd);
		free(err);
		return 1;
	}

	
	// unu 1op log -o logstd.nrrd
	Nrrd *logstd = nrrdNew();
	char *logstdPath = "./results/logstd.nrrd";
	if(nrrd_1op(logstd, nstd, 'l')
		| nrrdSave(logstdPath, logstd, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble calculating logstd or saving \"%s\":\n%s", logstdPath, err);
		nrrdNuke(logstd);
		free(err);
		return 1;
	}

 //   	unu resample -b weight -i logstd.nrrd     -s = x1 = -k gauss:$BLURY,4 -o blury.nrrd
	Nrrd *logstdBlury = nrrdNew();
	char *logstdbluryPath = "./results/logBlury.nrrd";
	if(nrrd_blur(logstdBlury, logstd, 0, BLURY, 4)
		| nrrdSave(logstdbluryPath, logstdBlury, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble calculating logblury or saving \"%s\":\n%s", logstdbluryPath, err);
		nrrdNuke(logstdBlury);
		free(err);
		return 1;
	}

 //   	unu resample -b weight -i blury.nrrd -s x1 = = -k gauss:$PASS1,4 | unu 2op - blury.nrrd - | unu 2op - logstd.nrrd - -o logstd-loX.nrrd
	Nrrd *logstd_loX = nrrdNew();
	char *logstdloxPath = "./results/logstd-loX.nrrd";
	if(nrrd_blur(temp1, logstdBlury, PASS1, 0, 4)
		| nrrd_2op(temp2, logstdBlury, temp1, 'm')
		| nrrd_2op(logstd_loX, logstd, temp2, 'm')
		| nrrdSave(logstdloxPath, logstd_loX, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble calculating logstdlox or saving \"%s\":\n%s", logstdloxPath, err);
		nrrdNuke(logstd_loX);
		nrrdNuke(temp1);
		nrrdNuke(temp2);
		free(err);
		return 1;
	}

 //   	unu resample -b weight -i logstd.nrrd -s x1 = = -k gauss:$BLURX,4 -o blurx.nrrd
	Nrrd *logstdBlurx = nrrdNew();
	char *logstdblurxPath = "./results/logBlurx.nrrd";
	if(nrrd_blur(logstdBlurx, logstd, 0, BLURX, 4)
		| nrrdSave(logstdblurxPath, logstdBlurx, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble calculating logblurx or saving \"%s\":\n%s", logstdblurxPath, err);
		nrrdNuke(logstdBlurx);
		free(err);
		return 1;
	}

 //   	unu resample -b weight -i blurx.nrrd -s = x1 = -k gauss:$PASS1,4 | unu 2op - blurx.nrrd - | unu 2op - logstd-loX.nrrd - -o logstd-loXloY.nrrd
	Nrrd *logstd_loXloY = nrrdNew();
	char *logstdloxloyPath = "./results/logstd-loXloY.nrrd";
	if(nrrd_blur(temp1, logstdBlurx, 0, PASS1, 4)
		| nrrd_2op(temp2, logstdBlurx, temp1, 'm')
		| nrrd_2op(logstd_loXloY, logstd_loX, temp2, 'm')
		| nrrdSave(logstdloxloyPath, logstd_loXloY, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble calculating logstdloxloy or saving \"%s\":\n%s", logstdloxloyPath, err);
		nrrdNuke(logstd_loXloY);
		nrrdNuke(temp1);
		nrrdNuke(temp2);
		free(err);
		return 1;
	}

 //   	unu resample -b weight -i logstd-loXloY.nrrd -s x1 x1 = -k gauss:$PASS2,4 -o logstd-lo.nrrd
	Nrrd *logstd_lo = nrrdNew();
	char *logstdloPath = "./results/logstd-lo.nrrd";
	if(nrrd_blur(logstd_lo, logstd_loXloY, PASS2, PASS2, 4)
		| nrrdSave(logstdloPath, logstd_lo, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", logstdloPath, err);
	    nrrdNuke(logstd_lo);
	    free(err);
	    return 1;
	}

 //   	unu 2op - logstd.nrrd logstd-lo.nrrd -o logstd-noise.nrrd
	Nrrd *logstd_noise = nrrdNew();
	char *logstdnoisePath = "./results/logstd-noise.nrrd";
	if(nrrd_2op(logstd_noise, logstd, logstd_lo, 'm')
		| nrrdSave(logstdnoisePath, logstd_noise, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", logstdnoisePath, err);
	    nrrdNuke(logstd_noise);
	    free(err);
	    return 1;
	}

 //   	unu 1op exp -i logstd-noise.nrrd
 //   	unu 1op r -o scale-hi.nrrd
	Nrrd *scale_hi = nrrdNew();
	char *scalehiPath = "./results/scale-hi.nrrd";
	if(nrrd_1op(temp1, logstd_noise, 'e')
		| nrrd_1op(scale_hi, temp1, 'r')
		| nrrdSave(scalehiPath, scale_hi, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", scalehiPath, err);
	    nrrdNuke(scale_hi);
	    free(err);
	    return 1;
	}

	Nrrd *std_noise_canceled = nrrdNew();
	char *stdncanceledPath = "./results/std-noise-canceled.nrrd";
	if(nrrd_2op(std_noise_canceled, nstd, scale_hi, 'm')
		| nrrdSave(stdncanceledPath, std_noise_canceled, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble calculating or saving %s :\n%s", stdncanceledPath, err);
	    nrrdNuke(std_noise_canceled);
	    free(err);
	    return 1;
	}

	



	return 0;
}