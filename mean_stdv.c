#include <stdio.h>
#include <teem/nrrd.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <math.h>

void mean_sd_from_list(char *basePath) {
	char meanname[]="mean.nrrd", sdname[] = "sd.nrrd" ,*err;
    DIR *dir;
    struct dirent *ptr;
    char base[1024];
 
    int first = 1;
    Nrrd* nmean = nrrdNew();
    Nrrd* nsd = nrrdNew();
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
			nrrdAlloc_nva(nsd, nrrdTypeDouble, 3, sizes);
			first = 0;
		}

		// calculate the sum of the square of the difference
		sdVal = (double*)nsd->data;
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

	if (nrrdSave(sdname, nsd, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", sdname, err);
	    free(err);
	    return;
	}

    return;
}



int nrrd_scale(Nrrd *nin, Nrrd *nout, double scale) {
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

int nrrd_1op(Nrrd *nin, Nrrd *nout, char nop) {
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

int nrrd_2op(Nrrd *nin1, Nrrd *nin2, Nrrd *nout, char nop) {
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
    char basePath[1000];

    ///get the current absoulte path
    memset(basePath,'\0',sizeof(basePath));
    getcwd(basePath, 999);
    printf("the current dir is : %s\n",basePath);

    ///get the file list
    memset(basePath,'\0',sizeof(basePath));
    strcpy(basePath,"./data");
    mean_sd_from_list(basePath);

	
    Nrrd *nin1, *nin2, *nout;
    char *err;
    char *in = "mean.nrrd";
    nin1 = nrrdNew();
    if (nrrdLoad(nin1, in, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in, err);
		nrrdNuke(nin1);
		free(err);
		return 1;
	}
    nin2 = nrrdNew();
	if (nrrdLoad(nin2, in, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in, err);
		nrrdNuke(nin2);
		free(err);
		return 1;
	}
	
    nout = nrrdNew();

	nrrd_2op(nin1, nin2, nout, 'm');


	if(nrrdSave("2op_test.nrrd", nout, NULL)) {
		err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing :\n%s", err);
	    free(err);
	    return 1;
	}
	nrrdNuke(nout);
	nrrdNuke(nin1);
	nrrdNuke(nin2);
	return 0;
}