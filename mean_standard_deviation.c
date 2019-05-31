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
		  				// printf("i:%d j:%d c:%d k:%d\n", i, j, c, k);
  						sum += ((double)(data[i + sx * (j + sy * (c + sc * k))]) - meanVal[i + sx * (j + sy * c)]) * ((double)(data[i + sx * (j + sy * (c + sc * k))]) - meanVal[i + sx * (j + sy * c)]);
		  				// printf("current sum: %lf\n", sum);
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

void nrrd_pm(char* in1, char* in2, char* out, char op) {
	double sign = 1;
	if(op == '-') {
		sign = -1;
	}
	else if(op != '+' && op != '-') {
		fprintf(stderr, "wrong operator over %s and %s\n", in1, in2);
		return;
	}
	char *err;
	Nrrd *nin1 = nrrdNew();
	Nrrd *nin2 = nrrdNew();
	Nrrd *nout = nrrdNew();
	if (nrrdLoad(nin1, in1, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in1, err);
		free(err);
		return;
	}
	if (nrrdLoad(nin2, in22, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in2, err);
		free(err);
		return;
	}
	unsigned long sx, sy, sc, sz;
	sx = nin1->axis[0].size, sy = nin1->axis[1].size, sc = nin1->axis[2].size;
	if(sx != nin2->axis[0].size || sy != nin2->axis[1].size || sc != nin2->axis[2].size) {
		fprintf(stderr, "sizes of %s and %s don't match\n", in1, in2);
	}
	const size_t sizes[3] = {sx, sy, sc};
	nrrdAlloc_nva(nout, nrrdTypeDouble, 3, sizes);
	double *data1 = nin1->data, *data2 = nin2->data, odata = nout->data;
	for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	 			int pos = i + sx * (j + sy * c)
	 			odata[pos] = data1[pos] + data2[pos] * sign;
	  		}
	  	}
	}
	if (nrrdSave(out, nout, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", out, err);
	    free(err);
	    return;
	}
	return;
}

void nrrd_scale(char *in1, char *out, double scale) {
	char *err;
	Nrrd *nin1 = nrrdNew();
	Nrrd *nout = nrrdNew();
	if (nrrdLoad(nin1, in1, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in1, err);
		free(err);
		return;
	}
	unsigned long sx, sy, sc, sz;
	sx = nin1->axis[0].size, sy = nin1->axis[1].size, sc = nin1->axis[2].size;
	const size_t sizes[3] = {sx, sy, sc};
	nrrdAlloc_nva(nout, nrrdTypeDouble, 3, sizes);
	double *data1 = nin1->data, odata = nout->data;
	for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	 			int pos = i + sx * (j + sy * c)
	 			odata[pos] = data1[pos] * scale;
	  		}
	  	}
	}
	nrrdNuke(nin1);
	if (nrrdSave(out, nout, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", out, err);
	    free(err);
	    return;
	}
	return;
}

void nrrd_1op(char *in1, char *out, char op) {
	char *err;
	Nrrd *nin1 = nrrdNew();
	Nrrd *nout = nrrdNew();
	if (nrrdLoad(nin1, in1, NULL)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "trouble reading \"%s\":\n%s", in1, err);
		free(err);
		return;
	}
	unsigned long sx, sy, sc, sz;
	sx = nin1->axis[0].size, sy = nin1->axis[1].size, sc = nin1->axis[2].size;
	const size_t sizes[3] = {sx, sy, sc};
	nrrdAlloc_nva(nout, nrrdTypeDouble, 3, sizes);
	double *data1 = nin1->data, odata = nout->data;
	for(int i = 0; i < sx; i++) {
		for(int j = 0; j < sy; j++) {
	 		for(int c = 0; c < sc; c++) {
	 			int pos = i + sx * (j + sy * c)
	 			if(op == 'r') odata[pos] = 1 / data1[pos];
	 			else if(op == 'e') odata[pos] = exp(data1[pos]);
	 			else if(op == 'l') odata[pos] = log(data1[pos]);
	 			else {
	 				fprintf(stderr, "wrong operator on %s\n", in1);
	 				nrrdNuke(nin1);
	 				nrrdNuke(nout);
	 				return;
	 			}
	  		}
	  	}
	}
	if (nrrdSave(out, nout, NULL)) {
	    err = biffGetDone(NRRD);
	    fprintf(stderr, "trouble writing \"%s\":\n%s", out, err);
	    free(err);
	    return;
	}
	return;
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
	return 0;
}