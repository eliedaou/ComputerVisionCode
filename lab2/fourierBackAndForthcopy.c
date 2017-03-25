// Alexander Leccese
// Computer Vission, Assignment 2, part 1
// Purpose: apply a 2-d fourier transform of an image and display its spectrum

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

int const n = 256;

four1(data,nn,isign)
float data[];
int nn,isign;
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

// get the minimum and maximum of the array
float getMaximumOfArray(float **inputtedArray);
float getMinimumOfArray(float **inputtedArray);

//find the magnitude
void findMagnitude(float **inputArray, float **outputArray);
// do fourier transform in 2d
void four2d(float **inArray);
// do inverse fourier in 2d
void four2dInverse(float **inArray);
// normalize the data based on 255 being the max
void normalize(float **inputArray, float **outputArray, float min, float max);
// function to output the files
void outputFile(char nameOfFile[256], float **dataArray);

// main function
void main(int argc, char *argv[]){
    int row, column;
    float **data, **spectrum, **newSpectrum, **realValuesOnly, **butter, **newButter, **butteredSpectrum, **centeredSpectrum;
    char outputFileName[256];
    unsigned char pixel;                // unsigned character to represent one pixel value
    FILE *input;
    int D0 = 10;
    int smoothness = 3;
  
  float minimum;           // minimum pixel value
  float maximum;          // maximum pixel value


    // 2-d array to hold original image data
    data =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        data[row] = (float *) calloc (n, sizeof (float));
    }
    
    // spectrum array to hold fourier transform image data
    spectrum =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        spectrum[row] = (float *) calloc (2*n, sizeof (float));
    }
    
    // new spectrum to hold the combined spectrum data
    newSpectrum =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        newSpectrum[row] = (float *) calloc (n, sizeof (float));
    }

    centeredSpectrum =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        centeredSpectrum[row] = (float *) calloc (n, sizeof (float));
    }
    
    // holds real data to output
    realValuesOnly =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        realValuesOnly[row] = (float *) calloc (n, sizeof (float));
    }
    butter =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        butter[row] = (float *) calloc (n, sizeof (float));
    }
    
    newButter =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        newButter[row] = (float *) calloc (n, sizeof (float));
    }
    butteredSpectrum =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        butteredSpectrum[row] = (float *) calloc (2*n, sizeof (float));
    }
    // put input image pixels into origData array

  input = fopen(argv[1], "rb");               // open input image
    for (row = 0; row < n; row++) {
        for (column = 0; column < n; column++) {
            fread(&pixel, sizeof(char), 1, input);
            data[row][column]= (float)pixel;
        }
    }
    
  fclose(input);                      // close input image after all data has been copied
  
    // store data into spectrum because you can't use data without losing the info
    for (row = 0; row < n; row++) {
        for (column = 0; column < n; column++) {
            spectrum[row][column*2] = data[row][column];
        }
    }
    
    four2d(spectrum);
    
    // find the magnitude of the fourier transform
    findMagnitude(spectrum, newSpectrum);

    maximum = getMaximumOfArray(newSpectrum);
    minimum = getMinimumOfArray(newSpectrum);
    
    // normalize the data based on the max and min values with the new max being 255
    normalize(newSpectrum, newSpectrum, minimum, maximum);
    
    for (row=0; row<n/2; row++) {
        for (column=0; column<n/2; column++) {
            float temp = newSpectrum[row][column];
            centeredSpectrum[row][column] = newSpectrum[row+n/2][column+n/2];
            centeredSpectrum[row+n/2][column+n/2] = temp;
            
            temp = newSpectrum[row + n/2][column];
            centeredSpectrum[row + n/2][column] = newSpectrum[row][column+n/2];
            centeredSpectrum[row][column+n/2] = temp;
        }
    }
    
    sprintf(outputFileName, "part_1_spectrum");
    outputFile(outputFileName, centeredSpectrum);
    
    // make butterworth filter
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            float center = sqrt(((column-128)*(column-128)) + ((row-128)*(row-128)));
            butter[row][column] = 1/(1+(sqrt(2)-1)*pow((center/D0),smoothness));
        }
    }
    
    maximum = getMaximumOfArray(butter);
    minimum = getMinimumOfArray(butter);
    
    normalize(butter, newButter, minimum, maximum);
    
    sprintf(outputFileName, "part_2_butterworth_filter");
    outputFile(outputFileName, newButter);
    
    // apply butterworth
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            butteredSpectrum[row][column*2] = spectrum[row][column*2] * newButter[row][column];
            butteredSpectrum[row][column*2+1] = spectrum[row][(column*2)+1] * newButter[row][column];
        }
    }
    
    four2dInverse(spectrum);
    four2dInverse(butteredSpectrum);
    
    // ignore all the imaginary numbers and only grab the real values
    for (row=0; row<n; row++) {
        for (column=0; column<2*n; column++) {
            if (column%2 == 0) {
                realValuesOnly[row][column/2] = spectrum[row][column];
            }
        }
    }
    
    // get the max and min to be able to normalize
    maximum = getMaximumOfArray(realValuesOnly);
    minimum = getMinimumOfArray(realValuesOnly);
    
    normalize(realValuesOnly, realValuesOnly, minimum, maximum);
    
    sprintf(outputFileName, "part_1_inverted_fourier");
    outputFile(outputFileName, realValuesOnly);
    
    for (row=0; row<n; row++) {
        for (column=0; column<2*n; column++) {
            if (column%2 == 0) {
                realValuesOnly[row][column/2] = butteredSpectrum[row][column];
            }
        }
    }

    maximum = getMaximumOfArray(realValuesOnly);
    minimum = getMinimumOfArray(realValuesOnly);
    
    normalize(realValuesOnly, realValuesOnly, minimum, maximum);

    sprintf(outputFileName, "part_2_filtered_image");
    outputFile(outputFileName, realValuesOnly);
    
    
    // free up all dynamically allocated memory
    for (row=0; row< n; row++) {
        free(data[row]);
        free(spectrum[row]);
        free(newSpectrum[row]);
        free(realValuesOnly[row]);
        free(butter[row]);
        free(newButter[row]);
        free(butteredSpectrum[row]);
        free(centeredSpectrum[row]);
    }
    free(data);
    free(spectrum);
    free(newSpectrum);
    free(realValuesOnly);
    free(butter);
    free(newButter);
    free(butteredSpectrum);
    free(centeredSpectrum);

}

float getMaximumOfArray(float **inputtedArray){
    int row, column;
    int n = 256;
    float max;
    max = inputtedArray[0][0];
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            if (inputtedArray[row][column] > max) {
                max = inputtedArray[row][column];
            }
        }
    }
    return max;
}
float getMinimumOfArray(float **inputtedArray){
    int row, column;
    int n = 256;
    float min;
    min = inputtedArray[0][0];
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            if (inputtedArray[row][column] < min) {
                min = inputtedArray[row][column];
            }
        }
    }
    return min;
}
void findMagnitude(float **inputArray, float **outputArray){
    int n = 256;
    int row, column;
    for ( row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            outputArray[row][column] = sqrt(pow(inputArray[row][2*column], 2) + pow(inputArray[row][(2*column)+1], 2));
        }
    }
}

void normalize(float **inputArray, float **outputArray, float min, float max){
    int n = 256;
    int row, column;
    for ( row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            outputArray[row][column] = 255 * (inputArray[row][column] - min)/(max - min);
            outputArray[row][column] = ceil(outputArray[row][column]);
        }
    }
}

void four2d(float **inArray){
    int row, column;
    int n = 256;
    float *temp;
    
    temp =(float*) calloc (2*n, sizeof (float ));
    // do one fourier transform on the rows of data
    for (row = 0; row < n; row ++) {
        for (column = 0; column <2*n; column++) {
            temp[column] = inArray[row][column];
        }
        four1(temp-1, n, 1);
        for (column=0; column<2*n; column++) {
            inArray[row][column] = temp[column];
        }
    }
    
    // now column by column
    for (column = 0; column<n; column++) {
        for (row=0; row<n; row++) {
            temp[row*2] = inArray[row][column*2];
            temp[(row*2)+1] = inArray[row][(column*2)+1];
        }
        four1(temp-1, n, 1);
        for (row=0; row<n; row++) {
            inArray[row][column*2] = temp[row*2];
            inArray[row][(column*2)+1] = temp[(row*2)+1];
        }
    }
    free(temp);

}

void four2dInverse(float **inArray){
    int row, column;
    int n = 256;
    float *temp;
    
    temp =(float*) calloc (2*n, sizeof (float ));
    // time for inverse fourier transform
    //start with columns, since you started fourier with rows
    for (column = 0; column<n; column++) {
        for (row=0; row<n; row++) {
            temp[row*2] = inArray[row][column*2];
            temp[(row*2)+1] = inArray[row][(column*2)+1];
        }
        four1(temp-1, n, -1);
        for (row=0; row<n; row++) {
            inArray[row][column*2] = temp[row*2];
            inArray[row][(column*2)+1] = temp[(row*2)+1];
        }
    }
    
    // next do the rows because that's what you did first
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            temp[(column*2)] = inArray[row][(column*2)];
            temp[(column*2)+1] = inArray[row][(column*2)+1];
        }
        four1(temp-1, n, -1);
        for (column=0; column<2*n; column++) {
            inArray[row][column] = temp[column];
        }
    }
    free(temp);

}

void outputFile(char nameOfFile[256], float **dataArray){
    unsigned char pixel;                // unsigned character to represent one pixel value
    FILE *outputFile;                   // output image file pointer (inverted fourier)
    int column, row;
    int n = 256;
    outputFile = fopen(nameOfFile, "wb");
    for (row=0; row<n; row++) {
        for (column=0; column<n; column++) {
            pixel = (unsigned char)(dataArray[row][column]);
            fwrite(&pixel, sizeof(char), 1, outputFile);
        }
    }
    fclose(outputFile);
}



