#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[]){
    int row, column, n;
    float **data;
    FILE *output1;
    FILE *input;
    n=256;
    
    //input image file
    input = fopen(argv[1], "rb");
    if (input == NULL){
        return 1;
    }
    
    output1 = fopen(argv[2], "wb");
    if (output1 == NULL){
        return 1;
    }
    

    data =(float**) calloc (n, sizeof (float *));
    for (row = 0; row < n; row++) {
        data[row] = (float *) calloc (n, sizeof (float));
    }
    
    //same as static 2-d array
    for (row = 0; row < n; row ++) {
        for (column = 0; column < n; column++) {
            fread(&data[row][column], sizeof(char), 1, input);

        }
    }
    for (row = 0; row < n; row ++) {
        for (column = 0; column < n; column++) {
            fwrite(&data[row][column], sizeof(char), 1, output1);
        }
    }
    

    for (row=0; row< n; row++) {
        free(data[row]);
    }
    free(data);
    
    fclose(output1);
    fclose(input);
    
    return 0;
}

