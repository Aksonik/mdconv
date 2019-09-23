#include <stdio.h>
#include <stdlib.h>

char main(){

    FILE *myFile;
    myFile = fopen("test.ndx", "r");

    //read file into array
    char numberArray[7];
    int i;

    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);
    }

    for (i = 0; i < 7; i++){
        fscanf(myFile, "%s,", &numberArray[i] );
    }

    for (i = 0; i < 7; i++){
        printf("Number is: %s\n", numberArray[i]);
    }

    fclose(myFile);

    return 0;
}
