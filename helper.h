//NORMALIZES A VECTOR
void NORMALIZE(float *rayDirection){
    float norm = sqrt(rayDirection[0]*rayDirection[0] + rayDirection[1]*rayDirection[1] +rayDirection[2]*rayDirection[2] );
    rayDirection[0]=rayDirection[0]/norm;
    rayDirection[1]=rayDirection[1]/norm;
    rayDirection[2]=rayDirection[2]/norm;
}

//FINDS THE PIXEL INDEX IN A SINGLE ARRAY BASED ON ROW AND COLUM 
int PIXEL_INDEX(int row, int column,int color){
    int index = row * 3 * 512 ;
    index  += column*3;
    index  += color;
    return index;
}

//FINDS THE DOT PRODUCT OF TWO VECTORS
float DOT_PRODUCT(float vec1[3],float vec2[3]){
    float output=0;
    output += (vec1[0]*vec2[0]);
    output += (vec1[1]*vec2[1]);
    output += (vec1[2]*vec2[2]);
    return output;
}

//ADDS TWO VECTORS 
void ADD_VEC3(float vec1[3],float vec2[3],float output[3]){
    output[0]=vec1[0]+vec2[0];
    output[1]=vec1[1]+vec2[1];
    output[2]=vec1[2]+vec2[2];
}

//SUBTRACTS TWO VECTORS
void SUB_VEC3(float vec1[3],float vec2[3],float output[3]){
    output[0]=vec1[0]-vec2[0];
    output[1]=vec1[1]-vec2[1];
    output[2]=vec1[2]-vec2[2];
}

//MULT A VECTOR BY A SINGLE VALUE
void SCALER_MULT(float scaler,float vec1[3],float output[3]){
    output[0] = vec1[0]*scaler;
    output[1] = vec1[1]*scaler;
    output[2] = vec1[2]*scaler;
}
//CALCULATES THE CROSS PRODUCT OF TWO VECTORS
void CROSS_PRODUCT(float vec1[3],float vec2[3],float output[3]){
    output[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];  
    output[1] = vec1[2]*vec2[1] - vec1[0]*vec2[2];
    output[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}
