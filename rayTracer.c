#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <math.h>
#include <stdio.h>
#include "helper.h"

#define RED (0) 
#define GREEN (1)
#define BLUE (2)
#define WIDTH (512)
#define HEIGHT (512)
int numberOfTriangles;
int numberOfSpheres;

/* RAY STRUCT */
typedef struct{  
    float start_Position[3];    //STARTING POSITION OF THE RAY
    float direction_Vector[3];  //DIRECTION OF THE RAY
} Ray;

/* PERSPECTIVE STRUCT */
typedef struct{
  float camera_Position[3];     //CAMERA POSITION
  float distance_To_Screen;     //DISTANCE TO SCREEN
  int Screen_Width_World;       //SCREEN WIDTH IN WORLD COORDINATES
  int Screen_Width_Pixels;      //SCREEN WITHD IN PIXELS
  float Light_Position[3];      //LIGHT POSITION IN WORLD COORDINATES
} Perspective;

/*BUILDS PERSPECTIVE STRUCT*/
Perspective persp = {.camera_Position ={0,0,0}, .Light_Position={3,5,-15}, .Screen_Width_Pixels=512, .Screen_Width_World=2, .distance_To_Screen = -2,};

/* TRIANGLE STRUCT */
typedef struct{
    float point_A[3];           //POINT 1 OF THE TRIANGLE
    float point_B[3];           //POINT 2 OF THE TRIANGLE
    float point_C[3];           //POINT 3 OF THE TRIANGLE
    float Color[3];             //COLOR OF TRIANGLE
    int reflective;             //IF TRIANGLE IS REFLECTIVE
} Triangle;
Triangle *TRIANGLES;

/* SPHERE STRUCT */
typedef struct{
    float Center[3];            //SPHERE CENTER IN WORLD COORDINATES
    float Radius;               //RADIUS OF THE SPHERE
    float Color[3];             //THE COLOR OF THE SPHERE, NOT USED IF REFLECTIVE
    int Reflective;             //IF THE SPHERE IS REFLECTIVE
} Sphere;
Sphere *SPHERES; 

/* RAYHIT STRUCT */
typedef struct{  
    int MISS;                   //IF THERE WAS A HIT
    float time;                 //TIME TO HIT
    float normal[3];            //NORMAL OF HIT RAY
    int reflective;             //IF THE HIT HIT A REFLECTIVE OBJ
    float color[3];             //COLOR OF OBJ HIT
    float incomingRay[3];       //INCOMING RAY OF HIT OBJ
    float hitLocation[3];       //HIT LOCATION OF HIT OBJ              
} RayHit;

/*BUILDS A RAY STRUCT BASED ON A POSITION AND DIRECTION*/
void setRayDirection(float *ScreenPosition, Ray *RAY){
    RAY->direction_Vector[0]=ScreenPosition[0];
    RAY->direction_Vector[1]=ScreenPosition[1];
    RAY->direction_Vector[2]=-2;
}

/*SEES IF THE RAY INTERSECNT WITH A SPHERE RETURN 0 FOR NO 1 FOR YES */
RayHit sphereIntersect(Sphere *sphere ,Ray *ray){
    RayHit rayhit;

    /*CALCULATES IF THERE WAS A HIT FIRST FINDS T */
    float zero_vec[]={0,0,0};
    float neg_d[3]; 
    SUB_VEC3(zero_vec,ray->direction_Vector,neg_d);
    float e_min_c[3];
    SUB_VEC3(ray->start_Position,sphere->Center,e_min_c); 
    float e_min_c_dot_e_min_c = DOT_PRODUCT(e_min_c,e_min_c);
    float d_dot_d = DOT_PRODUCT(ray->direction_Vector,ray->direction_Vector);
    float neg_d_dot_e_min_c = DOT_PRODUCT(neg_d,e_min_c);
    float discriminant = (neg_d_dot_e_min_c*neg_d_dot_e_min_c) - d_dot_d * (e_min_c_dot_e_min_c - (sphere->Radius*sphere->Radius));
    
    /*NO HIT THERE WAS A MISS*/
    if(discriminant < 0){    
        rayhit.MISS=1;
        return rayhit; 
    }
  
    /*CALCULATES THE SMALLEST T VALUE THATS POSITIVE*/
    float t1 = ((neg_d_dot_e_min_c + sqrt(discriminant))/d_dot_d);
    float t2 = ((neg_d_dot_e_min_c - sqrt(discriminant))/d_dot_d);

    if(t2 > 0)
        rayhit.time=t2;
    else 
         rayhit.time=t1;

    /*FINDS THE POSITION OF THE HIT*/
    float HitPos[3];
    float scaler[3];
    SCALER_MULT(rayhit.time,ray->direction_Vector,scaler);
    ADD_VEC3(ray->start_Position,scaler,HitPos);
    rayhit.hitLocation[0]=HitPos[0];
    rayhit.hitLocation[1]=HitPos[1];
    rayhit.hitLocation[2]=HitPos[2];

    //SETS THE NORMAL
    float tmp[3];
    SUB_VEC3(rayhit.hitLocation,sphere->Center,tmp);
    NORMALIZE(tmp);
    rayhit.normal[0] = tmp[0];
    rayhit.normal[1] = tmp[1];
    rayhit.normal[2] = tmp[2];

    /*SETS THE REST OF THE STRUCT*/
    rayhit.MISS=0;
    rayhit.color[0] = sphere->Color[0];
    rayhit.color[1] = sphere->Color[1];
    rayhit.color[2] = sphere->Color[2];
    rayhit.reflective = sphere->Reflective;
    rayhit.incomingRay[0]= ray->direction_Vector[0];
    rayhit.incomingRay[1]= ray->direction_Vector[1];
    rayhit.incomingRay[2]= ray->direction_Vector[2];

    return rayhit;
}

/*FINDS IF THERE WAS A HIT WITH A TRIANGLE*/
RayHit triangleIntersect(Triangle *triangle, Ray *ray){
    /*MATH TO FIND T VALUE AND IF A HIT OCCURED*/
    RayHit rayhit;
    rayhit.MISS=0;
    float T,A,B,C,D,E,F,G,H,I,J,K,L;
    int X=0,Y=1,Z=2;

    A = triangle->point_A[X] -  triangle->point_B[X];
    B = triangle->point_A[Y] -  triangle->point_B[Y];
    C = triangle->point_A[Z] -  triangle->point_B[Z];

    D = triangle->point_A[X] -  triangle->point_C[X];
    E = triangle->point_A[Y] -  triangle->point_C[Y];
    F = triangle->point_A[Z] -  triangle->point_C[Z];

    G = ray->direction_Vector[X];
    H = ray->direction_Vector[Y];
    I = ray->direction_Vector[Z];

    J = triangle->point_A[X] - ray->start_Position[X];
    K = triangle->point_A[Y] - ray->start_Position[Y];
    L = triangle->point_A[Z] - ray->start_Position[Z];

    float M = A*(E*I - H*F) + B*(G*F - D*I) + C*(D*H - E*G);

    T = -(F*(A*K - J*B) + E*(J*C - A*L) + D*(B*L - K*C)) / M;

    /*CHECKS T FOR NO HIT*/
    if(T < 0){  
        rayhit.MISS = 1;
        return rayhit;
    }
    /*CHECKS GAMMA FOR NOT HIT */
    float gamma = (I*(A*K - J*B) + H*(J*C - A*L) + G*(B*L - K*C)) / M;
    if(gamma < 0 || gamma > 1){
        rayhit.MISS = 1;
        return rayhit;
    }
    /*CHECKS BETA FOR NO HIT*/
    float beta = (J*(E*I - H*F) + K*(G*F - D*I) + L*(D*H - E*G)) / M;
    if(beta < 0 || beta > 1-gamma){
        rayhit.MISS = 1;
        return rayhit;
    }

    /*FILLS RAYHIT STRUCT*/
    rayhit.time = T;
    rayhit.color[0] = triangle->Color[0];
    rayhit.color[1] = triangle->Color[1];
    rayhit.color[2] = triangle->Color[2];
    rayhit.reflective = triangle->reflective;
    
    /*FINDS THE INTERSECTION POINT*/
    float HitPos[3];
    float scaler[3];
    SCALER_MULT(rayhit.time,ray->direction_Vector,scaler);
    ADD_VEC3(ray->start_Position,scaler,HitPos);
    rayhit.hitLocation[0]=HitPos[0];
    rayhit.hitLocation[1]=HitPos[1];
    rayhit.hitLocation[2]=HitPos[2];

    /*FINDS THE NORMAL OF THE HIT POINT*/
    float vec1[3];
    SUB_VEC3(triangle->point_A,triangle->point_B,vec1);
    float vec2[3];
    SUB_VEC3(triangle->point_C,triangle->point_B,vec2);
    float normal[3];
    CROSS_PRODUCT(vec2,vec1,normal);
    NORMALIZE(normal);
    rayhit.normal[0]= normal[0];
    rayhit.normal[1]= normal[1];
    rayhit.normal[2]= normal[2];
    rayhit.incomingRay[0]= ray->direction_Vector[0];
    rayhit.incomingRay[1]= ray->direction_Vector[1];
    rayhit.incomingRay[2]= ray->direction_Vector[2];

    return rayhit;
}

/*BUILDS THE SPHERE STRUCT WITH PASSED IN INFORMATION*/
 void FillSphereStruct(Sphere *sphere ,float centerX,float centerY,float centerZ,float radius,float colorR,float colorG,float colorB,int refl){
     sphere->Center[0]=centerX;
     sphere->Center[1]=centerY;
     sphere->Center[2]=centerZ;
     sphere->Color[0]=colorR;
     sphere->Color[1]=colorG;
     sphere->Color[2]=colorB;
     sphere->Radius=radius;
     sphere->Reflective=refl;
 }

 /*BUILDS THE TRIANGLE STRUCT WITH PASSED IN INFORMATION*/
 void FillTrianlgeStruct(Triangle *triangle ,float pointA_X,float pointA_Y,float pointA_Z,float pointB_X,float pointB_Y,float pointB_Z,float pointC_X,float pointC_Y,float pointC_Z,float R,float G,float B,int reflective){
     triangle->point_A[0] = pointA_X;
     triangle->point_A[1] = pointA_Y;
     triangle->point_A[2] = pointA_Z;

     triangle->point_B[0] = pointB_X;
     triangle->point_B[1] = pointB_Y;
     triangle->point_B[2] = pointB_Z;

     triangle->point_C[0] = pointC_X;
     triangle->point_C[1] = pointC_Y;
     triangle->point_C[2] = pointC_Z;

     triangle->reflective = reflective;
     triangle->Color[0] = R;
     triangle->Color[1] = G;
     triangle->Color[2] = B;
 }

/*GETS THE RAY FROM THE SCREEN COORDINATE*/
void getRay(int x, int y,Ray *ray){
    float CamPos[3];
    CamPos[0] = persp.camera_Position[0];
    CamPos[1] = persp.camera_Position[1];
    CamPos[2] = persp.camera_Position[2];

    float PixelPos[3];
    PixelPos[0] = -.998047 + (x*.003906);
    PixelPos[1] = .998047 - (y*.003906);
    PixelPos[2] = persp.distance_To_Screen;

    float output[3];
    SUB_VEC3(PixelPos,CamPos,output);
    NORMALIZE(output);

    ray->start_Position[0] = persp.camera_Position[0];
    ray->start_Position[1] = persp.camera_Position[1];
    ray->start_Position[2] = persp.camera_Position[2];
    ray->direction_Vector[0]= output[0];
    ray->direction_Vector[1]= output[1];
    ray->direction_Vector[2]= output[2];
}

/*GETS THE COLOR OF THE CURRENT RAYHIT*/
 void getColor(RayHit rayhit,float *NewColor,int reflectCounter){
    int z,SHADOW=0;
    float toLight[3];
    SUB_VEC3(persp.Light_Position,rayhit.hitLocation, toLight);
    NORMALIZE(toLight);

    /*MAX NUMBER OF REFLECTIVE OBJECTS HIT SETS TO BLACK*/
    if(rayhit.reflective){
        if(reflectCounter > 10){
            NewColor[0] = 0;
            NewColor[1] = 0;
            NewColor[2] = 0;
            return;
        }
        /*GETS REFLECTING OUT RAY*/
        float reflectRay[3];
        float dot = DOT_PRODUCT(rayhit.incomingRay,rayhit.normal);
        dot = 2 * dot;
        float tmp[3];
        SCALER_MULT(dot,rayhit.normal,tmp);
        SUB_VEC3(rayhit.incomingRay,tmp,reflectRay);
        NORMALIZE(reflectRay);
        Ray newRay = { .start_Position = {rayhit.hitLocation[0]+reflectRay[0]*.01,rayhit.hitLocation[1]+reflectRay[1]*.01,rayhit.hitLocation[2]+reflectRay[2]*.01},.direction_Vector={reflectRay[0],reflectRay[1],reflectRay[2]}};
        /*RUNS THROUGH ALL OF THE TRIANGLES TO FIND CLOSEST HIT*/
        RayHit nextRayhit;
        nextRayhit.MISS=1;
        nextRayhit.time=10000000;
        /*RUNS THROUGH ALL OF THE SPHERES TO FIND CLOSEST HIT*/
        for(z=0;z<numberOfSpheres;z++){
            RayHit tmpRayHit = sphereIntersect(&SPHERES[z] ,&newRay);
            if(!tmpRayHit.MISS && tmpRayHit.time > 0  && tmpRayHit.time < nextRayhit.time)
                nextRayhit = tmpRayHit; 
        }
        /*RUNS THROUGH ALL OF THE TRIANGLES TO FIND CLOSEST HIT*/
        for(z=0;z<numberOfTriangles;z++){
            RayHit tmpRayHit = triangleIntersect(&TRIANGLES[z] ,&newRay);
            if(!tmpRayHit.MISS && tmpRayHit.time > 0  && tmpRayHit.time < nextRayhit.time)
                nextRayhit = tmpRayHit;       
        }
        /*NO HIT SETS COLOR TO BLACK*/
        if(nextRayhit.MISS){
            NewColor[0] = 0;
            NewColor[1] = 0;
            NewColor[2] = 0;
        }else{
            /*RUNS AGAIN FOR REFLECTIVE OBJECTS*/
            getColor(nextRayhit,NewColor,reflectCounter+1);    
        }
        return;
    }else{
        /*CALCULATES DISTANCE TO LIGHT*/
        float Xpart = (rayhit.hitLocation[0]-persp.Light_Position[0]) * (rayhit.hitLocation[0]-persp.Light_Position[0]);
        float Ypart = (rayhit.hitLocation[1]-persp.Light_Position[1]) * (rayhit.hitLocation[1]-persp.Light_Position[1]);
        float Zpart =(rayhit.hitLocation[2]-persp.Light_Position[2]) * (rayhit.hitLocation[2]-persp.Light_Position[2]);
        float distanceToLight = sqrt( Xpart + Ypart + Zpart );
        /*NEW REFLECTIVE RAY*/
        Ray ray = {.start_Position={rayhit.hitLocation[0]+toLight[0]*.0001,rayhit.hitLocation[1]+toLight[1]*.0001,rayhit.hitLocation[2]+toLight[2]*.0001},
                    .direction_Vector={toLight[0],toLight[1],toLight[2]}};
        /*CHECK IF THE PIXEL IS IN VIEW OF LIGHT SOURCE*/
        /*RUNS THROUGH ALL OF THE TRIANGLES TO FIND CLOSEST HIT*/
        for(z=0;z<numberOfTriangles;z++){
            RayHit Hit = triangleIntersect(&TRIANGLES[z] , &ray);   
                if(Hit.MISS == 0 ){
                    if(Hit.time < distanceToLight && Hit.time > 0)
                        SHADOW = 1;
                }
        }
        /*RUNS THROUGH ALL OF THE SPHERES TO FIND CLOSEST HIT*/
        for(z=0;z<numberOfSpheres;z++){
            RayHit Hit = sphereIntersect(&SPHERES[z] , &ray);   
                if(Hit.MISS == 0 ){
                    if(Hit.time < distanceToLight && Hit.time > 0)
                        SHADOW = 1;
                }
        }
        /*CALCULATES SHADOW COLOR*/
        if(SHADOW){
            NewColor[0] = .2 * rayhit.color[0];
            NewColor[1] = .2 * rayhit.color[1];
            NewColor[2] = .2 * rayhit.color[2];
        }else{
            /*ADDS DIFUSE SHADING*/
            float dotProduct = DOT_PRODUCT(rayhit.normal,toLight);
            if(dotProduct < .2)
                dotProduct = .2;
            NewColor[0] = dotProduct * rayhit.color[0];
            NewColor[1] = dotProduct * rayhit.color[1];
            NewColor[2] = dotProduct * rayhit.color[2];
        }
    }
 }


int main (int argc, char** argv){
    if(argc != 2){
        printf("error: must pass argument for file to be created\n");
    }
    int i,x,y,z;
    /*LOADS OBJECTS FOR REFERENCE.PNG*/
    if(0== strcmp(argv[1], "reference")){
        //*BUILDS SPHERES*/
        numberOfSpheres=3;
        SPHERES = (Sphere *)malloc(numberOfSpheres * sizeof(Sphere));
        FillSphereStruct(&SPHERES[0],0,0,-16,2,255,0,255,1);
        FillSphereStruct(&SPHERES[1],3,-1,-14,1,0,0,255,1);
        FillSphereStruct(&SPHERES[2],-3,-1,-14,1,255,0,0,0);
        //*BUILDS TRIANGLES*/ 
        numberOfTriangles=5;
        TRIANGLES = (Triangle *)malloc(numberOfTriangles * sizeof(Triangle));
        //BACK
        FillTrianlgeStruct(&TRIANGLES[0],  -8,-2,-20,  8,-2,-20,  8,10,-20,  0,0,255 , 0);
        FillTrianlgeStruct(&TRIANGLES[1],  -8,-2,-20,  8,10,-20,  -8,10,-20,  0,0,255 , 0);
        //FLOOR
        FillTrianlgeStruct(&TRIANGLES[2],   8,-2,-10,   8,-2,-20, -8,-2,-20,  255,255,255 , 0);
        FillTrianlgeStruct(&TRIANGLES[3],  -8,-2,-20,  -8,-2,-10,  8,-2,-10,  255,255,255 , 0);
        //RIGHT
        FillTrianlgeStruct(&TRIANGLES[4],  8,-2,-20,  8,-2,-10,  8,10,-20,  255,0,0 , 0);
    } else if(0== strcmp(argv[1], "custom")){ /*LOADS OBJECTS FOR CUSTOM.PNG*/
        //SETS LIGHT
        persp.Light_Position[0] = 0;
        persp.Light_Position[1] = 0;
        persp.Light_Position[2] = -15;

        //*BUILDS SPHERES*/
        numberOfSpheres=7;
        SPHERES = (Sphere *)malloc(numberOfSpheres * sizeof(Sphere));
        //SNOWMAN1
        FillSphereStruct(&SPHERES[0],-3,-4,-20,2,200,200,200,1);
        FillSphereStruct(&SPHERES[1],-3,-1,-20,1.5,200,200,200,1);
        FillSphereStruct(&SPHERES[2],-3,1.25,-20,1,200,200,200,1);
        //SNOWMAN2
        FillSphereStruct(&SPHERES[3],2,-4,-20,2,200,200,200,1);
        FillSphereStruct(&SPHERES[4],2,-1,-20,1.5,200,200,200,1);
        FillSphereStruct(&SPHERES[5],2,1.25,-20,1,200,200,200,1);
        FillSphereStruct(&SPHERES[6],3,-5,-14,2,0,156,98,0);

        //*BUILDS TRIANGLES*/ 
        numberOfTriangles=10;
        TRIANGLES = (Triangle *)malloc(numberOfTriangles * sizeof(Triangle));
        //BACK
        FillTrianlgeStruct(&TRIANGLES[0],   -15,-5,-30,   8,-5,-30,  8,10,-30,   255,0,0 , 0);
        FillTrianlgeStruct(&TRIANGLES[1],   -15,-5,-30,   8,10,-30,  -15,10,-30,  255,0,0 , 0);
        //FLOOR
        FillTrianlgeStruct(&TRIANGLES[2],   8,-5,-10,   8,-5,-30,   -15,-5,-30,  0,255,0, 0);
        FillTrianlgeStruct(&TRIANGLES[3],   -15,-5,-30,  -15,-5,-10,  8,-5,-10,   0,255,0 , 0);
        //RIGHT
        FillTrianlgeStruct(&TRIANGLES[4],   8,-5,-30,    8,-5,-10,  8,10,-30,   0,0,255 , 0);
        FillTrianlgeStruct(&TRIANGLES[5],   8,10,-30,    8,-5,-10,  8,10,-10,   0,0,255, 0);
        //LEFT
        FillTrianlgeStruct(&TRIANGLES[8],   -8,-5,-10   ,-8,-5,-30, -8,10,-30,  255,255,0 , 0);
        FillTrianlgeStruct(&TRIANGLES[9],   -8,10,-30,   -8,10,-10, -8,-5,-10,  255,255,0, 0);
        //TOP
        FillTrianlgeStruct(&TRIANGLES[6],    8,10,-30,  8,10,-10,   -30,10,-30,  0,255,255 , 0);
        FillTrianlgeStruct(&TRIANGLES[7],   -30,10,-30,  8,10,-10,   -30,10,-10,  0,255,255 , 0);
    }else{
        printf("Incorrect parameter: must be reference.png or custom.png\n");
        exit(1);
    }

    /*BUILDS AN ARRAY OF PIXELS, SETS OF 3 RGB ALL SET TO BLACK*/
    int NumberOfPixels =(512*512*3);
    unsigned char arrayContainingImage[NumberOfPixels];
    for(i=0 ;i<NumberOfPixels ; i++)
             arrayContainingImage[i] = 0;  
    /*RUNS THROUGH EACH PIXEL*/
    for(y=0; y<512; y++){
         for(x=0; x<512; x++) {
            /*GETS THE RAY GOING THOUGH THE PIXEL*/
            Ray RAY;
            getRay(x,y,&RAY);
            RayHit rayhit= {.MISS=1,.time=1000000000};

            /*RUNS THROUGH ALL OF THE SPHERES TO FIND CLOSEST HIT*/
            for(z=0;z<numberOfSpheres;z++){
                RayHit tmpRayHit = sphereIntersect(&SPHERES[z] ,&RAY);
                if(!tmpRayHit.MISS &&tmpRayHit.time < rayhit.time)
                    rayhit = tmpRayHit; 
            }
            /*RUNS THROUGH ALL OF THE TRIANGLES TO FIND CLOSEST HIT*/
            for(z=0;z<numberOfTriangles;z++){
                RayHit tmpRayHit = triangleIntersect(&TRIANGLES[z] ,&RAY);
                if(!tmpRayHit.MISS && tmpRayHit.time < rayhit.time)
                    rayhit = tmpRayHit; 
            }
           /*SETS THE COLOR OF HIT PIXELS*/
            if(rayhit.MISS == 0){   
                    float NewColor[3]={0,0,0};
                    getColor(rayhit,NewColor,0);
                    int index = PIXEL_INDEX(y,x,0);
                    arrayContainingImage[index]=NewColor[0];
                    arrayContainingImage[index+1]=NewColor[1];
                    arrayContainingImage[index+2]=NewColor[2];
            }   
         }
    }

    /*DRAWS THE FILE OUT*/
    if(0== strcmp(argv[1], "reference"))
        stbi_write_png("reference.png", 512, 512, 3, arrayContainingImage, 512*3);
    if(0== strcmp(argv[1], "custom"))
        stbi_write_png("custom.png", 512, 512, 3, arrayContainingImage, 512*3);

    /*FREES ALL MALLOC SPACE*/
    free(SPHERES);
    free(TRIANGLES);
    printf("PROGRAM COMPLETED: check for created file\n");
    return 0;
}

