#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <random>
#include "Vector.hpp"

//#include <boost/math/constants/constants.hpp>
#include <cmath>

// temp run parms
//#include "atdparms.hpp"


//TODO: This typedef is unnecessary. The code in main.cpp should be updated to
//use Vector straight. 
typedef Vector vec_T;
const float_T pi= 3.14159265358979323846;
//const float_T pi = M_PI;


//This will define our centrosome options. 
enum MTOC {M_CENTROSOME, D_CENTROSOME};

//Debug parameters. 
bool spitValues = false;

//Type Parameters
const bool motherSpringOn = false;
const bool daughterSpringOn = false;
const bool translation = true;
const bool ONLY_COMMA = true;

//Parameters
const float_T Duration = 60;         //duration in minutes
const float_T Tau      = 1.0/5000.0; //time step in minutes
const int max_rot_count= 50;
const float_T max_rot_mag = 0.995;

// MT Parameters
#define MT_numb_M 1000
#define MT_numb_D 1000

vec_T MT_Pos_M[MT_numb_M];
vec_T MT_Pos_D[MT_numb_D];
bool MT_Growing_M[MT_numb_M];
bool MT_Growing_D[MT_numb_D];
float_T MT_GrowthVel_M[MT_numb_M];
float_T MT_GrowthVel_D[MT_numb_D];



//  Contact Parameters:
float_T contact_length = 400*Tau; //This quantity is measured in min

float_T contact_length_dynein = 400*Tau;// This quantity measures the contact time for each MT pullig due to dynein 1/koff
float_T contact_length_push = 30*Tau; // This quantity measures the contact time for each MT pushing on the cortex =1/fcat.

float_T MT_Contact_M[MT_numb_M];
float_T MT_Contact_D[MT_numb_D];

float_T MT_Contact_M_push[MT_numb_M];
float_T MT_Contact_D_push[MT_numb_D];

//  MT Growth/Shrinking Parameters from the DifferentSpringsAndRateTables.pdf
const float_T Vg   = 1*60; //growth velocity in mum/min
//const float_T Vg   = __ATD_PARM_1__;
const float_T Vs_c = 0.5*60; //shortening velocity after contact in mum/min
//const float_T Vs_c = __ATD_PARM_2__;
const float_T Vs   = 0.5*60; //shortening velocity in mum/min
//const float_T Vs   = __ATD_PARM_3__;
const float_T kc   = .04*60;   //catastrophe frequency in min-1
const float_T kr   = .1*60;   //rescue frequency in min-1


float_T Pr_catastrophe = 1-exp(-kc*Tau);
float_T Pr_rescue      = 1-exp(-kr*Tau);


//   Force Parameters
const float_T F_MT_pull   = 1.0; //Force per MT in pN.
//const float_T F_MT_pull   = __ATD_PARM_1__;
const float_T F_MT_push   = 0; //Force per MT in pN.
//const float_T F_MT_push   = __ATD_PARM_2__;

const float_T Fratio = 1.0;// Env. killing: 0.32299363312512985;
vec_T force_M;
vec_T force_D;
vec_T force;
float_T torque_M;
float_T torque_D;
float_T torque;


//   Envelope Parameters
const float_T envWidthM = 2*pi/3;
const float_T envelopeM[2] = {pi/2.0 - envWidthM/2.0, pi/2.0 + envWidthM/2.0}; // The envelope in which MTs from M can grow. 
const float_T envWidthD = 2*pi/3;
const float_T envelopeD[2] = {pi/2.0 - envWidthD/2.0, pi/2.0 + envWidthD/2.0}; // The envelope in which MTs from M can grow.

//   Spring Parameters
const float_T kM = 3;
const float_T kD = 3;

// Pronucleus Parameters
const float_T R1_max = 25;       //Embryo width in mum
const float_T R2_max = 15;       //Embryo length in mum
const float_T Prad   = 4;        //Pronucleous radius (mum)
const float_T Eta    = 1.25;     //Cytoplasmic viscosity =1/60 (pN min/min), Drag coeff = 6 pi Prad/60
const float_T Mu     = 5*Eta;      //Rotational drag coeff ((pN mum)/min)
const float_T Eta2   = Eta/10;   //Trans. drag coeff of each MTOC (pN min/mum)
const float_T kbT    = .00414;   //pN mum
const float_T D      = kbT/Eta2; //Diffusion Coefficient for each MTOC free motion. (mum^2/min)

//  Starting Centered Coordinates:
//const float_T startPsi = pi/2.0;
//const float_T startX   = 0;
//const float_T startY   = 0;

//   Off-center Coordinates
const float_T startPsi = pi/2.0;
const float_T startX   = 10.0;
const float_T startY   = 0;

//  General Coordinate Initializations: 
float_T psi;
vec_T basePosM;
vec_T basePosD;
vec_T proNucPos;

// Band parameters: 
//  regionAngles:
//   The regionAngles vector defines the endpoints of the angular span of
//   various regions along the cortex. It MUST start with 0 and it MUST end with
//   2*pi. Each region is presumed to be homogenous and all encompassing: In
//   particular, a single band of a protein on the cortex does not totally
//   define a region unless it is the only protein defining connectivity
//   probability in that region. If two or more bands intersect, you must make
//   your regions various homogenous intersections such that the probability is
//   constant throughout any given region. Note that technically speaking, in
//   implementation a region is defined by two endpoints \theta_1, \theta_2 \in
//   [0,2\pi] such that a point on the cortex at angle \alpha is in that region
//   if \theta_1 <= \alpha < \theta_2
//  regionProbabilities: 
//    This defines the probabilities associated with the regions defined by the
//    regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 
//  regionForceMultipliers:
//    This defines the force multipliers associated with the regions defined by
//    the regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 

//No Bands:
const int numRegions = 1;
const float_T regionAngles[numRegions + 1] = {0, 2*pi};
const float_T regionProbabilities[numRegions] = {1};
const float_T regionForceMultipliers[numRegions] = {1};

//Standard Bands:
//const int numRegions = 5;
//const float_T start = 0.45*pi;
//const float_T end = 0.5*pi;
//const float_T regionAngles[numRegions+1] = {0, start, end, 2*pi - end, 2*pi - start, 2*pi};
//const float_T regionProbabilities[numRegions] = {1,1,1,1,1};
//const float_T regionForceMultipliers[numRegions] = {1,1,1,1,1};

// MT density limitations
//  Only one contact per window, windows of length ~1
const size_t numberContactWindows = 128;
const float_T contactWindowAngles[numberContactWindows+1] =
//1.Initial cortical angles
/*{0., 0.039985, 0.079885, 0.119585, 0.159085, 0.198485, 0.237785, 0.277085,
 0.316485, 0.356085, 0.396085, 0.436485, 0.477585, 0.519485, 0.562285,
 0.606085, 0.651085, 0.697385, 0.745085, 0.794285, 0.845085, 0.897585,
 0.951785, 1.00768, 1.06538, 1.12468, 1.18558, 1.24798, 1.31168, 1.37648,
 1.44218, 1.50848, 1.57508, 1.64168, 1.70788, 1.77348, 1.83818, 1.90178,
 1.96398, 2.02468, 2.08378, 2.14119, 2.19689, 2.25089, 2.30309, 2.35369,
 2.40269, 2.45019, 2.49629, 2.54109, 2.58479, 2.62739, 2.66909, 2.71009,
 2.75049, 2.79039, 2.82999, 2.86939, 2.90869, 2.94799, 2.98739, 3.02699,
 3.06669, 3.10659, 3.14649, 3.18639, 3.22619, 3.26589, 3.30539, 3.34479,
 3.38409, 3.42339, 3.46279, 3.50239, 3.54239, 3.58289, 3.62409, 3.66609,
 3.70899, 3.75289, 3.79799, 3.84439, 3.89219, 3.94159, 3.99259, 4.04529,
 4.09969, 4.15579, 4.21359, 4.27309, 4.33419, 4.39679, 4.46069, 4.52569,
 4.59149, 4.65779, 4.72439, 4.79089, 4.85709, 4.92259, 4.98719, 5.05059,
 5.11269, 5.17329, 5.23219, 5.28939, 5.34489, 5.39869, 5.45069, 5.50109,
 5.54989, 5.59719, 5.64309, 5.68779, 5.73128, 5.77378, 5.81548, 5.85638,
 5.89668, 5.93658, 5.97608, 6.01538, 6.05468, 6.09398, 6.13338, 6.17298,
 6.21268, 6.25258, 6.28319};*/
//2. Corrected angles
/*{0., 0.0398, 0.0795, 0.1191, 0.1585, 0.1978, 0.237, 0.2762, 0.3155,
 0.355, 0.3948, 0.4351, 0.476, 0.5177, 0.5603, 0.6039, 0.6487, 0.6948,
 0.7423, 0.7913, 0.8418, 0.894, 0.9479, 1.0035, 1.0609, 1.1199,
 1.1805, 1.2426, 1.306, 1.3706, 1.436, 1.5021, 1.5685, 1.6349, 1.701,
 1.7665, 1.8311, 1.8946, 1.9568, 2.0176, 2.0768, 2.1343, 2.19, 2.244,
 2.2963, 2.347, 2.3961, 2.4437, 2.4899, 2.5348, 2.5785, 2.6212, 2.663,
 2.704, 2.7443, 2.7842, 2.8237, 2.863, 2.9022, 2.9414, 2.9807, 3.0201,
 3.0597, 3.0994, 3.1392, 3.179, 3.2187, 3.2583, 3.2977, 3.337, 3.3762,
 3.4154, 3.4546, 3.4941, 3.5339, 3.5742, 3.6151, 3.6568, 3.6993,
 3.7429, 3.7876, 3.8336, 3.881, 3.9299, 3.9803, 4.0324, 4.0862,
 4.1417, 4.199, 4.2579, 4.3185, 4.3805, 4.4439, 4.5084, 4.5738,
 4.6398, 4.7062, 4.7726, 4.8387, 4.9042, 4.9689, 5.0325, 5.0948,
 5.1556, 5.2148, 5.2724, 5.3283, 5.3824, 5.4348, 5.4855, 5.5347,
 5.5823, 5.6285, 5.6734, 5.7172, 5.7599, 5.8017, 5.8427, 5.8831,
 5.923, 5.9625, 6.0018, 6.041, 6.0802, 6.1195, 6.1589, 6.1985, 6.2382,
 6.278};*/
//3.Cortical angles split in 4 quadrants.
{0.0000, 0.0400, 0.0799, 0.1197, 0.1593, 0.1988, 0.2382, 0.2776, 0.3171, 0.3568,
    0.3968, 0.4373, 0.4785, 0.5205, 0.5634, 0.6073, 0.6524,0.6988, 0.7466,0.7959,
    0.8468,0.8994, 0.9537,1.0097,1.0674,1.1268,1.1878,1.2503,1.3141,1.3790,1.4447,
    1.5110, 1.5708, 1.6306,1.6969,1.7626,1.8275,1.8913,1.9538,2.0148,2.0742,2.1319,
    2.1879,2.2422,2.2948,2.3457,2.3950,2.4428,2.4892,2.5343,2.5782,2.6211,2.6631,2.7043,
    2.7448,2.7848,2.8245,2.8640,2.9034,2.9428,2.9823,3.0219,3.0617,3.1016,3.1416,3.1816,
    3.2215,3.2613,3.3009,3.3404,3.3798,3.4192,3.4587,3.4984,3.5384,3.5789,3.6201,3.6621,
    3.7050,3.7489,3.7940,3.8404,3.8882,3.9375,3.9884,4.0410,4.0953,4.1513,4.2090,4.2684,
    4.3294,4.3919,4.4557,4.5206,4.5863,4.6526,4.7124,4.7722,4.8385,4.9042,4.9691,5.0329,
    5.0954,5.1564,5.2158,5.2735,5.3295,5.3838,5.4364,5.4873,5.5366,5.5844,5.6308,5.6759,
    5.7198,5.7627,5.8047,5.8459,5.8864,5.9264,5.9661,6.0056,6.0450,6.0844,6.1239,6.1635,
    6.2033,6.2432,6.2832};
bool contacts[numberContactWindows];

// File Parameters
std::ofstream file;
std::string fileName  = "";
std::string fileDir   = "../data/";
std::string fileOrder = "t,proNucPos,psi,MT_Pos_M,MT_Pos_D,force_M,force_D,\
                         force,torque_M,torque_D,torque,basePosM,basePosD";

// Random Number Generation Parameters: 
std::normal_distribution<float_T> stdNormalDist;
std::uniform_real_distribution<float_T> stdUniformDist(0.0,1.0);
std::default_random_engine generator;
