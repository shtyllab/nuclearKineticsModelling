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
const float_T Duration = 160;         //duration in minutes
const float_T Tau      = 1.0/4000.0; //time step in minutes

// MT Parameters
// TODO THIS IS REALLY BAD!!! THIS CODE NEEDS TO BE IMPROVED. ... but it should
// work. In the meantime, until it gets really fixed, if you need to change
// this, be careful you know what you're doing with pre-processor macros.
#define MT_numb_M 1000
#define MT_numb_D 1000
// TODO We can't currently use the uncommented lines as there is no easy way to
// dynamically define functions for printing things nicely if these numbers are
// different, unless we use pre-processor macros. This is not acceptable, we
// need to just use object encapsulation. 
//const unsigned MT_numb_M = 100; //# of MTs from the M centrosome.
//const unsigned MT_numb_D = 100; //# of MTs from the D centrosome.
vec_T MT_Pos_M[MT_numb_M];
vec_T MT_Pos_D[MT_numb_D];
bool MT_Growing_M[MT_numb_M];
bool MT_Growing_D[MT_numb_D];
float_T MT_GrowthVel_M[MT_numb_M];
float_T MT_GrowthVel_D[MT_numb_D];


//  Contact Parameters:
float_T contact_length = 500*Tau; //This quantity is measured in min = 0.3750 min
float_T MT_Contact_M[MT_numb_M];
float_T MT_Contact_D[MT_numb_D];


//  MT Growth/Shrinking Parameters from the DifferentSpringsAndRateTables.pdf
const float_T Vg   = 0.3*60; //growth velocity in mum/min
const float_T Vs_c = 0.5*60; //shortening velocity after contact in mum/min
const float_T Vs   = 0.5*60; //shortening velocity in mum/min
const float_T kc   = .04*60;   //catastrophe frequency in min-1
const float_T kr   = .1*60;   //rescue frequency in min-1


float_T Pr_catastrophe = 1-exp(-kc*Tau);
float_T Pr_rescue      = 1-exp(-kr*Tau);


//   Force Parameters
const float_T F_MT   = 1.0; //Force per MT in pN.
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
const float_T Eta2   = Eta/5;   //Trans. drag coeff of each MTOC (pN min/mum)
const float_T kbT    = .00414;   //pN mum
const float_T D      = kbT/Eta2; //Diffusion Coefficient for each MTOC free motion. (mum^2/min)

//  Starting Centered Coordinates:
const float_T startPsi = pi/2.0;
const float_T startX   = 0;
const float_T startY   = 0;

//   Off-center Coordinates
//const float_T startPsi = pi/2.0;
//const float_T startX   = 10.0;
//const float_T startY   = 0;

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
//const float_T width = pi/8;
//const float_T centerPos = 1.24287;
//const float_T start = centerPos - width/2.0;
//const float_T end = centerPos + width/2.0;
//const float_T regionAngles[numRegions+1] = {0, start, end, 2*pi - end, 2*pi - start, 2*pi};
//const float_T regionProbabilities[numRegions] = {1,1,1,1,1};
//const float_T regionForceMultipliers[numRegions] = {1,-1,1,-1,1};

// MT density limitations
//  Only one contact per window, windows of length ~1
const size_t numberContactWindows = 128;
const float_T contactWindowAngles[numberContactWindows+1] =
  {0., 0.039985, 0.079885, 0.119585, 0.159085, 0.198485, 0.237785, 0.277085,
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
    6.21268, 6.25258, 6.28319};


 //Other parameters:
/* const size_t numberContactWindows = 643;
 const float_T contactWindowAngles[numberContactWindows] =
   {0., 0.0079, 0.0158, 0.0237, 0.0316, 0.0395, 0.0474, 0.0553, 0.0632, \
     0.0711, 0.079, 0.0869, 0.0948, 0.1027, 0.1106, 0.1185, 0.1264, \
       0.1343, 0.1422, 0.1501, 0.158, 0.1658, 0.1736, 0.1814, 0.1892, 0.197, \
       0.2048, 0.2126, 0.2204, 0.2282, 0.236, 0.2438, 0.2516, 0.2594, \
       0.2672, 0.275, 0.2828, 0.2906, 0.2984, 0.3062, 0.314, 0.3219, 0.3298, \
       0.3377, 0.3456, 0.3535, 0.3614, 0.3693, 0.3772, 0.3852, 0.3932, \
       0.4012, 0.4092, 0.4172, 0.4253, 0.4334, 0.4415, 0.4496, 0.4578, \
       0.466, 0.4742, 0.4825, 0.4908, 0.4991, 0.5074, 0.5158, 0.5242, \
       0.5327, 0.5412, 0.5497, 0.5583, 0.5669, 0.5756, 0.5843, 0.593, \
       0.6018, 0.6106, 0.6195, 0.6284, 0.6374, 0.6464, 0.6555, 0.6646, \
       0.6738, 0.683, 0.6923, 0.7016, 0.711, 0.7205, 0.73, 0.7396, 0.7492, \
       0.7589, 0.7687, 0.7785, 0.7884, 0.7983, 0.8083, 0.8184, 0.8285, \
       0.8387, 0.849, 0.8593, 0.8697, 0.8802, 0.8907, 0.9013, 0.912, 0.9227, \
       0.9335, 0.9444, 0.9554, 0.9664, 0.9775, 0.9887, 0.9999, 1.0112, \
       1.0226, 1.034, 1.0455, 1.0571, 1.0687, 1.0804, 1.0922, 1.104, 1.1159, \
       1.1279, 1.1399, 1.152, 1.1642, 1.1764, 1.1887, 1.201, 1.2134, 1.2259, \
       1.2384, 1.251, 1.2636, 1.2763, 1.289, 1.3018, 1.3146, 1.3275, 1.3404, \
       1.3533, 1.3663, 1.3793, 1.3923, 1.4054, 1.4185, 1.4316, 1.4448, \
       1.458, 1.4712, 1.4844, 1.4976, 1.5109, 1.5242, 1.5375, 1.5508, \
       1.5641, 1.5774, 1.5907, 1.604, 1.6173, 1.6306, 1.6439, 1.6571, \
       1.6703, 1.6835, 1.6967, 1.7099, 1.723, 1.7361, 1.7492, 1.7622, \
       1.7752, 1.7882, 1.8011, 1.814, 1.8269, 1.8397, 1.8525, 1.8652, \
       1.8779, 1.8905, 1.9031, 1.9156, 1.9281, 1.9405, 1.9528, 1.9651, \
       1.9773, 1.9895, 2.0016, 2.0136, 2.0256, 2.0375, 2.0493, 2.0611, \
       2.0728, 2.0844, 2.096, 2.1075, 2.1189, 2.1303, 2.1416, 2.1528, 2.164, \
       2.1751, 2.1861, 2.1971, 2.208, 2.2188, 2.2295, 2.2402, 2.2508, \
       2.2613, 2.2718, 2.2822, 2.2925, 2.3028, 2.313, 2.3231, 2.3332, \
       2.3432, 2.3531, 2.363, 2.3728, 2.3826, 2.3923, 2.4019, 2.4115, 2.421, \
       2.4305, 2.4399, 2.4492, 2.4585, 2.4677, 2.4769, 2.486, 2.4951, \
       2.5041, 2.5131, 2.522, 2.5309, 2.5397, 2.5485, 2.5572, 2.5659, \
       2.5746, 2.5832, 2.5918, 2.6003, 2.6088, 2.6173, 2.6257, 2.6341, \
       2.6424, 2.6507, 2.659, 2.6673, 2.6755, 2.6837, 2.6919, 2.7, 2.7081, \
       2.7162, 2.7243, 2.7323, 2.7403, 2.7483, 2.7563, 2.7643, 2.7722, \
       2.7801, 2.788, 2.7959, 2.8038, 2.8117, 2.8196, 2.8275, 2.8353, \
       2.8431, 2.8509, 2.8587, 2.8665, 2.8743, 2.8821, 2.8899, 2.8977, \
       2.9055, 2.9133, 2.9211, 2.9289, 2.9367, 2.9445, 2.9523, 2.9601, \
       2.9679, 2.9757, 2.9835, 2.9914, 2.9993, 3.0072, 3.0151, 3.023, \
       3.0309, 3.0388, 3.0467, 3.0546, 3.0625, 3.0704, 3.0783, 3.0862, \
       3.0941, 3.102, 3.1099, 3.1178, 3.1257, 3.1336, 3.1415, 3.1494, \
       3.1573, 3.1652, 3.1731, 3.181, 3.1889, 3.1968, 3.2047, 3.2126, \
       3.2205, 3.2284, 3.2363, 3.2442, 3.2521, 3.26, 3.2679, 3.2758, 3.2837, \
       3.2916, 3.2995, 3.3073, 3.3151, 3.3229, 3.3307, 3.3385, 3.3463, \
       3.3541, 3.3619, 3.3697, 3.3775, 3.3853, 3.3931, 3.4009, 3.4087, \
       3.4165, 3.4243, 3.4321, 3.4399, 3.4477, 3.4555, 3.4634, 3.4713, \
       3.4792, 3.4871, 3.495, 3.5029, 3.5108, 3.5187, 3.5267, 3.5347, \
       3.5427, 3.5507, 3.5587, 3.5668, 3.5749, 3.583, 3.5911, 3.5993, \
       3.6075, 3.6157, 3.624, 3.6323, 3.6406, 3.6489, 3.6573, 3.6657, \
       3.6742, 3.6827, 3.6912, 3.6998, 3.7084, 3.7171, 3.7258, 3.7345, \
       3.7433, 3.7521, 3.761, 3.7699, 3.7789, 3.7879, 3.797, 3.8061, 3.8153, \
       3.8245, 3.8338, 3.8431, 3.8525, 3.862, 3.8715, 3.8811, 3.8907, \
       3.9004, 3.9102, 3.92, 3.9299, 3.9398, 3.9498, 3.9599, 3.97, 3.9802, \
       3.9905, 4.0008, 4.0112, 4.0217, 4.0322, 4.0428, 4.0535, 4.0642, \
       4.075, 4.0859, 4.0969, 4.1079, 4.119, 4.1302, 4.1414, 4.1527, 4.1641, \
       4.1755, 4.187, 4.1986, 4.2102, 4.2219, 4.2337, 4.2455, 4.2574, \
       4.2694, 4.2814, 4.2935, 4.3057, 4.3179, 4.3302, 4.3425, 4.3549, \
       4.3674, 4.3799, 4.3925, 4.4051, 4.4178, 4.4305, 4.4433, 4.4561, \
       4.469, 4.4819, 4.4948, 4.5078, 4.5208, 4.5338, 4.5469, 4.56, 4.5731, \
       4.5863, 4.5995, 4.6127, 4.6259, 4.6391, 4.6524, 4.6657, 4.679, \
       4.6923, 4.7056, 4.7189, 4.7322, 4.7455, 4.7588, 4.7721, 4.7854, \
       4.7986, 4.8118, 4.825, 4.8382, 4.8514, 4.8645, 4.8776, 4.8907, \
       4.9037, 4.9167, 4.9297, 4.9426, 4.9555, 4.9684, 4.9812, 4.994, \
       5.0067, 5.0194, 5.032, 5.0446, 5.0571, 5.0696, 5.082, 5.0943, 5.1066, \
       5.1188, 5.131, 5.1431, 5.1551, 5.1671, 5.179, 5.1908, 5.2026, 5.2143, \
       5.2259, 5.2375, 5.249, 5.2604, 5.2718, 5.2831, 5.2943, 5.3055, \
       5.3166, 5.3276, 5.3386, 5.3495, 5.3603, 5.3711, 5.3818, 5.3924, \
       5.4029, 5.4134, 5.4238, 5.4341, 5.4444, 5.4546, 5.4647, 5.4748, \
       5.4848, 5.4947, 5.5046, 5.5144, 5.5242, 5.5339, 5.5435, 5.5531, \
       5.5626, 5.5721, 5.5815, 5.5908, 5.6001, 5.6093, 5.6185, 5.6276, \
       5.6367, 5.6457, 5.6547, 5.6636, 5.6725, 5.6813, 5.6901, 5.6988, \
       5.7075, 5.7162, 5.7248, 5.7334, 5.7419, 5.7504, 5.7589, 5.7673, \
       5.7757, 5.784, 5.7923, 5.8006, 5.8089, 5.8171, 5.8253, 5.8335, \
       5.8416, 5.8497, 5.8578, 5.8659, 5.8739, 5.8819, 5.8899, 5.8979, \
       5.9059, 5.9138, 5.9217, 5.9296, 5.9375, 5.9454, 5.9533, 5.9612, \
       5.9691, 5.9769, 5.9847, 5.9925, 6.0003, 6.0081, 6.0159, 6.0237, \
       6.0315, 6.0393, 6.0471, 6.0549, 6.0627, 6.0705, 6.0783, 6.0861, \
       6.0939, 6.1017, 6.1095, 6.1173, 6.1251, 6.133, 6.1409, 6.1488, \
       6.1567, 6.1646, 6.1725, 6.1804, 6.1883, 6.1962, 6.2041, 6.212, \
       6.2199, 6.2278, 6.2357, 6.2436, 6.2515, 6.2594, 6.2673, 6.2752, \
       6.28319};
*/  
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
