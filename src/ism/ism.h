/**
 * 
 *      @file ism.c
 *      @brief Implicit solvation model routines
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 16/08/2012
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation version 2 of the License.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 */



#define ISM_SCALE_RADIUS 0.60f
#define ISM_PROBE 1.40f
#define ISM_G_PARAM_PLUS 0.5f
#define ISM_MIN_PSURFACE 0.01f
#define ISM_APOLAR_PARAM_A 0.092f
#define ISM_APOLAR_PARAM_B 0.00542f
#define ISM_ALPHA 1.0367f
#define ISM_EPSILON 78.39f
#define ISM_EPSILON_PLUS_1 ISM_EPSILON + 1.0f
#define ISM_CONSTANT_DIELECTRIC_SCREENING (ISM_EPSILON - 1.0f) / 2.0f
#define ISM_LAMBDA1 0.013f
#define ISM_LAMBDA2 0.007f
#define ISM_PI 3.1416f
#define ISM_MAX_OVERLAPS 300
#define ISM_UNIT_ADJUST 166.0


/*                   C     O     N     H     P     S     I     Br    Cl    F     Unknown */
float ism_radii[11] = {1.70f,1.60f,1.65f,1.00f,1.90f,1.90f,1.98f,1.85f,1.80f,1.00f,1.70f};
/*                          ZN    MN    MG    K           C0    Fe    */
float ism_radii_metals[7] = {1.40f,1.70f,1.18f,1.70f,1.70f,1.70f,1.70f};


/* LCPO Parameters */

float LCPOParam1[31] = { 0.34947, 0.62334, 0.85418, 0.01911, 0.03747, 0.11704, 0.47838, 0.39092,
			 0.49584, 0.60049, 0.49720, 0.30541, 0.48891, 0.00150, 0.60586, 0.53086,
			 0.46679, 0.04221, 0.02967, 0.01094, 0.14807, 0.29507, 0.05376, 0.20262,
			 0.33431, 0.24876, 0.42885, 0.81175, 0.03898, 0.48748, 0.45878 };

float LCPOParam2[31] = {-0.06810,-0.13735,-0.22673,-0.00357,-0.00646,-0.02228,-0.11682,-0.10826,
			-0.11532,-0.14697,-0.12383,-0.06209,-0.09816,-0.00019,-0.11790,-0.11910,
			-0.07256, 0.03157,-0.00552,-0.00207,-0.03379,-0.06370,-0.00961,-0.03672,
			-0.05910,-0.04630,-0.09276,-0.17740,-0.00510, 0.02120,-0.12281};
	
float LCPOParam3[31] = {-0.00028,-0.00007,-0.00059, 0.00001, 0.00001, 0.00005,-0.00006, 0.00007,
			-0.00003,-0.00064,-0.00029,-0.00019,-0.00010, 0.00000,-0.00047,-0.00039,
			 0.00003,-0.00037, 0.00001, 0.00000,-0.00006,-0.00020, 0.00002,-0.00017,
			-0.00139,-0.00091,-0.00066,-0.00083, 0.00006,-0.00199,-0.00091};


float LCPOParam4[31] = { 0.00007, 0.00012, 0.00026, 0.00000, 0.00000, 0.00001, 0.00011, 0.00010,
			 0.00010, 0.00017, 0.00013, 0.00006, 0.00008, 0.00000, 0.00012, 0.00013,
			 0.00003,-0.00004, 0.00000, 0.00000, 0.00003, 0.00007, 0.00001, 0.00003,
			 0.00013, 0.00009, 0.00013, 0.00020, 0.00000,-0.00002, 0.00020};


/****
 *
 * ga rules
 *
 * HO +0.8
 * HN(N4) +0.765
 * HN(guanidinio) +0.560
 * HN(ring) 0.821
 * HN(amido) 0.327
 *
 * OH +0.1
 * O 0.00
 * O- -0.256
 * N(acceptor) -0.189
 *
 */


typedef struct{
  /* Desolvation */
  /* Protein */
  int *numOverlaps_prot;
  int *buried_prot;
  int **OverlapsMatrix_prot;
  float **AijMatrix_prot;
  float **lcpoComps_prot;
  float *total_surface_prot;
  float *atomSASA_prot;
  float *atomSASAComplex_prot;
  float *atomSASAH_prot;
  /* Ligand */
  float **matrix;
  int *numOverlaps;
  int *buried; 
  int **OverlapsMatrix;
  float **AijMatrix;
  float **lcpoComps;
  float *total_surface;
  float *atomSASA; 
  float *atomSASAComplex;
  float *atomSASAH;
  /* Complex */
  int *buried_lig_complex; 
  int *buried_prot_complex;
  float **AijMatrixLig; 
  float **AijMatrixProt;
  int *numOverlapsLigComplex;
  int *numOverlapsProtComplex;
  int **OverlapsMatrixLig;
  int **OverlapsMatrixProt;
  /* H-bonds */
  float *ga;
  float *ga_prot;
  /* Intra-molecular bonds */
  int *numOverlaps_hbond_prot;
  int *buried_hbond_prot;
  int **OverlapsMatrix_hbond_prot;
  float **AijMatrix_hbond_prot;
  int *numOverlaps_hbond;
  int *buried_hbond;
  int **OverlapsMatrix_hbond;
  float **AijMatrix_hbond;
  /* Inter-molecular bonds */
  int *buried_lig_complex_hbond;
  int *buried_prot_complex_hbond;
  float **AijMatrixLig_hbond; 
  float **AijMatrixProt_hbond;
  int *numOverlapsLigComplex_hbond;
  int *numOverlapsProtComplex_hbond;
  int **OverlapsMatrixLig_hbond;
  int **OverlapsMatrixProt_hbond;
} ISM_COMPLEX;
