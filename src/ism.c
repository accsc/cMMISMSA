/**
 *
 *      @file ism.c
 *      @brief Implicit solvation model routines
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @date 09/08/2012
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
 *
 *	30/10/2012. I finally have found the bug on desolvation calcs
 */

#include <stdio.h>
#include <stdlib.h>
#include <ism.h>

int ism_get_distance_matrix(MOL2 *mymol, float ***mymatrix);


/**
 *	@brief Calculate the apolar contribution for a complex using the Implicit Solvation Model
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param lig Ligand MOL2 structure
 *	@param prot Protein MOL2 structure
 *	@param atomSASA Atomic Accesible Surfaces for the free ligand
 *	@param atomSASAComplex Atomic Accessible Surfaces for the ligand in complex
 *	@param totalSurface Total surface for the free ligand without taking into account overlaps
 *	@param atomSASA_prot Accesible Surfaces for the free protein
 *	@param atomSASAComplex_prot Atomic Accessible Surfaces for the protein in complex
 *	@param totalSurface_prot Total surface for the free protein without taking into account overlaps
 *	@return Energy value
 *
 */
float ism_get_apolar_energy(MOL2 *lig, MOL2 *prot, float *atomSASA, float *atomSASAComplex, float *totalSurface,
			   float *atomSASA_prot, float *atomSASAComplex_prot, float *totalSurface_prot)
{

	int i = 0;
	float areaA = 0.0f, areaAB = 0.0f, areaB = 0.0f, tmpApolarA = 0.0f, tmpApolarB = 0.0f, tmpApolarAB = 0.0f;
	float res = 0.0f;
	
	for( i = 0; i < lig->n_atoms; i++)
	{
		areaB += atomSASA[i] * totalSurface[i];
		areaAB += atomSASAComplex[i] * totalSurface[i];
	}

        for( i = 0; i < prot->n_atoms; i++)
        {
                areaA += atomSASA_prot[i] * totalSurface_prot[i];
                areaAB += atomSASAComplex_prot[i] * totalSurface_prot[i];
        }

        tmpApolarA = ISM_APOLAR_PARAM_A + areaA * ISM_APOLAR_PARAM_B;
        tmpApolarB = ISM_APOLAR_PARAM_A + areaB * ISM_APOLAR_PARAM_B;
        tmpApolarAB = ISM_APOLAR_PARAM_A + areaAB * ISM_APOLAR_PARAM_B;

        res = tmpApolarAB - tmpApolarA - tmpApolarB;

	return res;

}

/**
 *	@brief Get polar contribution of desolvation for protein
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Protein MOL2 structure
 *	@param lig Ligand MOL2 structure
 *	@param atomSASA Atomic accesible surfaces for the free protein
 *	@param atomSASAComplex Atomic accesible surfaces for the protein in complex
 *	@param atomSASAHs Fraction of proton buried in acceptor for hbonds
 *	@return Energy value
 */
float ism_get_desolvation_energy_prot(MOL2 *mymol, MOL2 *lig, float *atomSASA, float *atomSASAComplex, float *atomSASAHs)
{

	float res = 0.0f;
	float lambda = 0.0f;
	int i = 0;

	float rsu = 0.0f, rsc = 0.0f, dRsu = 0.0f, dRsc = 0.0f, r1 = 0.0f;
	float tmpDesolvation = 0.0f;

	float Hpar = 0.0f, Gpar = 0.0f;

	lambda = ISM_LAMBDA2;

	for( i = 0; i < lig->n_atoms; i++)
	{
		if( lig->gaff_types[i] == N4)
		{
			lambda = ISM_LAMBDA1;
			break;
		}
	}


	for( i = 0; i < mymol->n_atoms; i++)
	{
	   if( mymol->ism_selection[i] == 1)
           {
                if( mymol->atoms[i] <= 100)
                  r1 = ism_radii[ mymol->atoms[i]-1];
                else
                  r1 = ism_radii_metals[ mymol->atoms[i]-100];

		r1 = r1 * ISM_SCALE_RADIUS; /* Covalent radius */

		if( mymol->pcharges[i] > 0)
		  Hpar = 0.85f;
		else
		  Hpar = 0.35f;

		Gpar = Hpar + ISM_G_PARAM_PLUS;


		rsu = ((r1 + Hpar) * atomSASA[i]) + ((r1 + Gpar) * (1.0f - atomSASA[i]));

		rsc = ((r1 + Hpar) * atomSASAComplex[i]) + (1.44f * atomSASAHs[i] ) + ((r1 + Gpar) * (1.0f - atomSASAComplex[i]- atomSASAHs[i]));

		dRsu = (ISM_EPSILON_PLUS_1 / (1 + ISM_CONSTANT_DIELECTRIC_SCREENING * exp(-lambda * ISM_EPSILON_PLUS_1 * rsu))) - 1.0f;
		dRsc = (ISM_EPSILON_PLUS_1 / (1 + ISM_CONSTANT_DIELECTRIC_SCREENING * exp(-lambda * ISM_EPSILON_PLUS_1 * rsc))) - 1.0f;

		tmpDesolvation = (mymol->pcharges[i] * mymol->pcharges[i]) * ((1.0f / (dRsc * rsc)) - (1.0f / (dRsu * rsu)) + (1.0f / rsu) - (1.0f / rsc));

		res = res + tmpDesolvation;

	   }
	}

	res *= ISM_UNIT_ADJUST;

	return res;

}


/**
 *     @brief Get polar contribution of desolvation for ligand
 *     @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *     @param mymol Ligand MOL2 structure
 *     @param lig Ligand MOL2 structure
 *     @param atomSASA Atomic accesible surfaces for the free ligand
 *     @param atomSASAComplex Atomic accesible surfaces for the ligand in complex
 *     @param atomSASAHs Fraction of proton buried in acceptor for hbonds
 *     @return Energy value
 */
float ism_get_desolvation_energy(MOL2 *mymol, MOL2 *lig, float *atomSASA, float *atomSASAComplex, float *atomSASAHs, float *ga)
{

	float res = 0.0f;
	float lambda = 0.0f;
	int i = 0;

	float rsu = 0.0f, rsc = 0.0f, dRsu = 0.0f, dRsc = 0.0f, r1 = 0.0f;
	float tmpDesolvation = 0.0f;

	float Hpar = 0.0f, Gpar = 0.0f;

	lambda = ISM_LAMBDA1;

/*	for( i = 0; i < lig->n_atoms; i++)
	{
		if( lig->gaff_types[i] == N4)
		{
			lambda = ISM_LAMBDA2;
			break;
		}
	}*/


	for( i = 0; i < mymol->n_atoms; i++)
	{
                if( mymol->atoms[i] <  100)
                  r1 = ism_radii[ mymol->atoms[i]-1];
                else
                  r1 = ism_radii_metals[ mymol->atoms[i]-100];

		r1 = r1 * ISM_SCALE_RADIUS; /* Covalent radius */

		if( mymol->pcharges[i] > 0)
		  Hpar = 0.85f;
		else
		  Hpar = 0.35f;

		Gpar = Hpar + ISM_G_PARAM_PLUS;

/*		if( atomSASAHs[i] > 0)
			fprintf(stderr,"%i %f %f %f | %f %f\n",i,atomSASA[i],atomSASAComplex[i],atomSASAHs[i],Hpar,Gpar);*/

		if( ga[i] <= 0.0f)
			atomSASAHs[i] = 0.0f;

		/* Effective Born radii */
                rsu = ((r1 + Hpar) * atomSASA[i]) + ((r1 + Gpar) * (1.0f - atomSASA[i]));
                rsc = ((r1 + Hpar) * atomSASAComplex[i]) + ((r1+ga[i]) * atomSASAHs[i] ) + ((r1 + Gpar) * (1.0f - atomSASAComplex[i] - atomSASAHs[i]));



/*		if( mymol->gaff_types[i] == N4)
			lambda = ISM_LAMBDA2;
		else
                        lambda = ISM_LAMBDA1;*/

		dRsu =  ((78.39f+1.0f) / (1.0f + ((78.39f-1.0f)/2.0f)*exp(-lambda *(78.39f+1.0f)*rsu)))-1.0f;
		dRsc = ((78.39f+1.0f) / (1.0f + ((78.39f-1.0f)/2.0f)*exp(-lambda *(78.39f+1.0f)*rsc)))-1.0f;

		tmpDesolvation = (mymol->pcharges[i] * mymol->pcharges[i]) * ((1.0f / (dRsc * rsc)) - (1.0f / (dRsu * rsu)) + (1.0f / rsu) - (1.0f / rsc));

		res = res + tmpDesolvation;

	}

	res *= ISM_UNIT_ADJUST;

	return res;

}

/**
 *	@brief Calculate SASA for the protein in complex
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Protein MOL2 structure
 *	@param total_surface Total surface without overlaps
 *	@param buried Buried atoms flag in the free protein
 *	@param buried_complex Buried atoms flag in the complex
 *	@param numOverlaps Overlaps count for the free protein
 *	@param numOverlapsComplex Overlaps count for the protein in complex
 *	@param numOverlapsB Overalps count for the ligand
 *	@param numOverlapsComplexB Overlaps count for the ligand in complex
 *	@param AijMatrix Overlaps SASA matrix for the free protein
 *	@param AijMatrixComplex Overlaps SASA matrix for the protein in complex
 *	@param AijMatrixB Overlaps SASA matrix for the ligand
 *	@param AijMatrixComplexB Overlaps SASA matrix for the ligand in complex
 *	@param OverlapsMatrix Overlaps indexes matrix for the free protein
 *	@param OverlapsMatrixComplex Overlaps indexes matrix for the protein in complex
 *	@param OverlapsMatrixB Overlaps indexes matrix for the ligand
 *	@param OverlapsMatrixComplexB Overlaps indexes matrix for the ligand in complex
 *	@param lcpoComps LCPO model components
 *	@param myatomSASA Resulting complex SASA
 *	@return 0 on success
 */
int ism_get_complex_sasa_prot(MOL2 *mymol, float *total_surface, int *buried, int *buried_complex,
		      int *numOverlaps, int *numOverlapsComplex, int *numOverlapsB, int *numOverlapsComplexB,
		      float **AijMatrix, float **AijMatrixComplex, float **AijMatrixB, float **AijMatrixComplexB,
		      int **OverlapsMatrix, int **OverlapsMatrixComplex, int **OverlapsMatrixB, int **OverlapsMatrixComplexB,
		      float **lcpoComps, float **myatomSASA)
{


	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float *atomSASA = NULL;

	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{

	  if( mymol->ism_selection != 0)
          {
		tmpSASA = total_surface[i];
		if( buried[i] == 1 || buried_complex[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] + numOverlapsComplex[i] == 1){
				if( numOverlaps[i] == 1)
				 tmpSASA = tmpSASA - AijMatrix[0][i];
				else
				 tmpSASA = tmpSASA - AijMatrixComplex[0][i];
			}else if( (numOverlaps[i] + numOverlapsComplex[i]) > 1){
				lcpo[0] = lcpoComps[i][0];
				lcpo[1] = lcpoComps[i][1];
				lcpo[2] = lcpoComps[i][2];
				lcpo[3] = lcpoComps[i][3];

				for( j = 0; j < numOverlaps[i]; j++)
				{
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlapsComplex[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixComplex[k][idxOverlap];
						for( l = 0; l < numOverlapsComplex[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrixComplex[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrixComplex[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrixComplex[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				for( j = 0; j < numOverlapsComplex[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrixComplex[j][i];
                                        scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrixComplex[j][i];
                                        for( k = 0; k < numOverlapsB[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixB[k][idxOverlap];
                                                for( l = 0; l < numOverlapsComplex[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrixComplex[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

					}

                                        for( k = 0; k < numOverlapsComplexB[idxOverlap]; k++)
                                        {
						subidxOverlap = OverlapsMatrixComplexB[k][idxOverlap];
                                                for( l = 0; l < numOverlaps[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrix[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

                                        }

                                        lcpo[3] = lcpo[3] + AijMatrixComplex[j][i] * scdlcpo4;

				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

/*				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];*/
			}



		}


		if( tmpSASA > 0)
			tmpSASA = tmpSASA / total_surface[i];
		else
		 	tmpSASA = 0.0f;
	

		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;
          }

	}


	*myatomSASA = atomSASA;
	return 0;

}

/**
 *      @brief Calculate SASA for the ligand in complex
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @param mymol Protein MOL2 structure
 *      @param total_surface Total surface without overlaps
 *      @param buried Buried atoms flag in the free ligand
 *      @param buried_complex Buried atoms flag in the complex
 *      @param numOverlaps Overlaps count for the free ligand
 *      @param numOverlapsComplex Overlaps count for the ligand in complex
 *      @param numOverlapsB Overalps count for the protein
 *      @param numOverlapsComplexB Overlaps count for the protein in complex
 *      @param AijMatrix Overlaps SASA matrix for the free ligand
 *      @param AijMatrixComplex Overlaps SASA matrix for the ligand in complex
 *      @param AijMatrixB Overlaps SASA matrix for the protein
 *      @param AijMatrixComplexB Overlaps SASA matrix for the protein in complex
 *      @param OverlapsMatrix Overlaps indexes matrix for the free ligand
 *      @param OverlapsMatrixComplex Overlaps indexes matrix for the ligand in complex
 *      @param OverlapsMatrixB Overlaps indexes matrix for the protein
 *      @param OverlapsMatrixComplexB Overlaps indexes matrix for the protein in complex
 *      @param lcpoComps LCPO model components
 *      @param myatomSASA Resulting complex SASA
 *      @return 0 on success
 */
int ism_get_complex_sasa(MOL2 *mymol, float *total_surface, int *buried, int *buried_complex,
		      int *numOverlaps, int *numOverlapsComplex, int *numOverlapsB, int *numOverlapsComplexB,
		      float **AijMatrix, float **AijMatrixComplex, float **AijMatrixB, float **AijMatrixComplexB,
		      int **OverlapsMatrix, int **OverlapsMatrixComplex, int **OverlapsMatrixB, int **OverlapsMatrixComplexB,
		      float **lcpoComps, float **myatomSASA)
{


	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float *atomSASA = NULL;

	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{

		tmpSASA = total_surface[i];
		if( buried[i] == 1 || buried_complex[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] + numOverlapsComplex[i] == 1){
				if( numOverlaps[i] == 1)
				 tmpSASA = tmpSASA - AijMatrix[0][i];
				else
				 tmpSASA = tmpSASA - AijMatrixComplex[0][i];
			}else if( (numOverlaps[i] + numOverlapsComplex[i]) > 1){
				lcpo[0] = lcpoComps[i][0];
				lcpo[1] = lcpoComps[i][1];
				lcpo[2] = lcpoComps[i][2];
				lcpo[3] = lcpoComps[i][3];

				for( j = 0; j < numOverlaps[i]; j++)
				{
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlapsComplex[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixComplex[k][idxOverlap];
						for( l = 0; l < numOverlapsComplex[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrixComplex[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrixComplex[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrixComplex[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				for( j = 0; j < numOverlapsComplex[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrixComplex[j][i];
                                        scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrixComplex[j][i];
                                        for( k = 0; k < numOverlapsB[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixB[k][idxOverlap];
                                                for( l = 0; l < numOverlapsComplex[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrixComplex[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

					}

                                        for( k = 0; k < numOverlapsComplexB[idxOverlap]; k++)
                                        {
						subidxOverlap = OverlapsMatrixComplexB[k][idxOverlap];
                                                for( l = 0; l < numOverlaps[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrix[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

                                        }

                                        lcpo[3] = lcpo[3] + (AijMatrixComplex[j][i] * scdlcpo4);

				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

/*				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];*/


			}



		}
		if( tmpSASA > 0)
			tmpSASA = tmpSASA / total_surface[i];
		else
		 	tmpSASA = 0.0f;
	


		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;

	}



	*myatomSASA = atomSASA;

	return 0;

}



/**
 *	@brief Calculate SASA for the free ligand
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Molecule MOL2 structure
 *	@param total_surface Total surface without overlaps
 *	@param buried Buried atom flag
 *	@param numOverlaps Number of overlaps per atom
 *	@param AijMatrix SASA contribution matrix
 *	@param OverlapsMatrix Overlaps index matrix
 *	@param mylcpoComps LCPO contributions
 *	@param myatomSASA Free SASA result
 *	@return 0 on success
 */

int ism_get_free_sasa(MOL2 *mymol, float *total_surface, int *buried, 
		      int *numOverlaps, float **AijMatrix, int **OverlapsMatrix,
		      float ***mylcpoComps, float **myatomSASA)
{
	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float **lcpoComps = NULL;
	float *atomSASA = NULL;

	lcpoComps = *mylcpoComps;
	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{

		tmpSASA = total_surface[i];
		if( buried[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] == 1){
				tmpSASA = tmpSASA - AijMatrix[0][i];
			}else if( numOverlaps[i] > 1){
				lcpo[0] = total_surface[i];
				lcpo[1] = 0.0f;
				lcpo[2] = 0.0f;
				lcpo[3] = 0.0f;

				for( j = 0; j < numOverlaps[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrix[j][i];
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlaps[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrix[k][idxOverlap];
						for( l = 0; l < numOverlaps[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrix[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrix[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrix[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];

			}



		}

		if( tmpSASA > 0)
			tmpSASA = tmpSASA / total_surface[i];
		else
		 	tmpSASA = 0.0f;
	

		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;

	}


	*myatomSASA = atomSASA;

	return 0;

}


/**
 *      @brief Calculate SASA for the free protein
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @param mymol Molecule MOL2 structure
 *      @param total_surface Total surface without overlaps
 *      @param buried Buried atom flag
 *      @param numOverlaps Number of overlaps per atom
 *      @param AijMatrix SASA contribution matrix
 *      @param OverlapsMatrix Overlaps index matrix
 *      @param mylcpoComps LCPO contributions
 *      @param myatomSASA Free SASA result
 *      @return 0 on success
 */
int ism_get_free_sasa_prot(MOL2 *mymol, float *total_surface, int *buried, 
		      int *numOverlaps, float **AijMatrix, int **OverlapsMatrix,
		      float ***mylcpoComps, float **myatomSASA)
{
	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float **lcpoComps = NULL;
	float *atomSASA = NULL;

	lcpoComps = *mylcpoComps;
	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{
	   if( mymol->ism_selection[i] != 0 )
           {

		tmpSASA = total_surface[i];
		if( buried[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] == 1){
				tmpSASA = tmpSASA - AijMatrix[0][i];
			}else if( numOverlaps[i] > 1){
				lcpo[0] = total_surface[i];
				lcpo[1] = 0.0f;
				lcpo[2] = 0.0f;
				lcpo[3] = 0.0f;

				for( j = 0; j < numOverlaps[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrix[j][i];
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlaps[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrix[k][idxOverlap];
						for( l = 0; l < numOverlaps[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrix[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrix[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrix[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];
			}



		}

		if( tmpSASA > 0)
			tmpSASA = tmpSASA / total_surface[i];
		else
		 	tmpSASA = 0.0f;
	

		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;

           }
	}

	*myatomSASA = atomSASA;

	return 0;

}


/**
 *	@brief Calculate protein-protein overlaps 
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Protein MOL2 structure
 *	@param mynumOverlaps Vector of number of atomic overlaps
 *	@param myburied Atom buried flag
 *	@param myAijMatrix Overlap contributions
 *	@param myOverlapsMatrix Overlap indices
 *	@return 0 on success
 */
int ism_getInternalOverlaps_prot(MOL2 *mymol, int **mynumOverlaps, int **myburied, float ***myAijMatrix, int ***myOverlapsMatrix)
{

        int i = 0, j = 0, k = 0;
        float r1 = 0.0f, r2 = 0.0f;
        float deltaRadii = 0.0f;
	float dx = 0.0f, dy = 0.0f, dz = 0.0, distance = 0.0f;

        int *numOverlaps = NULL, *buried = NULL, **OverlapsMatrix = NULL;
	float **AijMatrix = NULL;
        numOverlaps = *mynumOverlaps;
        buried = *myburied;
        AijMatrix = *myAijMatrix;
        OverlapsMatrix = *myOverlapsMatrix;


        for( i = 0; i < mymol->n_atoms; i++)
        {
		if( mymol->ism_selection[i] == 1)
		{
	                numOverlaps[i] = 0;
	                buried[i] = 0;
		}
        }

        for( i = 0; i < (mymol->n_atoms-1); i++)
        {
                if( mymol->atoms[i] <= 100)
                  r1 = ism_radii[ mymol->atoms[i]-1];
                else
                  r1 = ism_radii_metals[ mymol->atoms[i]-100];

                r1 += ISM_PROBE;
                for( j = i + 1; j < mymol->n_atoms; j++)
                {

		    if( mymol->ism_selection[i] != 0 && mymol->ism_selection[j] != 0)
		    {
                        if( mymol->atoms[j] <= 100)
                          r2 = ism_radii[ mymol->atoms[j]-1];
                        else
                          r2 = ism_radii_metals[ mymol->atoms[j]-100];

                        r2 += ISM_PROBE;

			dx = mymol->x[i] - mymol->x[j];
                        dy = mymol->y[i] - mymol->y[j];
                        dz = mymol->z[i] - mymol->z[j];
			distance = (dx * dx ) + (dy * dy) + (dz * dz);
			distance = sqrt(distance);


                        if( distance <= (r1 + r2))
                        {
                                numOverlaps[i] = numOverlaps[i] + 1;
                                if( numOverlaps[i] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                numOverlaps[j] = numOverlaps[j] + 1;
                                if( numOverlaps[j] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                OverlapsMatrix[numOverlaps[i]-1][i] = j;
                                OverlapsMatrix[numOverlaps[j]-1][j] = i;
                                AijMatrix[ numOverlaps[i]-1][i] = 2.0f * ISM_PI * r1 *
                                  (r1 - (distance / 2.0f) - ((r1*r1) - (r2*r2)) /
                                  ( 2.0f*distance));


                                AijMatrix[ numOverlaps[j]-1][j] = 2.0f * ISM_PI * r2 *
                                  (r2 - (distance / 2.0f) - ((r2*r2) - (r1*r1)) /
                                  ( 2.0f*distance));

                                deltaRadii = r1 - r2;

                                if( distance - fabs(deltaRadii) <= 0.0f)
                                {
                                        if( deltaRadii <= 0.0)
                                                buried[i] = 1;
                                        if( deltaRadii >= 0.0)
                                                buried[j] = 1;
                                }


                        }
		    }
                }
        }

        return 0;
}

/**
 *      @brief Calculate ligand-ligand overlaps 
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @param mymol Protein MOL2 structure
 *      @param mynumOverlaps Vector of number of atomic overlaps
 *      @param myburied Atom buried flag
 *      @param myAijMatrix Overlap contributions
 *      @param myOverlapsMatrix Overlap indices
 *      @return 0 on success
 */

int ism_getInternalOverlaps(MOL2 *mymol, float **distance, int **mynumOverlaps, int **myburied, float ***myAijMatrix, int ***myOverlapsMatrix)
{

	int i = 0, j = 0, k = 0;
	float r1 = 0.0f, r2 = 0.0f;
	float deltaRadii = 0.0f;

	int *numOverlaps = NULL, *buried = NULL, **OverlapsMatrix = NULL;
	float **AijMatrix = NULL;

	numOverlaps = *mynumOverlaps;
	buried = *myburied;
	AijMatrix = *myAijMatrix;
	OverlapsMatrix = *myOverlapsMatrix;

	for( i = 0; i < mymol->n_atoms; i++)
	{
		numOverlaps[i] = 0;
		buried[i] = 0;
	}
	
	for( i = 0; i < (mymol->n_atoms-1); i++)
	{
		if( mymol->atoms[i] <= 100)
		  r1 = ism_radii[ mymol->atoms[i]-1];
		else
		  r1 = ism_radii_metals[ mymol->atoms[i]-100];

		r1 += ISM_PROBE;

		for( j = i + 1; j < mymol->n_atoms; j++)
		{
	                if( mymol->atoms[j] <= 100)
	                  r2 = ism_radii[ mymol->atoms[j]-1];
	                else
	                  r2 = ism_radii_metals[ mymol->atoms[j]-100];
			
			r2 += ISM_PROBE;

			if( distance[j][i] <= (r1 + r2))
			{
				numOverlaps[i] = numOverlaps[i] + 1;
				if( numOverlaps[i] > ISM_MAX_OVERLAPS)
				{
					fprintf(stderr,"Number of overlaps exceed the allocated memroy\n");
					fflush(stderr);
					return -1;
				}

                                numOverlaps[j] = numOverlaps[j] + 1;
                                if( numOverlaps[j] > ISM_MAX_OVERLAPS)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memroy\n");
                                        fflush(stderr);
                                        return -1;
                                }

				OverlapsMatrix[numOverlaps[i]-1][i] = j;
                                OverlapsMatrix[numOverlaps[j]-1][j] = i;
				AijMatrix[ numOverlaps[i]-1][i] = 2.0f * ISM_PI * r1 * 
				  (r1 - (distance[j][i] / 2.0f) - ((r1*r1) - (r2*r2)) / 
				  ( 2.0f*distance[j][i]));

                                AijMatrix[ numOverlaps[j]-1][j] = 2.0f * ISM_PI * r2 * 
                                  (r2 - (distance[i][j] / 2.0f) - ((r2*r2) - (r1*r1)) / 
                                  ( 2.0f*distance[i][j]));
	
				deltaRadii = r1 - r2;

				if( distance[j][i] - fabs(deltaRadii) <= 0.0f)
				{
					if( deltaRadii <= 0.0)
						buried[i] = 1;
					if( deltaRadii >= 0.0)
						buried[j] = 1;
				}


			}	

		}
	}

	return 0;

}


/**
 *	@brief Calculate ligand-protein overlaps
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param lig Ligand MOL2 structure
 *	@param prot Protein MOL2 structure
 *	@param myburied_lig_complex Atom buried flag for ligand
 *	@param myburied_prot_complex Atom buried flag for protein
 *	@param mynumOverlapsLigComplex Vector of number of overlaps for ligand in complex
 *	@param mynumOverlapsProtComplex Vector of number of overlaps for protein in complex
 *	@param myOverlapsMatrixLig Overlaps indices for ligand in complex
 *	@param myOverlapsMatrixProt Overlaps indices for protein in complex
 *	@param myAijMatrixLig Matrix of overlap contributions for ligand in complex
 *	@param myAijMatrixProt Matrix of overlap contributions for protein in complex
 *	@return 0 on success
 */
int ism_getExternalOverlaps(MOL2 *lig, MOL2 *prot, int **myburied_lig_complex, int **myburied_prot_complex, int **mynumOverlapsLigComplex, int **mynumOverlapsProtComplex, int ***myOverlapsMatrixLig, int ***myOverlapsMatrixProt, float ***myAijMatrixLig, float ***myAijMatrixProt)
{

	int *buried_lig_complex = NULL, *buried_prot_complex = NULL;
        int *numOverlapsLigComplex = NULL, *numOverlapsProtComplex = NULL;
	int **OverlapsMatrixLig = NULL, **OverlapsMatrixProt;
	float **AijMatrixLig = NULL, **AijMatrixProt = NULL;

        float r1 = 0.0f, r2 = 0.0f;
        float deltaRadii = 0.0f;
        float dx = 0.0f, dy = 0.0f, dz = 0.0, distance = 0.0f;


	int i = 0, j = 0;

	buried_lig_complex = *myburied_lig_complex;
	buried_prot_complex = *myburied_prot_complex;
	numOverlapsLigComplex = *mynumOverlapsLigComplex;
	numOverlapsProtComplex = *mynumOverlapsProtComplex;

	OverlapsMatrixLig = *myOverlapsMatrixLig;
	OverlapsMatrixProt = *myOverlapsMatrixProt;
	AijMatrixLig = *myAijMatrixLig;
	AijMatrixProt = *myAijMatrixProt;

	for( i = 0; i < lig->n_atoms; i++)
	{
		buried_lig_complex[i] = 0;
		numOverlapsLigComplex[i] = 0;
		
	}

        for( i = 0; i < prot->n_atoms; i++)
        {
                buried_prot_complex[i] = 0;
                numOverlapsProtComplex[i] = 0;

        }

	for( i = 0; i < prot->n_atoms; i++)
	{
		if( prot->ism_selection[i] != 1)
			continue;
		
                if( prot->atoms[i] < 100)
                  r1 = ism_radii[ prot->atoms[i]-1];
                else
                  r1 = ism_radii_metals[ prot->atoms[i]-100];

                r1 += ISM_PROBE;


		for( j = 0; j < lig->n_atoms; j++)
		{
                        if( lig->atoms[j] < 100)
                          r2 = ism_radii[ lig->atoms[j]-1];
                        else
                          r2 = ism_radii_metals[ lig->atoms[j]-100];

                        r2 += ISM_PROBE;

                        dx = prot->x[i] - lig->x[j];
                        dy = prot->y[i] - lig->y[j];
                        dz = prot->z[i] - lig->z[j];
                        distance = (dx * dx ) + (dy * dy) + (dz * dz);
                        distance = sqrt(distance);

                        if( distance < (r1 + r2))
                        {
                                numOverlapsLigComplex[j] = numOverlapsLigComplex[j] + 1;
                                if( numOverlapsLigComplex[j] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                numOverlapsProtComplex[i] = numOverlapsProtComplex[i] + 1;
                                if( numOverlapsProtComplex[i] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                OverlapsMatrixProt[numOverlapsProtComplex[i]-1][i] = j;
                                OverlapsMatrixLig[numOverlapsLigComplex[j]-1][j] = i;
                                AijMatrixProt[ numOverlapsProtComplex[i]-1][i] = 2.0f * ISM_PI * r1 *
                                  (r1 - (distance / 2.0f) - ((r1*r1) - (r2*r2)) /
                                  ( 2.0f*distance));


                                AijMatrixLig[ numOverlapsLigComplex[j]-1][j] = 2.0f * ISM_PI * r2 *
                                  (r2 - (distance / 2.0f) - ((r2*r2) - (r1*r1)) /
                                  ( 2.0f*distance));

                                deltaRadii = r1 - r2;


                                if( distance - fabs(deltaRadii) <= 0.0f)
                                {
                                        if( deltaRadii <= 0.0)
                                                buried_prot_complex[i] = 1;
                                        if( deltaRadii >= 0.0)
                                                buried_lig_complex[j] = 1;
                                }
			}
		}
	}

	return 0;

}

/**
 *	@brief Select only atoms in the grid for Implicit solvent calculation
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Pointer to protein MOL2 structure
 *	@param max_grid Max grid coordinates
 *	@param min_grid Min grid coordinates
 *	@return 0 on success
 */
int ism_select_atoms(MOL2 **mymol, float *max_grid, float *min_grid)
{
	MOL2 *mol = NULL;
	int i = 0;
	mol = *mymol;


	for( i = 0; i < mol->n_atoms; i++)
	{
		mol->ism_selection[i] = 0;

	        if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
		{
			mol->ism_selection[i] = 1;
		}
	}


	*mymol = mol;
	return 0;
}

/**
 *	@brief Calculate the distance matrix for a molecule
 *	@param mymol Molecule MOL2 structure
 *	@param mymatrix Matrix
 *	@return 0 on success
 */
int ism_get_distance_matrix(MOL2 *mymol, float ***mymatrix)
{
	float **matrix = NULL;
	int i = 0, j = 0;
	float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;

	matrix = *mymatrix;


	for( i = 0; i < mymol->n_atoms; i++)
	{
		matrix[i][i] = 0.0f;
		for( j = i + 1; j < mymol->n_atoms; j++)
		{
			dx = mymol->x[i] - mymol->x[j];
                        dy = mymol->y[i] - mymol->y[j];
                        dz = mymol->z[i] - mymol->z[j];
			dist = ( dx * dx) + ( dy * dy ) + ( dz * dz);
			dist = sqrtf(dist);
			matrix[i][j] = dist;
			matrix[j][i] = dist;
		}
	}

	*mymatrix = matrix;
	return 0;
}


/**
 *	@brief Atom typing for implicit solvation model
 *	@param mymol Pointer to molecule MOL2 structure
 *	@return 0 on success
 */
int ism_typing(MOL2 **mymol)
{
	MOL2 *mols = NULL;
        int vecinos[8];
        int vecinos2[8];
	int k = -1, k_heavy = -1, i = 0, j = 0;

	mols = *mymol;

	
        for (i = 0; i < mols->n_atoms; ++i)
               mols->ism_types[i] = -1;


        for (i = 0; i < mols->n_atoms; ++i) {

                k = 0;
		k_heavy = 0;
                for (j = 0; j < mols->n_bonds; ++j) {
			if ( mols->bond_a1[j] == (i + 1)) {
                                k++;
				if( mols->atoms[mols->bond_a2[j] - 1] != 4 )
					k_heavy++;
                        }else if ( mols->bond_a2[j] == (i + 1)) {
                                k++;
                                if( mols->atoms[mols->bond_a1[j] - 1] != 4 )
                                        k_heavy++;
                        }
                }


                if ( mols->atoms[i] == 4 ) {  /* H */
                	mols->ism_types[i] = 1;
                }else if ( mols->atoms[i] == 9 ) {  /* Cl */
                        mols->ism_types[i] = 2;
                }else if ( mols->atoms[i] == 105 ) {  /* Ca or C0 */
                        mols->ism_types[i] = 3;
                }else if ( mols->atoms[i] == 5 ) {  /* P */
                        mols->ism_types[i] = 4;
                }else if ( mols->atoms[i] == 3 ) {  /* N */
			if ( k >= 4)
			{
				if( k_heavy == 1)
				  mols->ism_types[i] = 7;
				else if ( k_heavy == 2)
                                  mols->ism_types[i] = 6;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 5;
				else
                                  mols->ism_types[i] = 12;

			}else if( k == 3){

                               if( k_heavy == 1)
                                  mols->ism_types[i] = 10;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 9;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 8;
				else
                                  mols->ism_types[i] = 10;

			}else if( k == 2){
                              if( k_heavy == 1)
                                  mols->ism_types[i] = 12;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 11;
                                else
                                  mols->ism_types[i] = 12;


			}else{
				mols->ism_types[i] = 12;
			}

		}else if ( mols->atoms[i] == 2 ) {  /* O */
                        if ( k >= 2)
                        {
                                if( k_heavy == 1)
                                  mols->ism_types[i] = 14;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 13;
                                else
                                  mols->ism_types[i] = 14;

                        }else if( k == 1){
                                mols->ism_types[i] = 15;


                        }else{
                                mols->ism_types[i] = 15;
                        }


		}else if ( mols->atoms[i] == 6 ) {  /* S */
	                if( k_heavy == 2)
	                        mols->ism_types[i] = 18;
			else if( k_heavy == 3 )
                                mols->ism_types[i] = 17;
                        else if( k_heavy == 4 )
                                mols->ism_types[i] = 16;
			else
                                mols->ism_types[i] = 16;
		}else if ( mols->atoms[i] == 1 ) {  /* C */

                        if ( k >= 4)
                        {
                                if( k_heavy == 1)
                                  mols->ism_types[i] = 22;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 21;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 20;
                                else if ( k_heavy >= 4)
                                  mols->ism_types[i] = 19;
				else
                                  mols->ism_types[i] = 24;
			}else if ( k == 3){
                               if( k_heavy == 1)
                                  mols->ism_types[i] = 25;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 24;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 23;
				else
                                  mols->ism_types[i] = 24;
			}else if ( k == 2){
                               if( k_heavy == 1)
                                  mols->ism_types[i] = 26;
                                else
                                  mols->ism_types[i] = 24;
			}else{
                               mols->ism_types[i] = 24;
			}

		}else if ( mols->atoms[i] == 10) {
                        mols->ism_types[i] = 27;
                }else if ( mols->atoms[i] == 8) {
                        mols->ism_types[i] = 28;
                }else if ( mols->atoms[i] == 7) {
                        mols->ism_types[i] = 30;
                }else if ( mols->atoms[i] ==102 ) {
                        mols->ism_types[i] = 31;
                }else if ( mols->atoms[i] == 100) {
                        mols->ism_types[i] = 29;
                }else{
			fprintf(stderr,"LCPO model has no parameters for atom no %i. Assuming Csp2. Atom type: %i. GAFF type: %i\n",i+1,mols->atoms[i],mols->gaff_types[i]);
                        mols->ism_types[i] = 24;
		}

	}

	*mymol = mols;
	return 0;
}
int ism_select_atoms_by_ligand_noconformer(MOL2 **mymol, MOL2 *lig)
{
        MOL2 *mol = NULL;
        int i = 0, j = 0;
        float max_grid[3], min_grid[3];
        mol = *mymol;

        max_grid[0] = max_grid[1] = max_grid[2] = -9999.9f;
        min_grid[0] = min_grid[1] = min_grid[2] = 9999.9f;

        for( j = 0; j < lig->n_atoms; j++)
        {
                if( lig->x[j] > max_grid[0])
                 max_grid[0] = lig->x[j];
                if( lig->y[j] > max_grid[1])
                 max_grid[1] = lig->y[j];
                if( lig->z[j] > max_grid[2])
                 max_grid[2] = lig->z[j];

                if( lig->x[j] < min_grid[0])
                 min_grid[0] = lig->x[j];
                if( lig->y[j] < min_grid[1])
                 min_grid[1] = lig->y[j];
                if( lig->z[j] < min_grid[2])
                 min_grid[2] = lig->z[j];
        }
        for( i = 0; i < 3 ; i++)
        {
                max_grid[i] += 10.0f;
                min_grid[i] -= 10.0f;
        }


        for( i = 0; i < mol->n_atoms; i++)
        {
		if( mol->vdw_selection[i] == 2)
		   continue;

                mol->ism_selection[i] = 0;
                mol->vdw_selection[i] = 0;

                if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
                {
                        mol->ism_selection[i] = 1;
                        mol->vdw_selection[i] = 1;

                }
        }


        *mymol = mol;
        return 0;

}


int ism_select_atoms_by_ligand(MOL2 **mymol, MOL2 *lig)
{
        MOL2 *mol = NULL;
        int i = 0, j = 0;
        float max_grid[3], min_grid[3];
        mol = *mymol;

        max_grid[0] = max_grid[1] = max_grid[2] = -9999.9f;
        min_grid[0] = min_grid[1] = min_grid[2] = 9999.9f;


        for( i = 0; i < lig->nconformers; i++)
        {
                set_default_conformer(&lig, i);
                for( j = 0; j < lig->n_atoms; j++)
                {
                        if( lig->x[j] > max_grid[0])
                         max_grid[0] = lig->x[j];
                        if( lig->y[j] > max_grid[1])
                         max_grid[1] = lig->y[j];
                        if( lig->z[j] > max_grid[2])
                         max_grid[2] = lig->z[j];

                        if( lig->x[j] < min_grid[0])
                         min_grid[0] = lig->x[j];
                        if( lig->y[j] < min_grid[1])
                         min_grid[1] = lig->y[j];
                        if( lig->z[j] < min_grid[2])
                         min_grid[2] = lig->z[j];
                }
        }

        for( i = 0; i < 3 ; i++)
        {
                max_grid[i] += 10.0f;
                min_grid[i] -= 10.0f;
        }


        for( i = 0; i < mol->n_atoms; i++)
        {
                mol->ism_selection[i] = 0;
                mol->vdw_selection[i] = 0;

                if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
                {
                        mol->ism_selection[i] = 1;
                        mol->vdw_selection[i] = 1;

                }
        }


        *mymol = mol;
        return 0;

}


/**
 *	@brief Calculate ligand-protein overlaps
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param lig Ligand MOL2 structure
 *	@param prot Protein MOL2 structure
 *	@param myburied_lig_complex Atom buried flag for ligand
 *	@param myburied_prot_complex Atom buried flag for protein
 *	@param mynumOverlapsLigComplex Vector of number of overlaps for ligand in complex
 *	@param mynumOverlapsProtComplex Vector of number of overlaps for protein in complex
 *	@param myOverlapsMatrixLig Overlaps indices for ligand in complex
 *	@param myOverlapsMatrixProt Overlaps indices for protein in complex
 *	@param myAijMatrixLig Matrix of overlap contributions for ligand in complex
 *	@param myAijMatrixProt Matrix of overlap contributions for protein in complex
 *	@return 0 on success
 */
int ism_getExternalOverlaps_forHbonds(MOL2 *lig, MOL2 *prot, int **myburied_lig_complex, int **myburied_prot_complex, int **mynumOverlapsLigComplex, int **mynumOverlapsProtComplex, int ***myOverlapsMatrixLig, int ***myOverlapsMatrixProt, float ***myAijMatrixLig, float ***myAijMatrixProt)
{

	int *buried_lig_complex = NULL, *buried_prot_complex = NULL;
        int *numOverlapsLigComplex = NULL, *numOverlapsProtComplex = NULL;
	int **OverlapsMatrixLig = NULL, **OverlapsMatrixProt;
	float **AijMatrixLig = NULL, **AijMatrixProt = NULL;

        float r1 = 0.0f, r2 = 0.0f;
        float deltaRadii = 0.0f;
        float dx = 0.0f, dy = 0.0f, dz = 0.0, distance = 0.0f;


	int i = 0, j = 0;

	buried_lig_complex = *myburied_lig_complex;
	buried_prot_complex = *myburied_prot_complex;
	numOverlapsLigComplex = *mynumOverlapsLigComplex;
	numOverlapsProtComplex = *mynumOverlapsProtComplex;

	OverlapsMatrixLig = *myOverlapsMatrixLig;
	OverlapsMatrixProt = *myOverlapsMatrixProt;
	AijMatrixLig = *myAijMatrixLig;
	AijMatrixProt = *myAijMatrixProt;

	for( i = 0; i < lig->n_atoms; i++)
	{
		buried_lig_complex[i] = 0;
		numOverlapsLigComplex[i] = 0;
		
	}

        for( i = 0; i < prot->n_atoms; i++)
        {
                buried_prot_complex[i] = 0;
                numOverlapsProtComplex[i] = 0;

        }

	for( i = 0; i < prot->n_atoms; i++)
	{
		if( prot->ism_selection[i] != 1)
			continue;


		if( !( prot->gaff_types[i] == O || prot->gaff_types[i] == OH || prot->gaff_types[i] == OW || ( prot->atoms[i] == 3 &&  has_hydrogens(prot,i) == 0)))
   			continue;
		
                if( prot->atoms[i] < 100)
                  r1 = ism_radii[ prot->atoms[i]-1];
                else
                  r1 = ism_radii_metals[ prot->atoms[i]-100];

                r1 += ISM_PROBE;


		for( j = 0; j < lig->n_atoms; j++)
		{

			if( !(lig->gaff_types[j] == HO || lig->gaff_types[j] == HW || lig->gaff_types[j] == HN))
				continue;

                        if( lig->atoms[j] < 100)
                          r2 = ism_radii[ lig->atoms[j]-1];
                        else
                          r2 = ism_radii_metals[ lig->atoms[j]-100];

                        r2 += ISM_PROBE;

                        dx = prot->x[i] - lig->x[j];
                        dy = prot->y[i] - lig->y[j];
                        dz = prot->z[i] - lig->z[j];
                        distance = (dx * dx) + (dy * dy) + (dz * dz);
                        distance = sqrt(distance);

			if( distance > 2.5)
			  continue;

			/* We have a ga modulator for hbond geomtry */

                        if( distance < (r1 + r2))
                        {
                                numOverlapsLigComplex[j] = numOverlapsLigComplex[j] + 1;
                                if( numOverlapsLigComplex[j] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                numOverlapsProtComplex[i] = numOverlapsProtComplex[i] + 1;
                                if( numOverlapsProtComplex[i] > ISM_MAX_OVERLAPS-1)
                                {
                                        fprintf(stderr,"Number of overlaps exceed the allocated memory\n");
                                        fflush(stderr);
                                        return -1;
                                }

                                OverlapsMatrixProt[numOverlapsProtComplex[i]-1][i] = j;
                                OverlapsMatrixLig[numOverlapsLigComplex[j]-1][j] = i;
                                AijMatrixProt[ numOverlapsProtComplex[i]-1][i] = 2.0f * ISM_PI * r1 *
                                  (r1 - (distance / 2.0f) - ((r1*r1) - (r2*r2)) /
                                  ( 2.0f*distance));


                                AijMatrixLig[ numOverlapsLigComplex[j]-1][j] = 2.0f * ISM_PI * r2 *
                                  (r2 - (distance / 2.0f) - ((r2*r2) - (r1*r1)) /
                                  ( 2.0f*distance));

                                deltaRadii = r1 - r2;


                                if( distance - fabs(deltaRadii) <= 0.0f)
                                {
                                        if( deltaRadii <= 0.0)
                                                buried_prot_complex[i] = 1;
                                        if( deltaRadii >= 0.0)
                                                buried_lig_complex[j] = 1;
                                }
			}
		}
	}

	return 0;

}


/**
 *      @brief Calculate SASA for the ligand in complex only hbonds
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @param mymol Protein MOL2 structure
 *      @param total_surface Total surface without overlaps
 *      @param buried Buried atoms flag in the free ligand
 *      @param buried_complex Buried atoms flag in the complex
 *      @param numOverlaps Overlaps count for the free ligand
 *      @param numOverlapsComplex Overlaps count for the ligand in complex
 *      @param numOverlapsB Overalps count for the protein
 *      @param numOverlapsComplexB Overlaps count for the protein in complex
 *      @param AijMatrix Overlaps SASA matrix for the free ligand
 *      @param AijMatrixComplex Overlaps SASA matrix for the ligand in complex
 *      @param AijMatrixB Overlaps SASA matrix for the protein
 *      @param AijMatrixComplexB Overlaps SASA matrix for the protein in complex
 *      @param OverlapsMatrix Overlaps indexes matrix for the free ligand
 *      @param OverlapsMatrixComplex Overlaps indexes matrix for the ligand in complex
 *      @param OverlapsMatrixB Overlaps indexes matrix for the protein
 *      @param OverlapsMatrixComplexB Overlaps indexes matrix for the protein in complex
 *      @param lcpoComps LCPO model components
 *      @param myatomSASA Resulting complex SASA
 *      @return 0 on success
 */
int ism_get_complex_sasa_hbond(MOL2 *mymol, float *total_surface, int *buried, int *buried_complex,
		      int *numOverlaps, int *numOverlapsComplex, int *numOverlapsB, int *numOverlapsComplexB,
		      float **AijMatrix, float **AijMatrixComplex, float **AijMatrixB, float **AijMatrixComplexB,
		      int **OverlapsMatrix, int **OverlapsMatrixComplex, int **OverlapsMatrixB, int **OverlapsMatrixComplexB,
		      float **lcpoComps, float **myatomSASA)
{


	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float *atomSASA = NULL;

	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{

		if( !(mymol->gaff_types[i] == HO || mymol->gaff_types[i] == HW || mymol->gaff_types[i] == HN))
			continue;

		tmpSASA = total_surface[i];
		if( buried[i] == 1 || buried_complex[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] + numOverlapsComplex[i] == 1){
				if( numOverlaps[i] == 1)
				 tmpSASA = tmpSASA - AijMatrix[0][i];
				else
				 tmpSASA = tmpSASA - AijMatrixComplex[0][i];
			}else if( (numOverlaps[i] + numOverlapsComplex[i]) > 1){
				lcpo[0] = lcpoComps[i][0];
				lcpo[1] = lcpoComps[i][1];
				lcpo[2] = lcpoComps[i][2];
				lcpo[3] = lcpoComps[i][3];

				for( j = 0; j < numOverlaps[i]; j++)
				{
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlapsComplex[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixComplex[k][idxOverlap];
						for( l = 0; l < numOverlapsComplex[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrixComplex[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrixComplex[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrixComplex[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				for( j = 0; j < numOverlapsComplex[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrixComplex[j][i];
                                        scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrixComplex[j][i];
                                        for( k = 0; k < numOverlapsB[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixB[k][idxOverlap];
                                                for( l = 0; l < numOverlapsComplex[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrixComplex[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

					}

                                        for( k = 0; k < numOverlapsComplexB[idxOverlap]; k++)
                                        {
						subidxOverlap = OverlapsMatrixComplexB[k][idxOverlap];
                                                for( l = 0; l < numOverlaps[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrix[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

                                        }

                                        lcpo[3] = lcpo[3] + (AijMatrixComplex[j][i] * scdlcpo4);

				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];
			}



		}

/*		fprintf(stderr,"%f/%f = %f\n",tmpSASA,total_surface[i],tmpSASA / total_surface[i]);*/

		if( tmpSASA > 0)
			tmpSASA = 1.0f - (tmpSASA / total_surface[i]);
		else
		 	tmpSASA = 0.0f;
	


		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;

	}



	*myatomSASA = atomSASA;

	return 0;

}


/**
 *	@brief Calculate SASA for the protein in complex
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Protein MOL2 structure
 *	@param total_surface Total surface without overlaps
 *	@param buried Buried atoms flag in the free protein
 *	@param buried_complex Buried atoms flag in the complex
 *	@param numOverlaps Overlaps count for the free protein
 *	@param numOverlapsComplex Overlaps count for the protein in complex
 *	@param numOverlapsB Overalps count for the ligand
 *	@param numOverlapsComplexB Overlaps count for the ligand in complex
 *	@param AijMatrix Overlaps SASA matrix for the free protein
 *	@param AijMatrixComplex Overlaps SASA matrix for the protein in complex
 *	@param AijMatrixB Overlaps SASA matrix for the ligand
 *	@param AijMatrixComplexB Overlaps SASA matrix for the ligand in complex
 *	@param OverlapsMatrix Overlaps indexes matrix for the free protein
 *	@param OverlapsMatrixComplex Overlaps indexes matrix for the protein in complex
 *	@param OverlapsMatrixB Overlaps indexes matrix for the ligand
 *	@param OverlapsMatrixComplexB Overlaps indexes matrix for the ligand in complex
 *	@param lcpoComps LCPO model components
 *	@param myatomSASA Resulting complex SASA
 *	@return 0 on success
 */
int ism_get_complex_sasa_hbond_prot(MOL2 *mymol, float *total_surface, int *buried, int *buried_complex,
		      int *numOverlaps, int *numOverlapsComplex, int *numOverlapsB, int *numOverlapsComplexB,
		      float **AijMatrix, float **AijMatrixComplex, float **AijMatrixB, float **AijMatrixComplexB,
		      int **OverlapsMatrix, int **OverlapsMatrixComplex, int **OverlapsMatrixB, int **OverlapsMatrixComplexB,
		      float **lcpoComps, float **myatomSASA)
{


	int i = 0, j = 0, k = 0, l = 0, idxOverlap = 0, subidxOverlap = 0;
	float tmpSASA = 0.0f;
	float totalSASA = 0.0f;
	float lcpo[4], scdlcpo4 = 0.0f;
	float *atomSASA = NULL;

	atomSASA = *myatomSASA;

	for( i = 0; i < mymol->n_atoms; i++)
	{

	  if( mymol->ism_selection != 0)
          {
		if( !(mymol->gaff_types[j] == HO || mymol->gaff_types[j] == HW || mymol->gaff_types[j] == HN))
			continue;

		tmpSASA = total_surface[i];
		if( buried[i] == 1 || buried_complex[i] == 1)
		{
			tmpSASA = 0.0f;
		}else{

			if( numOverlaps[i] + numOverlapsComplex[i] == 1){
				if( numOverlaps[i] == 1)
				 tmpSASA = tmpSASA - AijMatrix[0][i];
				else{
				 tmpSASA = tmpSASA - AijMatrixComplex[0][i];
				}

			}else if( (numOverlaps[i] + numOverlapsComplex[i]) > 1){
				lcpo[0] = lcpoComps[i][0];
				lcpo[1] = lcpoComps[i][1];
				lcpo[2] = lcpoComps[i][2];
				lcpo[3] = lcpoComps[i][3];

				for( j = 0; j < numOverlaps[i]; j++)
				{
					scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrix[j][i];
					for( k = 0; k < numOverlapsComplex[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixComplex[k][idxOverlap];
						for( l = 0; l < numOverlapsComplex[i]; l++)
						{
							if(subidxOverlap==OverlapsMatrixComplex[l][i])
							{
							  lcpo[2] = lcpo[2] + 
							      AijMatrixComplex[k][idxOverlap];
							  scdlcpo4 = scdlcpo4 + 
							      AijMatrixComplex[k][idxOverlap];
							  break;
							}
						}
					}

					lcpo[3] = lcpo[3] + AijMatrix[j][i] * scdlcpo4;
				}

				for( j = 0; j < numOverlapsComplex[i]; j++)
				{
					lcpo[1] = lcpo[1] + AijMatrixComplex[j][i];
                                        scdlcpo4 = 0.0f;
					idxOverlap = OverlapsMatrixComplex[j][i];
                                        for( k = 0; k < numOverlapsB[idxOverlap]; k++)
					{
						subidxOverlap = OverlapsMatrixB[k][idxOverlap];
                                                for( l = 0; l < numOverlapsComplex[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrixComplex[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

					}

                                        for( k = 0; k < numOverlapsComplexB[idxOverlap]; k++)
                                        {
						subidxOverlap = OverlapsMatrixComplexB[k][idxOverlap];
                                                for( l = 0; l < numOverlaps[i]; l++)
                                                {
                                                        if(subidxOverlap==OverlapsMatrix[l][i])
                                                        {
                                                          lcpo[2] = lcpo[2] +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          scdlcpo4 = scdlcpo4 +
                                                              AijMatrixComplexB[k][idxOverlap];
                                                          break;
                                                        }
                                                }

                                        }

                                        lcpo[3] = lcpo[3] + AijMatrixComplex[j][i] * scdlcpo4;

				}

				tmpSASA = lcpo[0] * LCPOParam1[ mymol->ism_types[i]-1] + 
					  lcpo[1] * LCPOParam2[ mymol->ism_types[i]-1] + 
					  lcpo[2] * LCPOParam3[ mymol->ism_types[i]-1] + 
					  lcpo[3] * LCPOParam4[ mymol->ism_types[i]-1] ;

/*				lcpoComps[i][0] = lcpo[0];
                                lcpoComps[i][1] = lcpo[1];
                                lcpoComps[i][2] = lcpo[2];
	                        lcpoComps[i][3] = lcpo[3];*/
			}



		}


		if( tmpSASA > 0)
			tmpSASA = 1.0f - (tmpSASA / total_surface[i]);
		else
		 	tmpSASA = 0.0f;
	

		atomSASA[i] = tmpSASA;
	        totalSASA += tmpSASA;
          }

	}


	*myatomSASA = atomSASA;
	return 0;

}


void calculate_ga(MOL2 *lig, MOL2 *prot, float **myga, float **myga_prot, int debug)
{
	float *ga = NULL;
	float *ga_prot = NULL;

        int i = 0, out_flag = 0, l = 0; 
        register int u1,  v1,  w1; 
        float e1, e2, e3, x, y, z, xd, yd, zd; 

	float dx = 0.0f, dy = 0.0f, dz = 0.0f;
	float dist = 0.0f, distr = 0.0f, alfa = 0.0f, beta = 0.0f, gbeta = 0.0f;
	float heavy_nb[3], vector1[3], vector2[3], modulus1 = 0.0f, modulus2 = 0.0f;
	int neig_j = 0, neig_index = 0;

	float lp = 0.0f, hb = 0.0f, glp = 0.0f;


	ga = *myga;
	ga_prot = *myga_prot;

        for ( l = 0; l < lig->n_atoms; ++l) 
	{
                x = lig->x[l];
                y = lig->y[l];
                z = lig->z[l];
              if( lig->gaff_types[l] == HO || lig->gaff_types[l] == HN || lig->gaff_types[l] == HW)
              {

		for( i = 0; i < prot->n_atoms; i++)
		{

                        	if( (prot->gaff_types[i] == O || prot->gaff_types[i] == OH || (prot->atoms[i] == 3 && has_hydrogens(prot,i) == 0 && prot->gaff_types[i] != N4)))
                                {
                                        dx = x - prot->x[i];
                                        dy = y - prot->y[i];
                                        dz = z - prot->z[i];
                                        dist = (dx*dx)+(dy*dy)+(dz*dz);
                                        dist = sqrt(dist);
			
                                	if (dist > 2.5f)
                                 	  continue;

                                 	heavy_nb[0] = lig->x[l];
                                 	heavy_nb[1] = lig->y[l];
                                 	heavy_nb[2] = lig->z[l];

                                        for(neig_j=0; neig_j < lig->n_bonds; ++neig_j)
                                        {
                                        	if( lig->bond_a1[neig_j] == (l+1))
                                           	{
                                                	neig_index = lig->bond_a2[neig_j]-1;
		                                        if ( lig->atoms[neig_index] == 3 || lig->atoms[neig_index] == 2)
		                                        {
	                                                	heavy_nb[0] = lig->x[neig_index];
	                                                	heavy_nb[1] = lig->y[neig_index];
	                                                	heavy_nb[2] = lig->z[neig_index];
	                                        	}
                                                }else if( lig->bond_a2[neig_j] == (l+1)){
                                                	neig_index = lig->bond_a1[neig_j]-1;
                                      			if ( lig->atoms[neig_index] == 3 || lig->atoms[neig_index] == 2)
                                      			{
                                               			heavy_nb[0] = lig->x[neig_index];
                                               			heavy_nb[1] = lig->y[neig_index];
                                               			heavy_nb[2] = lig->z[neig_index];


                                       			}

                                           	}
                                         }


                                        /* [N,O]-H */
                                        vector1[0] = heavy_nb[0] - x;
                                        vector1[1] = heavy_nb[1] - y;
                                        vector1[2] = heavy_nb[2] - z;

                                        modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                        /* H...[N,O] */
                                        vector2[0] = prot->x[i] - x;
                                        vector2[1] = prot->y[i] - y;
                                        vector2[2] = prot->z[i] - z;

                                        modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                        alfa = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                        alfa = 180.0f*alfa/3.14159f;


					/* beta  */
                                         heavy_nb[0] = prot->x[i];
                                         heavy_nb[1] = prot->y[i];
                                         heavy_nb[2] = prot->z[i];
					 beta = 0.0f;
					 gbeta = 1.0f;
                                         for(neig_j=0; neig_j < prot->n_bonds; ++neig_j)
                                         {
                                           if( prot->bond_a1[neig_j] == (i+1))
                                           {
                                               neig_index = prot->bond_a2[neig_j]-1;
                                               heavy_nb[0] = prot->x[neig_index];
                                               heavy_nb[1] = prot->y[neig_index];
                                               heavy_nb[2] = prot->z[neig_index];

                                                /* [N,O]...H */
                                                vector1[0] = x - prot->x[i];
                                                vector1[1] = y - prot->y[i];
                                                vector1[2] = z - prot->z[i];
                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                /* [N,O]-X */
                                                vector2[0] = heavy_nb[0] - prot->x[i];
                                                vector2[1] = heavy_nb[1] - prot->y[i];
                                                vector2[2] = heavy_nb[2] - prot->z[i];
                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
						beta = 180.0f*beta/3.14159f;
						gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                           }else if( prot->bond_a2[neig_j] == (i+1)){
                                               neig_index = prot->bond_a1[neig_j]-1;
                                               heavy_nb[0] = prot->x[neig_index];
                                               heavy_nb[1] = prot->y[neig_index];
                                               heavy_nb[2] = prot->z[neig_index];

                                                /* [N,O]...H */
                                                vector1[0] = x - prot->x[i];
                                                vector1[1] = y - prot->y[i];
                                                vector1[2] = z - prot->z[i];
                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                /* [N,O]-X */
                                                vector2[0] = heavy_nb[0] - prot->x[i];
                                                vector2[1] = heavy_nb[1] - prot->y[i];
                                                vector2[2] = heavy_nb[2] - prot->z[i];
                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                beta = 180.0f*beta/3.14159f;
                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                           }
                                         }

                                        e1 = 0.0f;
                                        e1 = (eval_block( (dist-1.85f), 0.25f, 0.65f));
                                        e1 *= eval_block( fabs(180.0f-alfa), 30.0f, 80.0f);
					e1 *= gbeta;
					if( debug )
					{
						fprintf(stderr,"Putative hbond: %f %f %f %f %f\n",e1,alfa,gbeta,dist,beta);
						fflush(stderr);	
					}

					if( e1 > 0.25)
					{
						if( lig->gaff_types[l] == HO || lig->gaff_types[l] == HW )
						{
							if( prot->gaff_types[i] == O)
							{
								if( is_carboxylate(prot, i) == 0)
									ga[l] = 0.6f;
								else
									ga[l] = 0.8f;
							}else if( prot->gaff_types[i] == OH || prot->gaff_types[i] == OW){
									ga[l] = 1.1f;
							}else if( prot->atoms[i] == 3){
									ga[l] = 0.7f;
							}else{
								fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for atom %i\n",l);
							}
						}else if( lig->gaff_types[l] == HN){

							if( is_hn_amido(lig,l) == 1)
							{
                                                               if( prot->gaff_types[i] == O)
                                                                {
                                                                        if( is_carboxylate(prot, i) == 0)
                                                                                ga[l] = 0.1f;
                                                                        else
                                                                                ga[l] = 0.38f;
                                                                }else if( prot->gaff_types[i] == OH || prot->gaff_types[i] == OW){
                                                                                ga[l] = 0.4f;
                                                                }else if( prot->atoms[i] == 3){
                                                                                ga[l] = 0.19f;
                                                                }else{
                                                                        fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for atom %i\n",l);
                                                                }

							}else{
	                                                        if( prot->gaff_types[i] == O)
	                                                        {
	                                                                if( is_carboxylate(prot, i) == 0)
	                                                                        ga[l] = 0.4f;
	                                                                else
	                                                                        ga[l] = 0.8f;
	                                                        }else if( prot->gaff_types[i] == OH || prot->gaff_types[i] == OW){
	                                                                        ga[l] = 1.0f;
	                                                        }else if( prot->atoms[i] == 3){
	                                                                        ga[l] = 0.5f;
	                                                        }else{
	                                                                fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for atom %i\n",l);
	                                                        }
							}

						}
					}else{
						ga[l] = 0.0f;
					}
					

/*                                                        fprintf(stderr,"Hydrogen bond D: %f %f %f %f\n",dist,alfa,gbeta,e1);
							fflush(stderr);*/
/*							fprintf(stderr,"End selection\n");
                                                        fflush(stderr);*/

				}
			}
		}else if( lig->gaff_types[l] == OH || lig->gaff_types[l] == OW || lig->gaff_types[l] == O || (lig->atoms[l] == 3 && has_hydrogens(lig,l) == 0) ){
			for( i = 0; i < prot->n_atoms; i++)
			{
                                              if( (prot->gaff_types[i] == HN || prot->gaff_types[i]  == HO || prot->gaff_types[i] == HW) )
                                              {
                                                        dx = x - prot->x[i];
                                                        dy = y - prot->y[i];
                                                        dz = z - prot->z[i];
                                                        dist = (dx*dx)+(dy*dy)+(dz*dz);
                                                        dist = sqrt(dist);

						if (dist > 2.5f)
						 continue;
						

                                                 heavy_nb[0] = prot->x[i];
                                                 heavy_nb[1] = prot->y[i];
                                                 heavy_nb[2] = prot->z[i];

                                                         for(neig_j=0; neig_j < prot->n_bonds; ++neig_j)
                                                         {
                                                           if( prot->bond_a1[neig_j] == (i+1))
                                                           {
                                                                neig_index = prot->bond_a2[neig_j]-1;
                                                      if ( prot->atoms[neig_index] == 3 || prot->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];

                                                      }
                                                           }else if( prot->bond_a2[neig_j] == (i+1)){
                                                               neig_index = prot->bond_a1[neig_j]-1;
                                                      if ( prot->atoms[neig_index] == 3 || prot->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];


                                                       }

                                                           }
                                                         }


                                                        vector1[0] = heavy_nb[0] - prot->x[i];
                                                        vector1[1] = heavy_nb[1] - prot->y[i];
                                                        vector1[2] = heavy_nb[2] - prot->z[i];

                                                        modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                        vector2[0] = dx;
                                                        vector2[1] = dy;
                                                        vector2[2] = dz;

                                                        modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                                        alfa = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
							alfa = 180.0f*alfa/3.14159f;



	                                                 heavy_nb[0] = lig->x[l];
	                                                 heavy_nb[1] = lig->y[l];
	                                                 heavy_nb[2] = lig->z[l];
							 beta = 0.0f;
							 gbeta = 1.0f;
                                                         for(neig_j=0; neig_j < lig->n_bonds; ++neig_j)
                                                         {
                                                           if( lig->bond_a1[neig_j] == (l+1))
                                                           {
                                                               neig_index = lig->bond_a2[neig_j]-1;
                                                               heavy_nb[0] = lig->x[neig_index];
                                                               heavy_nb[1] = lig->y[neig_index];
                                                               heavy_nb[2] = lig->z[neig_index];

                                                                vector1[0] = prot->x[i] - x;
                                                                vector1[1] = prot->y[i] - y;
                                                                vector1[2] = prot->z[i] - z;
                                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                                vector2[0] = heavy_nb[0] - x;
                                                                vector2[1] = heavy_nb[1] - y;
                                                                vector2[2] = heavy_nb[2] - z;
                                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                                beta = 180.0f*beta/3.14159f;
                                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);
                                                           }else if( lig->bond_a2[neig_j] == (l+1)){
                                                               neig_index = lig->bond_a1[neig_j]-1;
                                                               heavy_nb[0] = lig->x[neig_index];
                                                               heavy_nb[1] = lig->y[neig_index];
                                                               heavy_nb[2] = lig->z[neig_index];


                                                                vector1[0] = prot->x[i] - x;
                                                                vector1[1] = prot->y[i] - y;
                                                                vector1[2] = prot->z[i] - z;
                                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                                vector2[0] = heavy_nb[0] - x;
                                                                vector2[1] = heavy_nb[1] - y;
                                                                vector2[2] = heavy_nb[2] - z;
                                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                                beta = 180.0f*beta/3.14159f;
                                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                                           }
                                                         }


							e1 = 0.0f;
							e1 = (eval_block( (dist-1.85f), 0.25f, 0.65f));
							e1 *= eval_block( fabs(180.0f-alfa), 30.0f, 80.0f);
							e1 *= gbeta;

					if( debug )
					{
	                                        fprintf(stderr,"Putative prot-hbond: %f %f %f %f %f\n",e1,alfa,gbeta,dist,beta);
						fflush(stderr);
					}

                                        if( e1 > 0.25)
                                        {
                                                if( prot->gaff_types[i] == HO || prot->gaff_types[i] == HW )
                                                {
                                                        if( lig->gaff_types[l] == O)
                                                        {
                                                                if( is_carboxylate(lig, l) == 0)
                                                                        ga_prot[i] = 0.6f;
                                                                else
                                                                        ga_prot[i] = 0.8f;
                                                        }else if( lig->gaff_types[l] == OH || lig->gaff_types[l] == OW){
                                                                        ga_prot[i] = 1.1f;
                                                        }else if( lig->atoms[l] == 3){
                                                                        ga_prot[i] = 0.7f;
                                                        }else{
                                                                fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for prot atom %i\n",i);
                                                        }
                                                }else if( prot->gaff_types[i] == HN){

                                                        if( is_hn_amido(prot,i) == 1)
                                                        {
                                                               if( lig->gaff_types[l] == O)
                                                                {
                                                                        if( is_carboxylate(lig, l) == 0)
                                                                                ga_prot[i] = 0.1f;
                                                                        else
                                                                                ga_prot[i] = 0.38f;
                                                                }else if( lig->gaff_types[l] == OH || lig->gaff_types[l] == OW){
                                                                                ga_prot[i] = 0.4f;
                                                                }else if( lig->atoms[l] == 3){
                                                                                ga_prot[i] = 0.19f;
                                                                }else{
                                                                        fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for prot atom %i\n",i);
                                                                }

                                                        }else{
                                                                if( lig->gaff_types[l] == O)
                                                                {
                                                                        if( is_carboxylate(lig, l) == 0)
                                                                                ga_prot[i] = 0.4f;
                                                                        else
                                                                                ga_prot[i] = 0.8f;
                                                                }else if( lig->gaff_types[l] == OH || lig->gaff_types[l] == OW){
                                                                                ga_prot[i] = 1.0f;
                                                                }else if( lig->atoms[l] == 3){
                                                                                ga_prot[i] = 0.5f;
                                                                }else{
                                                                        fprintf(stderr,"Warning: SPC-ISM here. No hydrogen bond correction for prot atom %i\n",i);
                                                                }
                                                        }

                                                }
                                        }else{
                                                ga_prot[l] = 0.0f;
                                        }





							/* */


					       }

					}

				}

                }

	*myga = ga;
	*myga_prot = ga_prot;
        return;
}

