/**
 *
 *      @file ism.c
 *      @brief Implicit solvation model SASA routines
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @date 2021/03
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


