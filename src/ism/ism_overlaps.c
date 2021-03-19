/**
 *
 *      @file ism.c
 *      @brief Implicit solvation model - overlap routines
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


