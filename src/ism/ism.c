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
#include <ism/ism.h>
#include <ism/ism_overlaps.c>
#include <ism/ism_typing.c>
#include <ism/ism_sasa.c>
#include <ism/ism_selection.c>
#include <ism/ism_hbond.c>


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

