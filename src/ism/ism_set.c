/**
 *
 *      @file ism.c
 *      @brief Implicit solvation model - Setup routines
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

int ISM_init(ISM_COMPLEX **mycomplex, MOL2 *lig, MOL2 *prot)
{

	ISM_COMPLEX *cplex = NULL;
	int i = 0;
	float r1 = 0;

	if( (cplex = (ISM_COMPLEX *) calloc(sizeof(ISM_COMPLEX), 1)) == NULL)
	{
		fprintf(stderr,"Memory error.\n");
		return -1;
	}

  /* Initiallize vars for three layers */
  cplex->ga = (float *) calloc (sizeof (float), lig->n_atoms);
  cplex->ga_prot = (float *) calloc (sizeof (float), prot->n_atoms);

  cplex->buried = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->numOverlaps = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->AijMatrix = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrix = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  cplex->lcpoComps = (float **) calloc (sizeof (float *), lig->n_atoms);
  cplex->total_surface = (float *) calloc (sizeof (float), lig->n_atoms);
  cplex->atomSASA = (float *) calloc (sizeof (float), lig->n_atoms);
  cplex->atomSASAComplex = (float *) calloc (sizeof (float), lig->n_atoms);
  cplex->atomSASAH = (float *) calloc (sizeof (float), lig->n_atoms);

  cplex->matrix = (float **) calloc (sizeof (float *), lig->n_atoms);

  cplex->buried_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->numOverlaps_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->AijMatrix_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrix_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  for (i = 0; i < lig->n_atoms; i++)
    {
      cplex->lcpoComps[i] = (float *) calloc (sizeof (float), 4);
      cplex->matrix[i] = (float *) calloc (sizeof (float), lig->n_atoms);
    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      cplex->OverlapsMatrix[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      cplex->AijMatrix[i] = (float *) calloc (sizeof (float), lig->n_atoms);

      cplex->OverlapsMatrix_hbond[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      cplex->AijMatrix_hbond[i] = (float *) calloc (sizeof (float), lig->n_atoms);

    }

  /* Free SASA Hard spheres no overlapping model */
  for (i = 0; i < lig->n_atoms; i++)
    {
      if (lig->atoms[i] <= 100)
	r1 = ism_radii[lig->atoms[i] - 1];
      else
	r1 = ism_radii_metals[lig->atoms[i] - 1];

      r1 = r1 + ISM_PROBE;

      cplex->total_surface[i] = 4.0f * ISM_PI * (r1 * r1);
    }


  cplex->buried_lig_complex = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->buried_prot_complex = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->AijMatrixLig = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->AijMatrixProt = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->numOverlapsLigComplex = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->numOverlapsProtComplex = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->OverlapsMatrixLig = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrixProt = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  cplex->buried_lig_complex_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->buried_prot_complex_hbond = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->AijMatrixLig_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->AijMatrixProt_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->numOverlapsLigComplex_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  cplex->numOverlapsProtComplex_hbond = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->OverlapsMatrixLig_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrixProt_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      cplex->OverlapsMatrixLig[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      cplex->OverlapsMatrixProt[i] = (int *) calloc (sizeof (int), prot->n_atoms);
      cplex->AijMatrixLig[i] = (float *) calloc (sizeof (float), lig->n_atoms);
      cplex->AijMatrixProt[i] = (float *) calloc (sizeof (float), prot->n_atoms);
    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      cplex->OverlapsMatrixLig_hbond[i] =
	(int *) calloc (sizeof (int), lig->n_atoms);
      cplex->OverlapsMatrixProt_hbond[i] =
	(int *) calloc (sizeof (int), prot->n_atoms);
      cplex->AijMatrixLig_hbond[i] = (float *) calloc (sizeof (float), lig->n_atoms);
      cplex->AijMatrixProt_hbond[i] =
	(float *) calloc (sizeof (float), prot->n_atoms);
    }

  cplex->buried_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->numOverlaps_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->AijMatrix_prot = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrix_prot =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  cplex->lcpoComps_prot = (float **) calloc (sizeof (float *), prot->n_atoms);
  cplex->total_surface_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  cplex->atomSASA_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  cplex->atomSASAComplex_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  cplex->atomSASAH_prot = (float *) calloc (sizeof (float), prot->n_atoms);


  cplex->buried_hbond_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->numOverlaps_hbond_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  cplex->AijMatrix_hbond_prot =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  cplex->OverlapsMatrix_hbond_prot =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);


  for (i = 0; i < prot->n_atoms; i++)
    {
      cplex->lcpoComps_prot[i] = (float *) calloc (sizeof (float), 4);

    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      cplex->OverlapsMatrix_prot[i] = (int *) calloc (sizeof (int), prot->n_atoms);
      cplex->AijMatrix_prot[i] = (float *) calloc (sizeof (float), prot->n_atoms);

      cplex->OverlapsMatrix_hbond_prot[i] =
	(int *) calloc (sizeof (int), prot->n_atoms);
      cplex->AijMatrix_hbond_prot[i] =
	(float *) calloc (sizeof (float), prot->n_atoms);

    }

  for (i = 0; i < prot->n_atoms; i++)
    {
      if (prot->ism_selection[i] == 1)
	{
	  if (prot->atoms[i] <= 100)
	    r1 = ism_radii[prot->atoms[i] - 1];
	  else
	    r1 = ism_radii_metals[prot->atoms[i] - 1];

	  r1 = r1 + ISM_PROBE;

	  cplex->total_surface_prot[i] = 4.0f * ISM_PI * (r1 * r1);
	}
      else
	{
	  cplex->total_surface_prot[i] = 0.0f;
	}
    }

  *mycomplex = cplex;

  return 0;

}

void ISM_cleanup(ISM_COMPLEX **mycplex, MOL2 *lig, MOL2 *prot)
{
	ISM_COMPLEX *cplex = NULL;
	int i = 0;

	cplex = *mycplex;
  free (cplex->ga);
  free (cplex->ga_prot);
  free (cplex->buried);
  free (cplex->numOverlaps);
  free (cplex->atomSASA);
  free (cplex->atomSASAComplex);
  free (cplex->atomSASAH);
  free (cplex->total_surface);
  free (cplex->buried_hbond);
  free (cplex->numOverlaps_hbond);

  for (i = 0; i < lig->n_atoms; i++)
    {
      free (cplex->lcpoComps[i]);
      free (cplex->matrix[i]);
    }
  free (cplex->lcpoComps);
  free (cplex->matrix);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      free (cplex->OverlapsMatrix[i]);
      free (cplex->AijMatrix[i]);
      free (cplex->OverlapsMatrix_hbond[i]);
      free (cplex->AijMatrix_hbond[i]);
    }
  free (cplex->AijMatrix_hbond);
  free (cplex->OverlapsMatrix_hbond);
  free (cplex->AijMatrix);
  free (cplex->OverlapsMatrix);



  free (cplex->buried_lig_complex);
  free (cplex->buried_prot_complex);
  free (cplex->numOverlapsLigComplex);
  free (cplex->numOverlapsProtComplex);
  free (cplex->buried_lig_complex_hbond);
  free (cplex->buried_prot_complex_hbond);
  free (cplex->numOverlapsLigComplex_hbond);
  free (cplex->numOverlapsProtComplex_hbond);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      free (cplex->OverlapsMatrixLig[i]);
      free (cplex->OverlapsMatrixProt[i]);
      free (cplex->AijMatrixLig[i]);
      free (cplex->AijMatrixProt[i]);
      free (cplex->OverlapsMatrixLig_hbond[i]);
      free (cplex->OverlapsMatrixProt_hbond[i]);
      free (cplex->AijMatrixLig_hbond[i]);
      free (cplex->AijMatrixProt_hbond[i]);
      free (cplex->OverlapsMatrix_prot[i]);
      free (cplex->AijMatrix_hbond_prot[i]);
      free (cplex->OverlapsMatrix_hbond_prot[i]);
      free (cplex->AijMatrix_prot[i]);
    }
  free (cplex->AijMatrixLig);
  free (cplex->AijMatrixProt);
  free (cplex->OverlapsMatrixLig);
  free (cplex->OverlapsMatrixProt);
  free (cplex->OverlapsMatrixLig_hbond);
  free (cplex->OverlapsMatrixProt_hbond);
  free (cplex->AijMatrixLig_hbond);
  free (cplex->AijMatrixProt_hbond);
  free (cplex->OverlapsMatrix_prot);
  free (cplex->AijMatrix_hbond_prot);
  free (cplex->OverlapsMatrix_hbond_prot);
  free (cplex->AijMatrix_prot);
  free (cplex->buried_prot);
  free (cplex->numOverlaps_prot);
  free (cplex->total_surface_prot);
  free (cplex->atomSASA_prot);
  free (cplex->atomSASAComplex_prot);
  free (cplex->atomSASAH_prot);

  free (cplex->buried_hbond_prot);
  free (cplex->numOverlaps_hbond_prot);

  for (i = 0; i < prot->n_atoms; i++)
    {
      free (cplex->lcpoComps_prot[i]);

    }
  free (cplex->lcpoComps_prot);
  free(cplex);
}
