/**
 *
 *	@file main.c
 *	@brief Main file of the cMMISMSA scoring function tool.
 *
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 23/09/2015
 *
 * 	This program is free software; you can redistribute it and/or modify
 * 	it under the terms of the GNU General Public License as published by
 * 	the Free Software Foundation version 3 of the License.
 *
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 	GNU General Public License for more details.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <getopt.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

/*
 *  Rev 5 - February 2021
 *        - Fixed bugs in PDB atom typing and protein typing based on dictionary
 *        - Fixed bugs in PDB procesing (multichain systems)
 *
 *	Rev 4 - March 2018
 *        - Added full vdW/qq/solv residue decomposition in output for COMBINE-like methods
 *        - Added support for XTC GROMACS trajectories (still require AMBER topology)
 *        - Excluded atoms mask selects individual atoms and not full residues to allow
 *          better analysis of protein interfaces and Ala-scanning
 *        - LAPACK can be disabled on copilation time to enable a more portable code
 *
 *	Rev 3 - January 2016
 * 	  - MultiPDB support (experimental) for GROMACS exports
 *
 *	Rev 2 - November 2015
 *	  - Finished QH entropy approximation. Fixed bugs.
 *	  - Multiple files trajectories (experimental, limited naming scheme)
 *	  - Exclusion regions for protein regions analysis
 *	  - Fixed bugs regarding residue decomposition
 *
 *
 *	Rev 1 - January-April 2014
 *	 C reimplementation with fancy stuff:
 *	  - No dry topology needed. Waters/ions stripped on-the-fly
 *	  - "Advanced" atom masks to select ligand
 *	  - Automatic "box info" detection
 *	  - PDB with residue contributions in b-factor columns
 *	  - Much better memory management (no more SegFaults on < 2Gb boxes).
 *	    But OVERLAPS still limited to 300 
 *	  - Corrected Hydrogen bond model (no more ionic problems). 3 layers
 *	  - Simpler output (parseable CSV files)
 *	  - Debug info control
 *	  - Quasiharmonic entropy approximation thanks to the indexation of
 *	    the trajectory
 *	  - Faster/Slower depending on the case than Fortran MMISMSA :/
 *
 */

/** XDRFILE lib for XTC file support */
#include "xdrfile_xtc.h"
#include "xdrfile.c"
#include "xdrfile_xtc.c"
/* This is for reading molecules and other stuff */
#include <reader.c>
/* Memory clean up*/
#include <waste.c>
/* Implict Solvation Model support */
#include <ism.c>
/* I/O and energy function for AMBER files and force field */
#include <amber.c>
/* I/O extension to read XTC trajectories from GROMACS*/
#include <gromacs.c>
/* Atom masks routines */
#include <masks.c>
/* Protein trans-rot fitting routines */
#include <superimpose.c>
/* Diagonalization support throught LAPACKe C interface */
#ifdef _ENTROPY
#include <lapacke.h>
#endif

/* Definitions to avoid errors */
void show_some_pride ();
void show_help ();


/**
 *	@brief Main function
 *	@author Alvaro Cortes Cabrera
 *	@date 05/03/2013
 */
int
main (argc, argv)
     int argc;
     char *argv[];
{


  MOL2 *mol = NULL;
  TOP *topo = NULL;
  MOL2 *lig = NULL;
  MOL2 *prot = NULL;


  float MMenergy = 0, desolvEnergy = 0, totalEnergy = 0;

  /* Desolvation */
  /* Protein */
  int *numOverlaps_prot = NULL, *buried_prot = NULL, **OverlapsMatrix_prot =
    NULL;
  float **AijMatrix_prot = NULL;
  float **lcpoComps_prot = NULL, *total_surface_prot = NULL;
  float *atomSASA_prot = NULL, *atomSASAComplex_prot = NULL, r1 = 0.0f;
  float *atomSASAH_prot = NULL;
  /* Ligand */
  float **matrix = NULL;
  int *numOverlaps = NULL, *buried = NULL, **OverlapsMatrix = NULL;
  float **AijMatrix = NULL;
  float **lcpoComps = NULL, *total_surface = NULL;
  float *atomSASA = NULL, *atomSASAComplex = NULL;
  float *atomSASAH = NULL;
  /* Complex */
  int *buried_lig_complex = NULL, *buried_prot_complex = NULL;
  float **AijMatrixLig = NULL, **AijMatrixProt = NULL;
  int *numOverlapsLigComplex = NULL, *numOverlapsProtComplex = NULL;
  int **OverlapsMatrixLig = NULL, **OverlapsMatrixProt;
  /* H-bonds */
  float *ga = NULL, *ga_prot = NULL;
  /* Intra-molecular bonds */
  int *numOverlaps_hbond_prot = NULL, *buried_hbond_prot =
    NULL, **OverlapsMatrix_hbond_prot = NULL;
  float **AijMatrix_hbond_prot = NULL;
  int *numOverlaps_hbond = NULL, *buried_hbond =
    NULL, **OverlapsMatrix_hbond = NULL;
  float **AijMatrix_hbond = NULL;
  /* Inter-molecular bonds */
  int *buried_lig_complex_hbond = NULL, *buried_prot_complex_hbond = NULL;
  float **AijMatrixLig_hbond = NULL, **AijMatrixProt_hbond = NULL;
  int *numOverlapsLigComplex_hbond = NULL, *numOverlapsProtComplex_hbond =
    NULL;
  int **OverlapsMatrixLig_hbond = NULL, **OverlapsMatrixProt_hbond = NULL;
  /* General */
  float desolvR = 0.0f, desolvL = 0.0f, apolarcomplex = 0.0f;

  /* Residue decomposition */
  float *res_energies;
  float **res_energies_details = NULL;
  float **res_averages = NULL;
  int nres = 0, res_idx = 0, last_res = 0, max_res = 0;
  float current_qq = 0.0f, current_vdw = 0.0f;
  float qq_average = 0.0f, vdw_average = 0.0f, dl_average = 0.0f, dr_average =
    0.0f, apolar_average = 0.0f;
  float total_average = 0.0f;

  int i = 0, j = 0, k = 0, real_frames = 0, jj = 0, ii = 0, l = 0, ll =
    0, kk = 0;

  static int verbose_flag = 1;
  char c;

  int sframe = 0, eframe = -1, iframe = 1, reframe = -1;
  char *topology_file = NULL, *mdcrd_file = NULL, *rst_file = NULL;
  char *pdb_file = NULL;
  char *mask = NULL, *mask_exclude = NULL;
  char *output_prefix = NULL;
  int pdb_mode = 0, mdcrd_mode = 0, xtc_flag = 0, xtc_end = 0;

  FILE *f_energy_step = NULL, *f_energy_average = NULL;
  FILE *f_res_step = NULL, *f_res_average = NULL;
  FILE *f_pdb_byres = NULL;
  FILE *f_pdb_average = NULL;

  char *tmp_file_name = NULL;


  /* QH entropy approach */
  int entropy_flag = 1, ism_flag = 1, mm_flag = 1;
  float **reference_structure = NULL, xj = 0.0;
  float **average_structure = NULL;
  double *quasi_matrix = NULL, *quasi_vector = NULL, *quasi_values = NULL;
  int *entropy_selection = NULL;
  int entropy_natoms, entropy_idx = 0, entropy_idx2 = 0, entropy_frames = 0;
  float atom1_vec[3], atom2_vec[3], x1 = 0.0f;
#ifdef _ENTROPY
  lapack_int n = 0, m = 0;
  lapack_int il = 0, iu = 0;
  lapack_int *isuppz;
#endif
  double abstol = 0.0, vl = 0.0, vu = 0.0;
  double qh_entropy = 0.0;
  FILE *f_eigenvals = NULL;


  /* Multi-MDCRD support */
  int current_mdcrd = 0, n_mdcrds = 0;
  char *base_mdcrd = NULL;

  show_some_pride ();

  while (1)
    {
      static struct option long_options[] = {
	{"verbose", no_argument, &verbose_flag, 1},
	{"brief", no_argument, &verbose_flag, 0},
	{"topology", required_argument, 0, 't'},
	{"mdcrd", required_argument, 0, 'm'},
	{"xtc", required_argument, 0, 'x'},
	{"rst", required_argument, 0, 'r'},
	{"mask", required_argument, 0, 'n'},
	{"excluded", required_argument, 0, 'a'},
	{"sframe", required_argument, 0, 's'},
	{"eframe", required_argument, 0, 'e'},
	{"iframe", required_argument, 0, 'i'},
	{"output", required_argument, 0, 'o'},
	{"pdb", required_argument, 0, 'p'},
	{"multi", required_argument, 0, '4'},
	{"base", required_argument, 0, '5'},
	{"disable-ism", no_argument, 0, '2'},
	{"disable-mm", no_argument, 0, '6'},
	{"noentropy", no_argument, 0, '3'},
	{0, 0, 0, 0}
      };
      int option_index = 0;

      c = getopt_long (argc, argv, "p:l:t:m:r:s:e:i:x:23a",
		       long_options, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	case 0:
	  /* If this option set a flag, do nothing else now. */
	  if (long_options[option_index].flag != 0)
	    break;
	  printf ("option %s", long_options[option_index].name);
	  if (optarg)
	    printf (" with arg %s", optarg);
	  printf ("\n");
	  break;

	case '2':
	  ism_flag = 0;
	  break;

	case '6':
	  mm_flag = 0;
	  break;

	case '4':
	  n_mdcrds = atoi (optarg);
	  mdcrd_mode = 1;
	  break;

	case '5':
	  base_mdcrd = optarg;
	  mdcrd_mode = 1;
	  break;

	case '3':
	  entropy_flag = 0;
	  break;

	case 'p':
	  pdb_file = optarg;
	  pdb_mode = 2;
	  entropy_flag = 0;
	  break;

	case 't':
	  topology_file = optarg;
	  break;

	case 'o':
	  output_prefix = optarg;
	  break;
	case 'x':
	  xtc_flag = 1;
	  mdcrd_file = optarg;
	  mdcrd_mode = 1;
	  break;

	case 'm':
	  mdcrd_file = optarg;
	  mdcrd_mode = 1;
	  break;

	case 'r':
	  rst_file = optarg;
	  mdcrd_mode = 0;
	  break;

	case 'n':
	  mask = optarg;
	  break;

	case 'a':
	  mask_exclude = optarg;
	  break;

	case 's':
	  sframe = atoi (optarg);
	  break;

	case 'e':
	  eframe = atoi (optarg);
	  break;

	case 'i':
	  iframe = atoi (optarg);
	  break;

	case '?':
	  /* getopt_long already printed an error message. */
	  break;

	default:
	  abort ();
	}
    }


#ifndef _ENTROPY
  entropy_flag = 0;
  if (verbose_flag)
    fprintf (stderr,
	     "This binary was not compiled with entropy calculation support\n");
#endif


  if (xtc_flag != 0)
    entropy_flag = 0;

  /* Check parameters */

  if (rst_file == NULL && (mdcrd_file == NULL && base_mdcrd == NULL)
      && pdb_file == NULL)
    {
      fprintf (stderr,
	       "No valid input file. Please provide a AMBER topology/coordinates pair or a pdb file\n");
      fflush (stderr);
      show_help ();
      exit (-1);
    }

  if (topology_file == NULL && (rst_file != NULL || mdcrd_file != NULL))
    {
      fprintf (stderr,
	       "Topology provided but missing coordinates. Provide a valid mdcrd or rst/crd file\n");
      fflush (stderr);
      show_help ();
      exit (-1);
    }


  if (output_prefix == NULL)
    {
      fprintf (stderr, "No valid output prefix!\n");
      fflush (stderr);
      show_help ();
      exit (-1);
    }



  if (pdb_file != NULL)
    pdb_mode = 2;


  if (pdb_mode == 0)
    {

      /* Load TOP on MOL and TOPO */
      if (TOP_reader (&mol, topology_file, &topo, verbose_flag) < 0)
	{
	  fprintf (stderr, "Error loading topology file.\n");
	  fflush (stderr);
	  show_help ();
	  exit (-1);
	}

    }
  else
    {
      /*		PDB_reader(&mol,pdb_file,0);*//* Load and perciebe PDB */
      MultiPDB_reader (&mol, pdb_file, 0);	/* Load and perciebe PDB */
      fprintf (stderr, "PDB loaded\n");
      fflush (stderr);
    }

  /* Select everything for vdw+qq pairs */
  for (i = 0; i < mol->n_atoms; i++)
    {
      mol->vdw_selection[i] = 1;
    }

  /* Process the excluded mask selection. Like: "1-5,18,20-40,11,12,60-120" */
  if (mask_exclude != NULL)
    {
      if (verbose_flag)
	{
	  fprintf (stderr, "Exclusion ");
	}
      select_on_mask_atoms (&mol, mask_exclude, verbose_flag);
    }

  for (i = 0; i < mol->n_atoms; i++)
    {
      if (mol->backbone[i] == 1)
	{
	  mol->vdw_selection[i] = 2;
	}
    }
  /* Reset selection */
  for (i = 0; i < mol->n_atoms; i++)
    {
      mol->backbone[i] = 0;
    }

  if (verbose_flag)
    fprintf (stderr, "Ligand ");
  select_on_mask (&mol, mask, verbose_flag);


  if (pdb_mode == 0)
    {
      if (mdcrd_mode == 1)
	{
	  if (base_mdcrd != NULL && n_mdcrds > 0)
	    {

	      if (verbose_flag)
		{
		  fprintf (stderr, "Multi-MDCRD mode\n");
		  fflush (stderr);
		}
	      get_next_mdcrd_name (base_mdcrd, &mdcrd_file, current_mdcrd);
	      if (verbose_flag)
		{
		  fprintf (stderr, "First MDCRD file: %s\n", mdcrd_file);
		  fflush (stderr);
		}
	    }
	  /* Load mdcrd file info and make and index */
	  if (xtc_flag == 0)
	    {
	      if (load_mdcrd (&mol, mdcrd_file, &topo, 0, verbose_flag) != 0)
		{
		  fprintf (stderr,
			   "Fatal error loading information from coordinate file\n");
		  fflush (stderr);
		  exit (-1);
		}
	    }
	  else
	    {
	      if (load_xtc (&mol, mdcrd_file, &topo, verbose_flag) != 0)
		{
		  fprintf (stderr,
			   "Fatal error loading information from coordinate file\n");
		  fflush (stderr);
		  exit (-1);
		}
	    }

	  if (entropy_flag == 1)
	    {
	      entropy_natoms = 0;
	      entropy_selection = (int *) calloc (sizeof (int), mol->n_atoms);
	      for (i = 0; i < mol->n_atoms; i++)
		{
		  if (mol->atom_names[i][0] == 'C'
		      && mol->atom_names[i][1] == 'A'
		      && mol->atom_names[i][2] == ' '
		      && mol->internal_res_num[i] < 2000)
		    {
		      entropy_selection[i] = 1;
		      entropy_natoms++;
		    }
		}

	      average_structure =
		(float **) calloc (sizeof (float *), entropy_natoms);
	      reference_structure =
		(float **) calloc (sizeof (float *), entropy_natoms);
	      for (i = 0; i < entropy_natoms; i++)
		{
		  average_structure[i] = (float *) calloc (sizeof (float), 3);
		  reference_structure[i] =
		    (float *) calloc (sizeof (float), 3);
		}
	      entropy_idx = 0;
	      for (i = 0; i < mol->n_atoms; i++)
		{
		  if (entropy_selection[i] == 1)
		    {
		      reference_structure[entropy_idx][0] = mol->x[i];
		      reference_structure[entropy_idx][1] = mol->y[i];
		      reference_structure[entropy_idx][2] = mol->z[i];
		      entropy_idx++;
		    }
		}
	      quasi_matrix =
		(double *) calloc (sizeof (double),
				   3 * entropy_natoms * 3 * entropy_natoms);


	      if (verbose_flag)
		{
		  fprintf (stderr,
			   "Number of atoms involved in entropy calculation: %i\n",
			   entropy_natoms);
		  fflush (stderr);
		}
	    }

	}
      else
	{
	  /* Load crd/rst file */
	  load_crd (&mol, rst_file, topo);
	}
    }

  /* Split ligand and protein in two MOLs to use standard routines */
  split_amber_mol (&mol, &lig, &prot);

  /* Calculate last residue from protein */
  nres = 1;
  last_res = prot->internal_res_num[0];
  for (i = 0; i < prot->n_atoms; i++)
    {
      if (prot->internal_res_num[i] != last_res)
	++nres;

      last_res = prot->internal_res_num[i];
    }
  fprintf (stderr, "Residues in protein/last in protein: %i/%i\n", nres,
	   last_res);
  max_res = MAX (last_res, nres);
/*	nres = last_res;*/
  /* Initialize energy and average values */
/*	res_energies = (float *) calloc(sizeof(float),nres);*/
  /*res_energies_details = (float **) calloc(sizeof(float*),nres+1); */

  res_energies_details = (float **) calloc (sizeof (float *), max_res + 1);
  /*res_averages = (float **) calloc(sizeof(float*),nres+1); */
  res_averages = (float **) calloc (sizeof (float *), max_res + 1);
  /*for( i = 0; i <= nres; ++i) */
  for (i = 0; i <= max_res; ++i)
    {
      res_averages[i] = (float *) calloc (sizeof (float), 4);
      res_energies_details[i] = (float *) calloc (sizeof (float), 3);
    }


	/*************************
 	*	ISM BLOCK
	*************************/

  /* ISM atom typing */
  ism_typing (&lig);
  ism_typing (&prot);

  /* Select the entire ligand */
  for (i = 0; i < lig->n_atoms; i++)
    {
      lig->ism_selection[i] = 1;
    }

  /* Select protein near ligand */
  ism_select_atoms_by_ligand_noconformer (&prot, lig);


  /* Remove excluded atoms */
  for (i = 0; i < prot->n_atoms; i++)
    {
      if (prot->vdw_selection[i] == 2)
	{
	  prot->ism_selection[i] = 0;
	}
    }


  /* Initiallize vars for three layers */
  ga = (float *) calloc (sizeof (float), lig->n_atoms);
  ga_prot = (float *) calloc (sizeof (float), prot->n_atoms);

  buried = (int *) calloc (sizeof (int), lig->n_atoms);
  numOverlaps = (int *) calloc (sizeof (int), lig->n_atoms);
  AijMatrix = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrix = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  lcpoComps = (float **) calloc (sizeof (float *), lig->n_atoms);
  total_surface = (float *) calloc (sizeof (float), lig->n_atoms);
  atomSASA = (float *) calloc (sizeof (float), lig->n_atoms);
  atomSASAComplex = (float *) calloc (sizeof (float), lig->n_atoms);
  atomSASAH = (float *) calloc (sizeof (float), lig->n_atoms);

  matrix = (float **) calloc (sizeof (float *), lig->n_atoms);

  buried_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  numOverlaps_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  AijMatrix_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrix_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  for (i = 0; i < lig->n_atoms; i++)
    {
      lcpoComps[i] = (float *) calloc (sizeof (float), 4);
      matrix[i] = (float *) calloc (sizeof (float), lig->n_atoms);
    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      OverlapsMatrix[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      AijMatrix[i] = (float *) calloc (sizeof (float), lig->n_atoms);

      OverlapsMatrix_hbond[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      AijMatrix_hbond[i] = (float *) calloc (sizeof (float), lig->n_atoms);

    }

  /* Free SASA Hard spheres no overlapping model */
  for (i = 0; i < lig->n_atoms; i++)
    {
      if (lig->atoms[i] <= 100)
	r1 = ism_radii[lig->atoms[i] - 1];
      else
	r1 = ism_radii_metals[lig->atoms[i] - 1];

      r1 = r1 + ISM_PROBE;

      total_surface[i] = 4.0f * ISM_PI * (r1 * r1);
    }


  buried_lig_complex = (int *) calloc (sizeof (int), lig->n_atoms);
  buried_prot_complex = (int *) calloc (sizeof (int), prot->n_atoms);
  AijMatrixLig = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  AijMatrixProt = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  numOverlapsLigComplex = (int *) calloc (sizeof (int), lig->n_atoms);
  numOverlapsProtComplex = (int *) calloc (sizeof (int), prot->n_atoms);
  OverlapsMatrixLig = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrixProt = (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  buried_lig_complex_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  buried_prot_complex_hbond = (int *) calloc (sizeof (int), prot->n_atoms);
  AijMatrixLig_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  AijMatrixProt_hbond =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  numOverlapsLigComplex_hbond = (int *) calloc (sizeof (int), lig->n_atoms);
  numOverlapsProtComplex_hbond = (int *) calloc (sizeof (int), prot->n_atoms);
  OverlapsMatrixLig_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrixProt_hbond =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      OverlapsMatrixLig[i] = (int *) calloc (sizeof (int), lig->n_atoms);
      OverlapsMatrixProt[i] = (int *) calloc (sizeof (int), prot->n_atoms);
      AijMatrixLig[i] = (float *) calloc (sizeof (float), lig->n_atoms);
      AijMatrixProt[i] = (float *) calloc (sizeof (float), prot->n_atoms);
    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      OverlapsMatrixLig_hbond[i] =
	(int *) calloc (sizeof (int), lig->n_atoms);
      OverlapsMatrixProt_hbond[i] =
	(int *) calloc (sizeof (int), prot->n_atoms);
      AijMatrixLig_hbond[i] = (float *) calloc (sizeof (float), lig->n_atoms);
      AijMatrixProt_hbond[i] =
	(float *) calloc (sizeof (float), prot->n_atoms);
    }

  buried_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  numOverlaps_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  AijMatrix_prot = (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrix_prot =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);
  lcpoComps_prot = (float **) calloc (sizeof (float *), prot->n_atoms);
  total_surface_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  atomSASA_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  atomSASAComplex_prot = (float *) calloc (sizeof (float), prot->n_atoms);
  atomSASAH_prot = (float *) calloc (sizeof (float), prot->n_atoms);


  buried_hbond_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  numOverlaps_hbond_prot = (int *) calloc (sizeof (int), prot->n_atoms);
  AijMatrix_hbond_prot =
    (float **) calloc (sizeof (float *), ISM_MAX_OVERLAPS + 1);
  OverlapsMatrix_hbond_prot =
    (int **) calloc (sizeof (int *), ISM_MAX_OVERLAPS + 1);


  for (i = 0; i < prot->n_atoms; i++)
    {
      lcpoComps_prot[i] = (float *) calloc (sizeof (float), 4);

    }

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      OverlapsMatrix_prot[i] = (int *) calloc (sizeof (int), prot->n_atoms);
      AijMatrix_prot[i] = (float *) calloc (sizeof (float), prot->n_atoms);

      OverlapsMatrix_hbond_prot[i] =
	(int *) calloc (sizeof (int), prot->n_atoms);
      AijMatrix_hbond_prot[i] =
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

	  total_surface_prot[i] = 4.0f * ISM_PI * (r1 * r1);
	}
      else
	{
	  total_surface_prot[i] = 0.0f;
	}
    }

  /* End of ISM allocation */

	/********************************
  	*	END OF ISM SETUP
 	*********************************/

	/********************************
 	*	OUTPUT FILES
 	********************************/

  tmp_file_name =
    (char *) calloc (sizeof (char), strlen (output_prefix) + 40);

  if (mdcrd_mode == 1 || pdb_mode == 2)
    {
      sprintf (tmp_file_name, "%s_energy_averages", output_prefix);
      f_energy_average = fopen (tmp_file_name, "wb");
      sprintf (tmp_file_name, "%s_res_averages", output_prefix);
      f_res_average = fopen (tmp_file_name, "wb");
    }

  sprintf (tmp_file_name, "%s_energy", output_prefix);
  f_energy_step = fopen (tmp_file_name, "wb");
  sprintf (tmp_file_name, "%s_res", output_prefix);
  f_res_step = fopen (tmp_file_name, "wb");
  sprintf (tmp_file_name, "%s.pdb", output_prefix);
  f_pdb_byres = fopen (tmp_file_name, "wb");


  if (entropy_flag == 1)
    {
      sprintf (tmp_file_name, "%s_average.pdb", output_prefix);
      f_pdb_average = fopen (tmp_file_name, "wb");

      sprintf (tmp_file_name, "%s_eigenvals", output_prefix);
      f_eigenvals = fopen (tmp_file_name, "wb");
    }


	/**************************************
 	*	END OUTPUT FILES
 	**************************************/

  /* We do not cycle if there are only one set of coordinates */
  reframe = eframe;
  if ((pdb_mode == 1 || mdcrd_mode == 0) && pdb_mode != 2)
    {
      eframe = 1;
      sframe = 0;
      iframe = 1;
    }
  else if (pdb_mode == 2)
    {
      --sframe;
      eframe++;

      if (eframe <= 0)
	eframe = mol->nconformers - 1;

      fprintf (stderr, "Frames: %i\n", eframe);

      if (sframe < 0)
	sframe = 0;

      if (sframe >= mol->nconformers)
	sframe = 0;

    }
  else
    {
      /* Snapshot checks done */
      /* Humans count: 1-n snapshots, my loop goes 0 to (n-1) */
      --sframe;
      eframe++;

      if (eframe <= 0)
	eframe = topo->frames;

      if (sframe < 0)
	sframe = 0;

      if (sframe >= topo->frames)
	sframe = topo->frames - 1;
    }

  if (iframe <= 0)
    iframe = 1;			/* Reset step if it is wrong */

  fprintf (f_energy_step,
	   "Vdw;Colomb;DesolvLigand;DesolvReceptor;Apolar;Total\n");
  fflush (f_energy_step);

  /* We have no index and no number of frames in XTC files */
  /* A little hack is needed to move the loop till the end of file error in load_xtc */
  /* xtc_end flag is activated to let the loader know that it should iterate till that */
  /* point */
  if (xtc_flag != 0)
    {
      if (eframe == -1)
	{
	  eframe = 2;
	  xtc_end = 1;
	}
      else
	{
	  eframe -= 1;
	}

      /* At the moment, we do not support discarding frames or steps != 1 */
      sframe = 0;
      iframe = 1;
    }

  if (verbose_flag)
    {
      fprintf (stderr,
	       "Requested Starting frame: %i | End frame: %i | Step size: %i\n",
	       sframe, eframe, iframe);
      fprintf (stderr,
	       "--------------------------------------------------------\n");

      fflush (stderr);
    }

  for (i = sframe; i < eframe; i = i + iframe)
    {
      if (verbose_flag)
	fprintf (stderr, "Loading frame %i\n", i);

      if (mdcrd_mode == 1)
	{
	  if (xtc_flag == 0)
	    {
	      /* Load frame using index. Should be quick */
	      load_mdcrd (&mol, mdcrd_file, &topo, i, verbose_flag);
	    }
	  else
	    {
	      /* Load frame from XTC no indexing, no jumping */
	      if (load_xtc (&mol, mdcrd_file, &topo, verbose_flag) == 0)
		{
		  if (xtc_end == 1)
		    eframe++;
		}
	      else
		{		/* End frame */
		  eframe = i;
		  break;
		}
	    }
	  /* Update coordinates fro lig and prot structures */
	  split_amber_mol_crd (mol, &lig, &prot);

	  if (n_mdcrds > 0 && reframe <= 0)
	    eframe = topo->frames;	/* Reset max */
	}
      else if (pdb_mode == 2)
	{
	  set_default_conformer (&mol, i);
	  update_split_coordinates (mol, &lig, &prot);
	}

      /* ISM three-layers model */

      if (ism_flag == 1)
	{

	  ism_getInternalOverlaps_prot (prot, &numOverlaps_prot, &buried_prot,
					&AijMatrix_prot,
					&OverlapsMatrix_prot);
	  ism_get_free_sasa_prot (prot, total_surface_prot, buried_prot,
				  numOverlaps_prot, AijMatrix_prot,
				  OverlapsMatrix_prot, &lcpoComps_prot,
				  &atomSASA_prot);
	  ism_get_distance_matrix (lig, &matrix);
	  ism_getInternalOverlaps (lig, matrix, &numOverlaps, &buried,
				   &AijMatrix, &OverlapsMatrix);
	  ism_get_free_sasa (lig, total_surface, buried, numOverlaps,
			     AijMatrix, OverlapsMatrix, &lcpoComps,
			     &atomSASA);
	  ism_getExternalOverlaps (lig, prot, &buried_lig_complex,
				   &buried_prot_complex,
				   &numOverlapsLigComplex,
				   &numOverlapsProtComplex,
				   &OverlapsMatrixLig, &OverlapsMatrixProt,
				   &AijMatrixLig, &AijMatrixProt);
	  ism_get_complex_sasa (lig, total_surface, buried,
				buried_lig_complex, numOverlaps,
				numOverlapsLigComplex, numOverlaps_prot,
				numOverlapsProtComplex, AijMatrix,
				AijMatrixLig, AijMatrix_prot, AijMatrixProt,
				OverlapsMatrix, OverlapsMatrixLig,
				OverlapsMatrix_prot, OverlapsMatrixProt,
				lcpoComps, &atomSASAComplex);
	  ism_get_complex_sasa_prot (prot, total_surface_prot, buried_prot,
				     buried_prot_complex, numOverlaps_prot,
				     numOverlapsProtComplex, numOverlaps,
				     numOverlapsLigComplex, AijMatrix_prot,
				     AijMatrixProt, AijMatrix, AijMatrixLig,
				     OverlapsMatrix_prot, OverlapsMatrixProt,
				     OverlapsMatrix, OverlapsMatrixLig,
				     lcpoComps_prot, &atomSASAComplex_prot);
	  ism_getExternalOverlaps_forHbonds (lig, prot,
					     &buried_lig_complex_hbond,
					     &buried_prot_complex_hbond,
					     &numOverlapsLigComplex_hbond,
					     &numOverlapsProtComplex_hbond,
					     &OverlapsMatrixLig_hbond,
					     &OverlapsMatrixProt_hbond,
					     &AijMatrixLig_hbond,
					     &AijMatrixProt_hbond);
	  ism_get_complex_sasa_hbond (lig, total_surface, buried_hbond,
				      buried_lig_complex_hbond,
				      numOverlaps_hbond,
				      numOverlapsLigComplex_hbond,
				      numOverlaps_hbond_prot,
				      numOverlapsProtComplex_hbond,
				      AijMatrix_hbond, AijMatrixLig_hbond,
				      AijMatrix_hbond_prot,
				      AijMatrixProt_hbond,
				      OverlapsMatrix_hbond,
				      OverlapsMatrixLig_hbond,
				      OverlapsMatrix_hbond_prot,
				      OverlapsMatrixProt_hbond, lcpoComps,
				      &atomSASAH);
	  ism_getExternalOverlaps_forHbonds (prot, lig,
					     &buried_prot_complex_hbond,
					     &buried_lig_complex_hbond,
					     &numOverlapsProtComplex_hbond,
					     &numOverlapsLigComplex_hbond,
					     &OverlapsMatrixProt_hbond,
					     &OverlapsMatrixLig_hbond,
					     &AijMatrixProt_hbond,
					     &AijMatrixLig_hbond);
	  ism_get_complex_sasa_hbond (prot, total_surface_prot,
				      buried_hbond_prot,
				      buried_prot_complex_hbond,
				      numOverlaps_hbond_prot,
				      numOverlapsProtComplex_hbond,
				      numOverlaps_hbond,
				      numOverlapsLigComplex_hbond,
				      AijMatrix_hbond_prot,
				      AijMatrixProt_hbond, AijMatrix_hbond,
				      AijMatrixLig_hbond,
				      OverlapsMatrix_hbond_prot,
				      OverlapsMatrixProt_hbond,
				      OverlapsMatrix_hbond,
				      OverlapsMatrixLig_hbond, lcpoComps_prot,
				      &atomSASAH_prot);
	  calculate_ga (lig, prot, &ga, &ga_prot, verbose_flag);

	  desolvR =
	    ism_get_desolvation_energy (prot, lig, atomSASA_prot,
					atomSASAComplex_prot, atomSASAH_prot,
					ga_prot);
	  desolvL =
	    ism_get_desolvation_energy (lig, lig, atomSASA, atomSASAComplex,
					atomSASAH, ga);
	  apolarcomplex =
	    ism_get_apolar_energy (lig, prot, atomSASA, atomSASAComplex,
				   total_surface, atomSASA_prot,
				   atomSASAComplex_prot, total_surface_prot);


	  desolvEnergy = desolvL + desolvR + apolarcomplex;
	}
      /* End solvation */

      /* Get vdw and qq parts for all the selected atoms, which are all */

      if (mm_flag == 1)
	{
	  if (pdb_mode == 0)
	    MMenergy =
	      get_amber_vdw_selected (mol, topo, &res_energies_details,
				      &current_vdw, &current_qq);
	  else
	    MMenergy =
	      get_pdb_vdw_selected (mol, &res_energies_details, &current_vdw,
				    &current_qq);

	}

      if (verbose_flag)
	{
	  fprintf (stderr,
		   "vdW: %f | qq: %f | Desolvation (L/R/Apolar): %f %f %f\n",
		   current_vdw, current_qq, desolvL, desolvR, apolarcomplex);
	}

      /* Collect information for average structure */
      if (mdcrd_mode == 1 || pdb_mode == 2)
	{
	  if (entropy_flag == 1)
	    {
	      entropy_idx = 0;
	      fit_structure (&mol, reference_structure, entropy_selection,
			     entropy_natoms);
	      for (j = 0; j < mol->n_atoms; j++)
		{
		  if (entropy_selection[j] == 1)
		    {
		      average_structure[entropy_idx][0] += mol->x[j];
		      average_structure[entropy_idx][1] += mol->y[j];
		      average_structure[entropy_idx][2] += mol->z[j];
		      entropy_idx++;
		    }
		}
	      entropy_frames++;


	      if (verbose_flag)
		{
		  fprintf (stderr, "Adding structure %i to average\n", i);
		  fflush (stderr);
		}
	    }
	}

      /* Split desolvation by residue using MMres/MMall ratio */
      /* and do the averages also */
      for (j = 0; j < nres; j++)
	{
	  res_energies_details[j][2] +=
	    desolvEnergy *
	    ((res_energies_details[j][0] +
	      res_energies_details[j][1]) / MMenergy);
	  res_averages[j][0] +=
	    res_energies_details[j][0] + res_energies_details[j][1] +
	    res_energies_details[j][2];
	  res_averages[j][1] += res_energies_details[j][0];
	  res_averages[j][2] += res_energies_details[j][1];
	  res_averages[j][3] += res_energies_details[j][2];
	  fprintf (f_res_step, "%f;%f;%f;", res_energies_details[j][0],
		   res_energies_details[j][1], res_energies_details[j][2]);
	  res_energies_details[j][0] = 0.0f;
	  res_energies_details[j][1] = 0.0f;
	  res_energies_details[j][2] = 0.0f;
	}
      fprintf (f_res_step, "\n");

      totalEnergy = desolvEnergy + MMenergy;
      fprintf (f_energy_step, "%i,%f;%f;%f;%f;%f;%f\n", i + 1, current_vdw,
	       current_qq, desolvL, desolvR, apolarcomplex, totalEnergy);
      fflush (f_energy_step);

      qq_average += current_qq;
      vdw_average += current_vdw;
      dl_average += desolvL;
      dr_average += desolvR;
      apolar_average += apolarcomplex;
      total_average += totalEnergy;
      ++real_frames;


      if (n_mdcrds > 0 && ((i + iframe) >= eframe)
	  && current_mdcrd < n_mdcrds)
	{

	  current_mdcrd++;
	  free (mdcrd_file);
	  get_next_mdcrd_name (base_mdcrd, &mdcrd_file, current_mdcrd);
	  fprintf (stderr, "New MDCRD file: %s\n", mdcrd_file);
	  fflush (stderr);
	  reset_mdcrd (&topo);
	  i = sframe - iframe;
	}


    }

  if (xtc_flag != 0)
    xdrfile_close (topo->xdr);

  if (mdcrd_mode == 1 || pdb_mode == 2)
    {

      /* Finish the byres averages */
      fprintf (f_res_average, "ResidueNumber,TotalEnergy,Vdw,qq,ISM\n");
      fflush (f_res_average);

      for (j = 0; j < nres; j++)
	{
	  if (real_frames > 0)
	    {
	      res_averages[j][0] /= (float) real_frames;
	      res_averages[j][1] /= (float) real_frames;
	      res_averages[j][2] /= (float) real_frames;
	      res_averages[j][3] /= (float) real_frames;
	    }
	  else
	    {
	      res_averages[j][0] = 0.0f;
	      res_averages[j][1] = 0.0f;
	      res_averages[j][2] = 0.0f;
	      res_averages[j][3] = 0.0f;
	    }

	  fprintf (f_res_average, "%i,%f,%f,%f,%f\n", j + 1,
		   res_averages[j][0], res_averages[j][1], res_averages[j][2],
		   res_averages[j][3]);
	  fflush (f_res_average);
	}

      qq_average /= (float) real_frames;
      vdw_average /= (float) real_frames;
      dl_average /= (float) real_frames;
      dr_average /= (float) real_frames;
      apolar_average /= (float) real_frames;
      total_average /= (float) real_frames;

      if (entropy_flag == 0)
	{
	  fprintf (f_energy_average,
		   "VDW;QQ;DesolvLigand;DesolvReceptor;Apolar;Total\n");
	  fprintf (f_energy_average, "%f;%f;%f;%f;%f;%f\n", vdw_average,
		   qq_average, dl_average, dr_average, apolar_average,
		   total_average);
	}
    }

  last_res = -1;
  res_idx = 0;

  for (i = 0; i < mol->n_atoms; i++)
    {
      if (mol->backbone[i] == 1)
	{
	  mol->pcharges[i] = 0.0f;
	  continue;
	}


      if (pdb_mode == 2)
	{
	  if (last_res == -1)
	    {
	      last_res = mol->internal_res_num[i];
	      res_idx = 0;
	      mol->pcharges[i] = res_averages[res_idx][0];
	    }
	  else
	    {
	      if (last_res != mol->internal_res_num[i])
		{
		  ++res_idx;
		}
	      mol->pcharges[i] = res_averages[res_idx][0];
	      last_res = mol->internal_res_num[i];
	    }

	}
      else
	{
	  if (last_res == -1)
	    {
	      last_res = mol->internal_res_num[i];
/*			res_idx = 0;*/
	      res_idx = mol->internal_res_num[i] - 1;
	      mol->pcharges[i] = res_averages[res_idx][0];
	    }
	  else
	    {
	      if (last_res != mol->internal_res_num[i])
		{
/*			 ++res_idx;*/
		  res_idx = mol->internal_res_num[i] - 1;
		}
	      mol->pcharges[i] = res_averages[res_idx][0];
	      last_res = mol->internal_res_num[i];
	    }

	}
    }

  dump_pdb_conservative_file (mol, f_pdb_byres);
  for( j = 0; j < mol->n_atoms; j++)
  printf("%i %s %s %i %i %i %i %f\n", j, mol->atom_names[j], mol->res_names[j], mol->res_num[j], mol->ringer[j], mol->aromatic[j], mol->gaff_types[j], mol->pcharges[j] );

  if (mdcrd_mode == 1 || pdb_mode == 2)
    {
      if (entropy_flag == 1)
	{
	  for (j = 0; j < entropy_natoms; j++)
	    {
	      average_structure[j][0] /= (float) entropy_frames;
	      average_structure[j][1] /= (float) entropy_frames;
	      average_structure[j][2] /= (float) entropy_frames;
	    }
	  entropy_idx = 0;
	  for (j = 0; j < mol->n_atoms; j++)
	    {
	      if (entropy_selection[j] == 1)
		{
		  mol->x[j] = average_structure[entropy_idx][0];
		  mol->y[j] = average_structure[entropy_idx][1];
		  mol->z[j] = average_structure[entropy_idx][2];
		  entropy_idx++;
		}
	    }

	  fprintf (stderr, "Dumping average structure\n");
	  fflush (stderr);
	  dump_pdb_conservative_file (mol, f_pdb_average);


	  if (base_mdcrd != NULL && n_mdcrds > 0)
	    {
	      current_mdcrd = 0;
	      get_next_mdcrd_name (base_mdcrd, &mdcrd_file, current_mdcrd);
	      fprintf (stderr, "Current MDCRD file: %s\n", mdcrd_file);
	      fflush (stderr);
	      reset_mdcrd (&topo);

	      load_mdcrd (&mol, mdcrd_file, &topo, 0, verbose_flag);

	      if (reframe <= 0)
		eframe = topo->frames;	/* Reset max */
	    }

	  for (i = sframe; i < eframe; i = i + iframe)
	    {

	      if (verbose_flag)
		fprintf (stderr, "Loading frame %i\n", i);

	      if (pdb_mode == 2)
		{
		  set_default_conformer (&mol, i);
		  update_split_coordinates (mol, &lig, &prot);

		}
	      else
		{
		  load_mdcrd (&mol, mdcrd_file, &topo, i, verbose_flag);
		  split_amber_mol_crd (mol, &lig, &prot);
		}

	      if (n_mdcrds > 0 && reframe <= 0)
		eframe = topo->frames;	/* Reset max */


	      fit_structure (&mol, reference_structure, entropy_selection,
			     entropy_natoms);
	      entropy_idx = 0;
	      for (j = 0; j < mol->n_atoms; j++)
		{
		  if (entropy_selection[j] == 1)
		    {
		      mol->x[j] -= average_structure[entropy_idx][0];
		      mol->y[j] -= average_structure[entropy_idx][1];
		      mol->z[j] -= average_structure[entropy_idx][2];
		      entropy_idx++;
		    }
		}		/* Calculate delta */
	      entropy_idx = 0;
	      for (j = 0; j < mol->n_atoms; j++)
		{
		  if (entropy_selection[j] == 1)
		    {
		      atom1_vec[0] = mol->x[j];
		      atom1_vec[1] = mol->y[j];
		      atom1_vec[2] = mol->z[j];

		      for (ll = 0; ll < 3; ll++)
			{
			  k = 3 * entropy_natoms * (3 * entropy_idx + ll);
			  x1 = atom1_vec[ll];

			  entropy_idx2 = entropy_idx;	/* THIIIIIS WAS A BUUUUUUUG ( = 0 ) */
			  for (jj = j; jj < mol->n_atoms; jj++)
			    {
			      if (entropy_selection[jj] == 1)
				{
				  l = k + 3 * entropy_idx2;
				  atom2_vec[0] = mol->x[jj];
				  atom2_vec[1] = mol->y[jj];
				  atom2_vec[2] = mol->z[jj];
				  for (kk = 0; kk < 3; kk++)
				    {
/*									printf("%i,%i %i,%i %f*%f = %f\n",entropy_idx,entropy_idx2,ll,kk,atom2_vec[kk],x1,atom2_vec[kk]*x1);
									fflush(stdout);*/
				      quasi_matrix[l + kk] +=
					atom2_vec[kk] * x1;
				    }

				  entropy_idx2++;
				}
			    }
			}
		      entropy_idx++;
		    }
		}



	      if (n_mdcrds > 0 && ((i + iframe) >= eframe)
		  && current_mdcrd < n_mdcrds)
		{

		  current_mdcrd++;
		  free (mdcrd_file);
		  get_next_mdcrd_name (base_mdcrd, &mdcrd_file,
				       current_mdcrd);
		  fprintf (stderr, "New MDCRD file: %s\n", mdcrd_file);
		  fflush (stderr);
		  reset_mdcrd (&topo);
		  i = sframe - iframe;
		}


	    }			/* Loop over frames */
	}			/* Entropy ON check */
    }				/* Multi CRD check */


  if (mdcrd_mode == 1 || pdb_mode == 2)
    {
#ifdef _ENTROPY
      if (entropy_flag == 1)
	{

	  quasi_vector =
	    (double *) calloc (sizeof (double),
			       3 * entropy_natoms * 3 * entropy_natoms);
	  quasi_values =
	    (double *) calloc (sizeof (double), 3 * entropy_natoms);

	  for (i = 0; i < 3 * entropy_natoms; i++)
	    {
	      for (j = 0; j < 3 * entropy_natoms; j++)
		{
		  quasi_matrix[(i * 3 * entropy_natoms) +
			       j] /= (double) entropy_frames;
		}
	    }


	  for (j = 0; j < 3 * entropy_natoms; j++)
	    {
	      for (i = j; i < 3 * entropy_natoms; i++)
		{
		  quasi_matrix[3 * entropy_natoms * i + j] =
		    quasi_matrix[3 * entropy_natoms * j + i];
		}
	    }

	  if (verbose_flag)
	    {
	      fprintf (stderr, "Performing diagonalization ... ");
	      fflush (stderr);
	    }
	  n = 3 * entropy_natoms;

	  il = 1;
	  iu = n;
	  isuppz = (lapack_int *) calloc (sizeof (lapack_int), 2 * n);
	  vl = vu = 0;


	  LAPACKE_dsyevr (LAPACK_ROW_MAJOR, 'V', 'A', 'L', n, quasi_matrix, n,
			  0., 0., 0, 0, abstol, &m, quasi_values,
			  quasi_vector, n, isuppz);



	  if (verbose_flag)
	    {
	      fprintf (stderr, "done");
	      fflush (stderr);
	    }

	  for (i = 0; i < 3 * entropy_natoms; i++)
	    {
	      fprintf (f_eigenvals, "%i,%E\n", i, quasi_values[i]);
	    }

	  fflush (f_eigenvals);
	  fclose (f_eigenvals);

	  qh_entropy =
	    quasiharmonic_entropy (quasi_values, 3 * entropy_natoms);
	  fprintf (f_energy_average,
		   "VDW;QQ;DesolvLigand;DesolvReceptor;Apolar;Entropy;Total\n");
	  fprintf (f_energy_average, "%f;%f;%f;%f;%f;%f;%f\n", vdw_average,
		   qq_average, dl_average, dr_average, apolar_average,
		   qh_entropy, total_average + qh_entropy);

	}
#endif
    }

  /* Close output files */

  if (mdcrd_mode == 1 || pdb_mode == 2)
    {
      fclose (f_energy_average);
      fclose (f_res_average);


      if (entropy_flag == 1)
	{
	  fclose (f_eigenvals);
	  for (i = 0; i < entropy_natoms; i++)
	    {
	      free (average_structure[i]);
	      free (reference_structure[i]);
	    }
	  free (average_structure);
	  free (reference_structure);

	  free (entropy_selection);
	  free (quasi_matrix);
	  free (quasi_vector);
	  free (quasi_values);

	  fclose (f_pdb_average);
	}
    }
  fclose (f_energy_step);
  fclose (f_res_step);
  fclose (f_pdb_byres);
  free (tmp_file_name);


  free (ga);
  free (ga_prot);
  free (buried);
  free (numOverlaps);
  free (atomSASA);
  free (atomSASAComplex);
  free (atomSASAH);
  free (total_surface);
  free (buried_hbond);
  free (numOverlaps_hbond);

  for (i = 0; i < lig->n_atoms; i++)
    {
      free (lcpoComps[i]);
      free (matrix[i]);
    }
  free (lcpoComps);
  free (matrix);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      free (OverlapsMatrix[i]);
      free (AijMatrix[i]);
      free (OverlapsMatrix_hbond[i]);
      free (AijMatrix_hbond[i]);
    }
  free (AijMatrix_hbond);
  free (OverlapsMatrix_hbond);
  free (AijMatrix);
  free (OverlapsMatrix);



  free (buried_lig_complex);
  free (buried_prot_complex);
  free (numOverlapsLigComplex);
  free (numOverlapsProtComplex);
  free (buried_lig_complex_hbond);
  free (buried_prot_complex_hbond);
  free (numOverlapsLigComplex_hbond);
  free (numOverlapsProtComplex_hbond);

  for (i = 0; i < ISM_MAX_OVERLAPS; i++)
    {
      free (OverlapsMatrixLig[i]);
      free (OverlapsMatrixProt[i]);
      free (AijMatrixLig[i]);
      free (AijMatrixProt[i]);
      free (OverlapsMatrixLig_hbond[i]);
      free (OverlapsMatrixProt_hbond[i]);
      free (AijMatrixLig_hbond[i]);
      free (AijMatrixProt_hbond[i]);
      free (OverlapsMatrix_prot[i]);
      free (AijMatrix_hbond_prot[i]);
      free (OverlapsMatrix_hbond_prot[i]);
      free (AijMatrix_prot[i]);
    }
  free (AijMatrixLig);
  free (AijMatrixProt);
  free (OverlapsMatrixLig);
  free (OverlapsMatrixProt);
  free (OverlapsMatrixLig_hbond);
  free (OverlapsMatrixProt_hbond);
  free (AijMatrixLig_hbond);
  free (AijMatrixProt_hbond);
  free (OverlapsMatrix_prot);
  free (AijMatrix_hbond_prot);
  free (OverlapsMatrix_hbond_prot);
  free (AijMatrix_prot);
  free (buried_prot);
  free (numOverlaps_prot);
  free (total_surface_prot);
  free (atomSASA_prot);
  free (atomSASAComplex_prot);
  free (atomSASAH_prot);

  free (buried_hbond_prot);
  free (numOverlaps_hbond_prot);

  for (i = 0; i < prot->n_atoms; i++)
    {
      free (lcpoComps_prot[i]);

    }
  free (lcpoComps_prot);



  for (i = 0; i < nres + 1; i++)
    {
      free (res_energies_details[i]);
      free (res_averages[i]);
    }
  free (res_energies_details);
  free (res_averages);


  if (pdb_mode == 0)
    {
      clean_amber_mol (&mol);
      free (mol);
      clean_topology (&topo);
      free (topo);
    }
  else
    {
      cleanup (&mol);
    }

  /*clean_amber_split_mol (&lig);*/
  cleanup (&lig);
  /*clean_amber_split_mol (&prot);*/
  cleanup (&prot);

  /*free (lig);
  free (prot);*/


}

void
show_some_pride ()
{
  fprintf (stderr, "");
  fprintf (stderr,
	   "cMMISMSA - rev5 - Another program to quantify ligand-target interactions\n\n");
  fprintf (stderr, "Alvaro Cortes Cabrera <alvarocortesc@gmail.com>\n");
  fprintf (stderr, "Have fun!\n\n");
}

void
show_help ()
{
  fprintf (stderr, "List of options:\n\n");
  fprintf (stderr, "I/O:\n");
  fprintf (stderr, "--topology <AMBER topology> \n");
  fprintf (stderr, "--mdcrd <AMBER trajectory in ASCII MDCRD format> \n");
  fprintf (stderr, "--xtc <GROMACS trajectory in XTC> \n");
  fprintf (stderr, "--rst <AMBER coordinate file ASCII RST/CRD format> \n");
  fprintf (stderr, "--pdb <trajectory file in PDB format> \n");
  fprintf (stderr, "--output <prefix>\n\n");
  fprintf (stderr, "Misc:\n");
  fprintf (stderr,
	   "--mask <string> Residues to select as ligand. e.g. 1-6,10,20-30\n");
  fprintf (stderr, "--sframe <integer> Starting frame (from 1 to n)\n");
  fprintf (stderr, "--eframe <integer> Ending frame (from 1 to n)\n");
  fprintf (stderr, "--iframe <integer> Frame step\n");
  fprintf (stderr,
	   "--noentropy Disables Quasiharmonic approximation of entropy\n");
  fprintf (stderr, "--disable-ism Disables the Implicit Solvation Model\n");
  fprintf (stderr,
	   "--disable-mm Disables the AMBER vdw and electrostatic terms\n");
  fprintf (stderr, "--brief  Deactivate debugging information\n\n");


}
