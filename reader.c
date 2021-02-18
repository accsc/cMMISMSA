/**
 *
 *      @file reader.c
 *      @brief Read and allocation routine for PDB and MOL2 molecules
 *
 *        - Automatic parser: atom type, dihedrals, contacts, angles and bonds.
 *        - Import & export topology.
 *        - Modular energy terms.
 *        - Analytical gradients for most of the terms
 *        - Parser based on residue and atom name
 *        - Automatic interactions ligand-protein
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 01/10/2010
 *
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/* Molecular library definitions */
#include <libmol2.h>

/* GAFF force field definitions and constants */
#include <gaff.h>

#include <time.h>

/* General tools for reader */
#include "groups.c"
/* Ring detection routines */
#include "rings.c"

#define G_PI 3.14159265358979323846

/* Energetic calculation routines */
#include "conformers.c"
#include "energy.c"
#include "assign_bonds.c"
#include "assign_angles.c"
#include "assign_torsionals.c"
#include "assign_pairs.c"
#include "gaff_atom_types.c"
#include "residues.c"

#define MAX_BUFFER 1024


/**
 *
 *	@brief	Read and parse the contents of a PDB file 
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *
 *	@param mymol MOL2 structure to allocate and load
 *	@param finput_name Filename of the PDB
 *	@param import Import topology flag. 0. Do not import topology from files. 1. Import topology from top files
 *
 *      @return 0 if success
 *
 */
int
PDB_reader (MOL2 ** mymol, char *finput_name, int import)
{

  double elapsed;
  double treadab;
  double tring;
  double tdihe;
  double tiped;

  FILE *input;
  char *line = NULL;
  int molecules = 0;
  int i = 0;
  int i2 = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int step = 0;
  int cpair = 0;
  char tmp_atom[MAX_BUFFER];
  char tmp_bond[MAX_BUFFER];
  char tmp_radius[6];
  char tmp_charge[10];
  int rec_atoms = 0;
  MOL2 *mols = NULL;
  float pcharge = 0.0;
  int current_atom = 0;
  int current_bond = 0;
  int total_ppp = 0;
  char *adj = NULL;
  RING *test = NULL;
  char c = 0;
  char *input_name = NULL;
  int mode = 0;
  float x, y, z;
  int biggrad = 0;
  float lastforce = 0;
  float tmpforce = 0.0;
  float cbiggrad = 0.0f;
  float lastenergy = 999999.9f;
  float bestenergy = 0.0f;
  float tmpgrad = 0.0f;
  float dist = 0.0f;
  char name1[20];
  char name2[20];
  char name3[20];
  char name4[20];
  char name5[20];
  int ppp_count = 0;
  int my_type = 0;
  char myx[12];
  char myy[12];
  char myz[12];
  int bonds1, bonds2, bonds3;
  int flag = 0;
  int flag2 = 0;
  int flag3 = 0;
  int t1, t2, t3, t4;
  double bondE = 0;
  double angleE = 0;
  double torE = 0;
  double totalE = 0;
  double mytot;
  double vec1[3], vec2[3];
  double uvec1[3], uvec2[3];
  double *lvec1, *lvec2;
  double tor1;
  double tor2;
  double tor3;
  float thetha;
  int ttor = 0;
  int ewg = 0;
  int ewg2 = 0;
  int bondis = 0;

  RESIDUE *list_r = NULL;
  RESIDUE tmp_r;
  ROTAMER *tmp_rota = NULL;
  FILE *inres = NULL;
  FILE *rotlib = NULL;
  int tmp_resn_i = 0;
  char tmp_resn[10];
  int total_res = 0;
  int resn[1000];
  float target1[3], target2[3];
  float mx, my, mz;
  float brotE;
  int bestrot = -1;

  int cur_angle = 0;
  int *tor_type;

  float dx, dy, dz;
  float df;
  float ddf1;
  float df1;
  float e1;
  float dgx, dgy, dgz;
  float dtxi, dtxj, dtyi, dtyj, dtzi, dtzj;
  float dfx, dfy, dfz;
  float gaa, ra2, rb2, gbb;
  float fg, hg, fga, hgb;
  float r2, r;
  float RJR, RIR, cst;

  float *xc, *yc, *zc;
  float *xb, *yb, *zb;

  float idi, pn, pk, phi0;

  double A, B;
  double epsilon;
  double tmpE;
  double vdwE;
  float r6, r12;

  float rx = 0.0f, ry = 0.0f, rz = 0.0f;
  int matchs[3], nmatchs = 0;

  int conformers = 0;
  int current_conformer = 0;
  int no_count = 0;

  /* FLAGS */
  int gradflag;
  int verboseflag;

  char res_num[10];
  char res_type[10];

  int wats = 0;


/* CONTROL FLAGS */
  gradflag = 1;
  verboseflag = 0;



  if ((input = fopen (finput_name, "r")) == NULL)
    {
      fprintf (stderr, "Error. Cant open file %s.\n", finput_name);
      fflush (stderr);
      exit (-1);
    }



  if ((line = (char *) calloc (1, MAX_BUFFER + 2)) == NULL)
    {
      fprintf (stderr, "Error. Cant al memory.\n");
      fflush (stderr);
      exit (-2);
    }


  conformers = 1;

  while (fgets (line, MAX_BUFFER, input))
    {
      if ((strstr (line, "ATOM") != NULL || strstr (line, "HETATM") != NULL)
	  && no_count == 0)
	++molecules;

/*		if ( (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {
			conformers++;
			no_count = 1;
		}*/

/*                if ( (line[0] == 'T' && line[1] == 'E' && line[2] == 'R') || (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {
                        conformers++;
                        no_count = 1;
                }*/


    }

  rewind (input);

  if ((mols = (MOL2 *) calloc (sizeof (MOL2), 1)) == NULL)
    {
      fprintf (stderr, "Error. Cannot allocate big memory chunk.\n");
      fflush (stderr);
      exit (-2);
    }

  *mymol = mols;

  i = 0;
  mols->n_atoms = molecules;
  mols->x = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->y = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->z = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->pcharges = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->radius = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->atoms = (int *) calloc (sizeof (int), mols[i].n_atoms);
  mols->ringer = (int *) calloc (sizeof (int), mols[i].n_atoms);
  mols->aromatic = (int *) calloc (sizeof (int), mols[i].n_atoms);
  mols->bond_dist = (float *) calloc (sizeof (float), mols[i].n_atoms * 8);
  mols->grads_X = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->grads_Y = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->grads_Z = (float *) calloc (sizeof (float), mols[i].n_atoms);
  mols->backbone = (int *) calloc (sizeof (int), mols[i].n_atoms);
  mols->exclude = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->selection = (int *) calloc (sizeof (int), mols->n_atoms + 10);
  mols->fragment_flag = (int *) calloc (sizeof (int), mols[i].n_atoms);
  mols->conformers = (CONFORMER *) calloc (sizeof (CONFORMER), conformers);
  mols->n_fragments = 0;

  mols->res_num = (int *) calloc (sizeof (int *), mols->n_atoms + 10);
  mols->res_names = (char **) calloc (sizeof (char **), mols->n_atoms + 1);
  mols->atom_names = (char **) calloc (sizeof (char **), mols->n_atoms + 1);
  for (i = 0; i < mols->n_atoms; i++)
    {
      mols->res_names[i] = (char *) calloc (sizeof (char *), 4);
      mols->atom_names[i] = (char *) calloc (sizeof (char *), 4);
    }
  i = 0;
  mols->res_type = (int *) calloc (sizeof (int *), mols->n_atoms + 10);


  for (i = 0; i < conformers; ++i)
    {
      mols->conformers[i].x =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].y =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].z =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].pcharges =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].radius =
	(float *) calloc (sizeof (float), mols->n_atoms);
    }


  i = 0;

  mols->gaff_types = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);
  mols->ism_types = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);
  mols->hyde_types = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);

  mols->ism_selection = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);
  mols->hyde_selection = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);

  mols->vdw_selection = (int *) calloc (sizeof (int), mols[i].n_atoms + 1);


  if (verboseflag == 1)
    printf ("Cleaning up the gradients.\n");

  if (gradflag == 1)
    {
      for (j = 0; j < mols->n_atoms; ++j)
	{
	  mols->grads_X[j] = 0.0;
	  mols->grads_Y[j] = 0.0;
	  mols->grads_Z[j] = 0.0;
	}
    }
  j = 0;
  current_atom = 0;
  current_conformer = 0;
  mols->nconformers = conformers;
  mols->default_conformer = 0;

  while (fgets (line, MAX_BUFFER, input))
    {
#ifdef DEBUG
      fprintf (stderr, "Estado: Linea:%s", line);
      fflush (stderr);
#endif

/*		if ( (line[0] == 'T' && line[1] == 'E' && line[2] == 'R') || (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {*/
      if ((line[0] == 'E' && line[1] == 'N' && line[2] == 'D'
	   && line[3] != 'M'))
	{
/*			current_conformer++;
			current_atom = 0;*/
	}
      else if (strstr (line, "ATOM") != NULL
	       || strstr (line, "HETATM") != NULL)
	{
	  sscanf (line, "%*s %*d %4s", tmp_atom);
/*			strncpy(tmp_charge, &line[60], 8);
			strncpy(tmp_radius, &line[56], 4);*/

	  strncpy (tmp_radius, &line[62], 4);
	  strncpy (tmp_charge, &line[55], 7);

	  if (current_conformer == 0)
	    {
	      mols->pcharges[current_atom] = atof (tmp_charge);
	      mols->radius[current_atom] = atof (tmp_radius);

	      if (tmp_atom[0] == 'C'
		  && (tmp_atom[1] != 'l' && tmp_atom[1] != 'L'))
		mols->atoms[current_atom] = 1;
	      else if (tmp_atom[0] == 'O')
		mols->atoms[current_atom] = 2;
	      else if (tmp_atom[0] == 'N')
		mols->atoms[current_atom] = 3;
	      else if (tmp_atom[0] == 'H')
		mols->atoms[current_atom] = 4;
	      else if (tmp_atom[0] == 'P')
		mols->atoms[current_atom] = 5;
	      else if (tmp_atom[0] == 'S')
		mols->atoms[current_atom] = 6;
	      else if (tmp_atom[0] == 'I')
		mols->atoms[current_atom] = 7;
	      else if (tmp_atom[0] == 'B'
		       && (tmp_atom[1] == 'r' || tmp_atom[1] == 'R'))
		mols->atoms[current_atom] = 8;
	      else if (tmp_atom[0] == 'C'
		       && (tmp_atom[1] == 'l' || tmp_atom[1] == 'L'))
		mols->atoms[current_atom] = 9;
	      else if (tmp_atom[0] == 'F' && tmp_atom[1] != 'E'
		       && tmp_atom[1] != 'e')
		mols->atoms[current_atom] = 10;
	      else if (tmp_atom[0] == 'Z'
		       && (tmp_atom[1] == 'N' || tmp_atom[1] == 'n'))
		mols->atoms[current_atom] = 100;
	      else if (tmp_atom[0] == 'M'
		       && (tmp_atom[1] == 'N' || tmp_atom[1] == 'n'))
		mols->atoms[current_atom] = 101;
	      else if (tmp_atom[0] == 'M'
		       && (tmp_atom[1] == 'G' || tmp_atom[1] == 'g'))
		mols->atoms[current_atom] = 102;
	      else if (tmp_atom[0] == 'C' && tmp_atom[1] == '0')
		mols->atoms[current_atom] = 105;
	      else if (tmp_atom[0] == 'F'
		       && (tmp_atom[1] == 'E' || tmp_atom[1] == 'e'))
		mols->atoms[current_atom] = 106;
	      else if (tmp_atom[0] == 'K')
		mols->atoms[current_atom] = 103;
	      else if (tmp_atom[0] == 'N'
		       && (tmp_atom[1] == 'A' || tmp_atom[1] == 'a'))
		mols->atoms[current_atom] = 104;
	      else
		{

		  if (tmp_atom[1] == 'H')
		    mols->atoms[current_atom] = 4;
		  else
		    {
		      mols->atoms[current_atom] = 11;
		      mols->gaff_types[current_atom] = C3;
		      fprintf (stderr,
			       "Reader %i: Error. Do not know the atom type $%s$\n",
			       __LINE__, tmp_atom);
		    }
		}

	      strncpy (res_num, &line[23], 6);
	      mols->res_num[current_atom] = atoi (res_num);

	      strncpy (res_type, &line[17], 3);
	      res_type[3] = '\0';
	      mols->res_type[current_atom] =
		name_to_number_residue (res_type);

	      if (res_type[0] == 'W' && res_type[1] == 'A'
		  && res_type[2] == 'T')
		++wats;

	      strncpy (mols->res_names[current_atom], res_type, 3);
	      strncpy (mols->atom_names[current_atom], &line[12], 4);

	      if (mols->atoms[current_atom] >= 100)
		{
		  fprintf (stderr, "%s is metal",
			   mols->atom_names[current_atom]);
		}



	    }

	  strncpy (myx, &line[29], 10);
	  strncpy (myy, &line[38], 10);
	  strncpy (myz, &line[46], 10);

	  if (current_conformer == 0)
	    {
	      mols->x[current_atom] = atof (myx);
	      mols->y[current_atom] = atof (myy);
	      mols->z[current_atom] = atof (myz);
	    }

	  mols->conformers[current_conformer].x[current_atom] = atof (myx);
	  mols->conformers[current_conformer].y[current_atom] = atof (myy);
	  mols->conformers[current_conformer].z[current_atom] = atof (myz);
	  mols->conformers[current_conformer].pcharges[current_atom] =
	    atof (tmp_charge);
	  mols->conformers[current_conformer].radius[current_atom] =
	    atof (tmp_radius);

	  ++current_atom;


	}
    }


  free (line);
  fclose (input);

  if (mols->n_atoms > 20000)
    {
      fprintf (stderr,
	       "Warning: The system has more than 20000 atoms (%i) and more than 50 water molecules (%i).\nAutomatic topology detection will take ages or may fail to detect the correct bonds. Please remove all not needed water molecules of the system.\n",
	       mols->n_atoms, wats);
      fflush (stderr);
    }

  mol_percieve (&mols);


  return 0;
}

/**
 *
 *	@brief Percieve connectivity of the molecule (bonds, angles, torsionals, rings, aromaticity, etc.)
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *
 *	@param mymol MOL2 structure which contains the molecule
 *	@param import flag to import (1) the topology or not (0)
 *	@return 0 on success
 *
 */
int
mol_percieve (MOL2 ** mymol)
{

  MOL2 *mols = NULL;
  int *ringer = NULL, *aro = NULL, rnum = 0;
  int angleflag = 0, torflag = 0, vdwflag = 0;

  mols = *mymol;

  assign_bond_types (&mols);

  ringer = mols->ringer;
  aro = mols->aromatic;

  if (mols->n_atoms <= 300)
    get_number_of_rings2 (mols[0], &ringer, &aro);
  else{
      for( rnum = 0; rnum < get_number_of_residues(mols); rnum++)
         process_rings_residue(&mols, rnum);
  }

  GAFF_atom_typing (&mols);
  assign_gaff_vdw_types(&mols);

  if (angleflag)
    assign_angles (&mols);

  if (torflag)
    assign_torsionals (&mols);

  if (vdwflag)
    assign_pairs (&mols);

  *mymol = mols;

}

/**
 *
 *	@brief Initialize memory structures for a MOLECULE of n atoms and y conformers
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *
 *	@param mymol MOL2 structure which contains the molecule
 *	@param natoms Number of atoms
 *	@param conformers Number of conformers
 *	@return NULL
 *
 */
void
init_molecule (MOL2 ** mymol, int natoms, int conformers)
{
  MOL2 *mols = NULL;
  int i = 0;
  mols = *mymol;

  mols->n_atoms = natoms;
  mols->x = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->y = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->z = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->pcharges = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->radius = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->atoms = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->ringer = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->aromatic = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->bond_dist = (float *) calloc (sizeof (float), mols->n_atoms * 8);
  mols->grads_X = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->grads_Y = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->grads_Z = (float *) calloc (sizeof (float), mols->n_atoms);
  mols->backbone = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->selection = (int *) calloc (sizeof (int), mols->n_atoms + 10);
  mols->exclude = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->fragment_flag = (int *) calloc (sizeof (int), mols->n_atoms);
  mols->conformers = (CONFORMER *) calloc (sizeof (CONFORMER), conformers);
  mols->n_fragments = 0;

  mols->res_num = (int *) calloc (sizeof (int *), mols->n_atoms + 10);
  mols->internal_res_num = (int *) calloc (sizeof (int *), mols->n_atoms + 10);
  mols->res_names = (char **) calloc (sizeof (char **), mols->n_atoms + 1);
  mols->atom_names = (char **) calloc (sizeof (char **), mols->n_atoms + 1);
  for (i = 0; i < mols->n_atoms; i++)
    {
      mols->res_names[i] = (char *) calloc (sizeof (char *), 4);
      mols->atom_names[i] = (char *) calloc (sizeof (char *), 4);
    }

  mols->res_type = (int *) calloc (sizeof (int *), mols->n_atoms + 10);


  for (i = 0; i < conformers; ++i)
    {
      mols->conformers[i].x =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].y =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].z =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].pcharges =
	(float *) calloc (sizeof (float), mols->n_atoms);
      mols->conformers[i].radius =
	(float *) calloc (sizeof (float), mols->n_atoms);
    }

  mols->nconformers = conformers;
  mols->default_conformer = 0;

  mols->gaff_types = (int *) calloc (sizeof (int), mols->n_atoms + 1);
  mols->vdw_parm1 = (float *) calloc (sizeof (float), mols->n_atoms + 1);
  mols->vdw_parm2 = (float *) calloc (sizeof (float), mols->n_atoms + 1);
  mols->ism_types = (int *) calloc (sizeof (int), mols->n_atoms + 1);
  mols->hyde_types = (int *) calloc (sizeof (int), mols->n_atoms + 1);

  mols->ism_selection = (int *) calloc (sizeof (int), mols->n_atoms + 1);
  mols->hyde_selection = (int *) calloc (sizeof (int), mols->n_atoms + 1);

  mols->vdw_selection = (int *) calloc (sizeof (int), mols->n_atoms + 1);

  *mymol = mols;

}




/**
 *
 *	@brief	Read and parse the contents of a MultiPDB file 
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *
 *	@param mymol MOL2 structure to allocate and load
 *	@param finput_name Filename of the PDB
 *	@param import Import topology flag. 0. Do not import topology from files. 1. Import topology from top files
 *
 *      @return 0 if success
 *
 */
int
MultiPDB_reader (MOL2 ** mymol, char *finput_name, int import)
{

  FILE *input = NULL;
  MOL2 *mols = NULL;
  char *line = NULL;
  char tmp_atom[MAX_BUFFER], tmp_radius[6], tmp_charge[10];
  char res_num[10], res_type[10];
  char myx[12], myy[12], myz[12];

  int molecules = 0, i = 0, j = 0, wats = 0, prev_res_num = -999, internal_res_num = -1;
  int current_atom = 0, current_conformer = 0, no_count = 0, conformers = 0;
  float pcharge = 0.0;

  /* FLAGS */
  int gradflag = 1, verboseflag = 0;

  if ((input = fopen (finput_name, "r")) == NULL)
    {
      fprintf (stderr, "Error. Cant open file %s.\n", finput_name);
      fflush (stderr);
      exit (-1);
    }

  if ((line = (char *) calloc (1, MAX_BUFFER + 2)) == NULL)
    {
      fprintf (stderr, "Error. Cant al memory.\n");
      fflush (stderr);
      exit (-2);
    }


  conformers = 1;
  while (fgets (line, MAX_BUFFER, input))
    {
      if ((strstr (line, "ATOM") != NULL || strstr (line, "HETATM") != NULL)
	  && no_count == 0)
	++molecules;

      if ((line[0] == 'E' && line[1] == 'N' && line[2] == 'D'
	   && line[3] == 'M'))
	{
	  conformers++;
	  no_count = 1;
	}

    /*  if ((line[0] == 'T' && line[1] == 'E' && line[2] == 'R')
	  || (line[0] == 'E' && line[1] == 'N' && line[2] == 'D'
	      && line[3] != 'M'))
	{
	  conformers++;
	  no_count = 1;
	}*/
    }

  rewind (input);

  if ((mols = (MOL2 *) calloc (sizeof (MOL2), 1)) == NULL)
    {
      fprintf (stderr, "Error. Cannot allocate big memory chunk.\n");
      fflush (stderr);
      exit (-2);
    }

  init_molecule (&mols, molecules, conformers);
  *mymol = mols;


  if (gradflag == 1)
    {
      if (verboseflag == 1)
	printf ("Cleaning up the gradients.\n");

      for (j = 0; j < mols->n_atoms; ++j)
	{
	  mols->grads_X[j] = 0.0;
	  mols->grads_Y[j] = 0.0;
	  mols->grads_Z[j] = 0.0;
	}
    }

  current_atom = 0;
  current_conformer = 0;

  while (fgets (line, MAX_BUFFER, input))
    {
#ifdef DEBUG
      fprintf (stderr, "State: Line:%s", line);
      fflush (stderr);
#endif

	  if (line[0] == 'E' && line[1] == 'N' && line[2] == 'D'
	      && line[3] == 'M')
	{

/*		if ( (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {*/
	  current_conformer++;
	  current_atom = 0;
	}
      else if (strstr (line, "ATOM") != NULL
	       || strstr (line, "HETATM") != NULL)
	{
	  sscanf (line, "%*s %*d %4s", tmp_atom);
/*			strncpy(tmp_charge, &line[60], 8);
			strncpy(tmp_radius, &line[56], 4);*/

	  strncpy (tmp_radius, &line[62], 4);
	  strncpy (tmp_charge, &line[55], 7);

	  if (current_conformer == 0)
	    {
	      mols->pcharges[current_atom] = atof (tmp_charge);
	      mols->radius[current_atom] = atof (tmp_radius);

	      if (tmp_atom[0] == 'C'
		  && (tmp_atom[1] != 'l' && tmp_atom[1] != 'L'))
		mols->atoms[current_atom] = 1;
	      else if (tmp_atom[0] == 'O')
		mols->atoms[current_atom] = 2;
	      else if (tmp_atom[0] == 'N')
		mols->atoms[current_atom] = 3;
	      else if (tmp_atom[0] == 'H')
		mols->atoms[current_atom] = 4;
	      else if (tmp_atom[0] == 'P')
		mols->atoms[current_atom] = 5;
	      else if (tmp_atom[0] == 'S')
		mols->atoms[current_atom] = 6;
	      else if (tmp_atom[0] == 'I')
		mols->atoms[current_atom] = 7;
	      else if (tmp_atom[0] == 'B'
		       && (tmp_atom[1] == 'r' || tmp_atom[1] == 'R'))
		mols->atoms[current_atom] = 8;
	      else if (tmp_atom[0] == 'C'
		       && (tmp_atom[1] == 'l' || tmp_atom[1] == 'L'))
		mols->atoms[current_atom] = 9;
	      else if (tmp_atom[0] == 'F' && tmp_atom[1] != 'E'
		       && tmp_atom[1] != 'e')
		mols->atoms[current_atom] = 10;
	      else if (tmp_atom[0] == 'Z'
		       && (tmp_atom[1] == 'N' || tmp_atom[1] == 'n'))
		mols->atoms[current_atom] = 100;
	      else if (tmp_atom[0] == 'M'
		       && (tmp_atom[1] == 'N' || tmp_atom[1] == 'n'))
		mols->atoms[current_atom] = 101;
	      else if (tmp_atom[0] == 'M'
		       && (tmp_atom[1] == 'G' || tmp_atom[1] == 'g'))
		mols->atoms[current_atom] = 102;
	      else if (tmp_atom[0] == 'C' && tmp_atom[1] == '0')
		mols->atoms[current_atom] = 105;
	      else if (tmp_atom[0] == 'F'
		       && (tmp_atom[1] == 'E' || tmp_atom[1] == 'e'))
		mols->atoms[current_atom] = 106;
	      else if (tmp_atom[0] == 'K')
		mols->atoms[current_atom] = 103;
	      else if (tmp_atom[0] == 'N'
		       && (tmp_atom[1] == 'A' || tmp_atom[1] == 'a'))
		mols->atoms[current_atom] = 104;
	      else
		{

		  if (tmp_atom[1] == 'H')
		    mols->atoms[current_atom] = 4;
		  else
		    {
		      mols->atoms[current_atom] = 11;
		      mols->gaff_types[current_atom] = C3;
		      fprintf (stderr,
			       "Reader %i: Error. Do not know the atom type $%s$\n",
			       __LINE__, tmp_atom);
		    }
		}

	      strncpy (res_num, &line[23], 6);
	      mols->res_num[current_atom] = atoi (res_num);
          if (atoi(res_num) != prev_res_num)
          {
              internal_res_num++;
              prev_res_num = atoi(res_num);
          }
          /*fprintf(stderr,"Internal residue number %i\n, Real res number: %i\n", internal_res_num,atoi(res_num));*/
          mols->internal_res_num[current_atom] = internal_res_num;
    
	      strncpy (res_type, &line[17], 3);
	      res_type[3] = '\0';
	      mols->res_type[current_atom] =
		name_to_number_residue (res_type);

	      if (res_type[0] == 'W' && res_type[1] == 'A'
		  && res_type[2] == 'T')
		++wats;

	      strncpy (mols->res_names[current_atom], res_type, 3);
	      strncpy (mols->atom_names[current_atom], &line[12], 4);

	      if (mols->atoms[current_atom] >= 100)
		{
		  fprintf (stderr, "%s is metal",
			   mols->atom_names[current_atom]);
		}



	    }

	  strncpy (myx, &line[29], 10);
	  strncpy (myy, &line[38], 10);
	  strncpy (myz, &line[46], 10);

	  if (current_conformer == 0)
	    {
	      mols->x[current_atom] = atof (myx);
	      mols->y[current_atom] = atof (myy);
	      mols->z[current_atom] = atof (myz);
	    }

	  mols->conformers[current_conformer].x[current_atom] = atof (myx);
	  mols->conformers[current_conformer].y[current_atom] = atof (myy);
	  mols->conformers[current_conformer].z[current_atom] = atof (myz);
	  mols->conformers[current_conformer].pcharges[current_atom] =
	    atof (tmp_charge);
	  mols->conformers[current_conformer].radius[current_atom] =
	    atof (tmp_radius);

	  ++current_atom;


	}
    }


  free (line);
  fclose (input);

  if (mols->n_atoms > 20000)
    {
      fprintf (stderr,
	       "Warning: The system has more than 20000 atoms (%i) and more than 50 water molecules (%i).\nAutomatic topology detection will take ages or may fail to detect the correct bonds. Please remove all not needed water molecules of the system.\n",
	       mols->n_atoms, wats);
      fflush (stderr);
    }

  mol_percieve (&mols);

  AMBER_atom_typing(&mols);

  return 0;
}
