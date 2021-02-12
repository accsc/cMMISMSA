/**
 *
 *      @file conformers.c
 *      @brief Provides routines for handling conformers
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 2020/02/12
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




/**
 *
 *	@brief Set the default conformer of the MOL2 structure to the given one
 *	@param mymol MOL2 structure of the molecule
 *	@param conformer number of the conformer starting at 0
 *	@return 0 on success
 */
int
set_default_conformer (MOL2 ** mymol, int conformer)
{

  MOL2 *mols = NULL;
  int i = 0;

  mols = *mymol;

  if (mols->nconformers <= conformer)
    return -1;


  for (i = 0; i < mols->n_atoms; ++i)
    {
      mols->x[i] = mols->conformers[conformer].x[i];
      mols->y[i] = mols->conformers[conformer].y[i];
      mols->z[i] = mols->conformers[conformer].z[i];
/*		mols->pcharges[i] = mols->conformers[conformer].pcharges[i];
		mols->radius[i] = mols->conformers[conformer].radius[i];*/

    }

  mols->default_conformer = conformer;


  *mymol = mols;

  return 0;
}

/**
 *
 *	@brief Syncronizes the coordinates of the current default conformer with the PDB set
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol MOL2 structure of the molecule
 *
 */
int
update_default_conformer (MOL2 ** mymol)
{

  MOL2 *mols = NULL;
  int i = 0;

  mols = *mymol;

  for (i = 0; i < mols->n_atoms; ++i)
    {
      mols->conformers[mols->default_conformer].x[i] = mols->x[i];
      mols->conformers[mols->default_conformer].y[i] = mols->y[i];
      mols->conformers[mols->default_conformer].z[i] = mols->z[i];
    }

  *mymol = mols;

  return 0;
}


