/**
 *
 *      @file assign_pairs.c
 *      @brief Compute vdw pairs list based on bonds 
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 2021/02/12
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


void
assign_pairs (MOL2 ** mymol)
{
  MOL2 *mols = NULL;
  int j = 0, k = 0, cpair = 0;
  float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;


  mols = *mymol;

  mols->vdwpairs_a =
    (int *) realloc (mols->vdwpairs_a,
		     sizeof (int) * mols->n_atoms * mols->n_atoms);
  mols->vdwpairs_b =
    (int *) realloc (mols->vdwpairs_b,
		     sizeof (int) * mols->n_atoms * mols->n_atoms);
  mols->vdw_type =
    (int *) realloc (mols->vdw_type,
		     sizeof (int) * mols->n_atoms * mols->n_atoms);

  for (j = 0; j < mols->n_atoms; ++j)
    {


      for (k = j + 1; k < mols->n_atoms; ++k)
	{
	  dx = mols->x[j] - mols->x[k];
	  dy = mols->y[j] - mols->y[k];
	  dz = mols->z[j] - mols->z[k];
	  dist = dx * dx + dy * dy + dz * dz;
	  if (dist <= 64)
	    {
	      if (get_number_any_bond (mols[0], j + 1, k + 1) == 0)
		{
		  if (q13_bond (mols, j, k) == 0)
		    {
		      mols->vdwpairs_a[cpair] = j;
		      mols->vdwpairs_b[cpair] = k;

		      if (q14_bond (mols, j, k) != 0)
			mols->vdw_type[cpair] = 2;
		      else
			mols->vdw_type[cpair] = 1;

		      cpair++;
		    }
		}
	    }
	}
    }

  mols->n_pairs = cpair;

  *mymol = mols;
}
