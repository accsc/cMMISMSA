/**
 *
 *      @file assign_angles.c
 *      @brief Compute angle list based on bonds 
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
assign_angles (MOL2 ** mymol)
{
  MOL2 *mols = NULL;
  int cur_angle = 0, i = 0, j = 0, l = 0, step = 0, bonds1 = 0, bonds2 = 0;
  int *vecinos = NULL;
  int *vecinos2 = NULL;


  mols = *mymol;
  vecinos = (int *) calloc (sizeof (int *), mols->n_atoms + 1);
  vecinos2 = (int *) calloc (sizeof (int *), mols->n_atoms + 1);

  mols->ia =
    (int *) realloc (mols->ia, sizeof (int) * mols->n_bonds * mols->n_bonds);
  mols->ja =
    (int *) realloc (mols->ja, sizeof (int) * mols->n_bonds * mols->n_bonds);
  mols->ka =
    (int *) realloc (mols->ka, sizeof (int) * mols->n_bonds * mols->n_bonds);


  for (i = 0; i < mols->n_atoms; ++i)
    {
      bonds1 = 0;
      for (step = 0; step < 4; step++)
	vecinos[step] = 0;

      for (step = 0; step < mols->n_bonds; ++step)
	{
	  if (mols->bond_a1[step] == (i + 1))
	    {
	      vecinos[bonds1] = mols->bond_a2[step] - 1;
	      bonds1++;
	    }
	  else if (mols->bond_a2[step] == (i + 1))
	    {
	      vecinos[bonds1] = mols->bond_a1[step] - 1;
	      bonds1++;
	    }
	}
      for (j = 0; j < bonds1; ++j)
	{
	  bonds2 = 0;
	  for (step = 0; step < 4; step++)
	    vecinos2[step] = 0;

	  for (step = 0; step < mols->n_bonds; ++step)
	    {
	      if (mols->bond_a1[step] == (vecinos[j] + 1))
		{
		  vecinos2[bonds2] = mols->bond_a2[step] - 1;
		  bonds2++;
		}
	      else if (mols->bond_a2[step] == (vecinos[j] + 1))
		{
		  vecinos2[bonds2] = mols->bond_a1[step] - 1;
		  bonds2++;
		}
	    }

	  for (l = 0; l < bonds2; l++)
	    {


	      if (i != vecinos[j] && i != vecinos2[l]
		  && vecinos[j] != vecinos2[l]
		  && check_angle (&mols, cur_angle, i, vecinos[j],
				  vecinos2[l]) == 0)
		{
		  mols->ia[cur_angle] = i;
		  mols->ja[cur_angle] = vecinos[j];
		  mols->ka[cur_angle] = vecinos2[l];
		  cur_angle++;
		}
	    }


	}


    }

  mols->n_angles = cur_angle;
  free (vecinos);
  free (vecinos2);
  *mymol = mols;
}
