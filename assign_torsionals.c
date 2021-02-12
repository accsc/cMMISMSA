/**
 *
 *      @file assign_torsionals.c
 *      @brief Assign GAFF torsionals and impropers
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 2021/02/12
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
assign_torsionals (MOL2 ** mymol)
{
  int i = 0, j = 0, k = 0, l = 0;
  MOL2 *mols = NULL;
  int *vecinos = NULL;
  int *vecinos2 = NULL;
  int matchs[3], nmatchs = 0, ttor = 0;
  int step = 0, bonds1 = 0, bonds2 = 0, bonds3 = 0;
  int vecinos4[16];
  mols = *mymol;

  vecinos = (int *) calloc (sizeof (int *), mols->n_atoms + 1);
  vecinos2 = (int *) calloc (sizeof (int *), mols->n_atoms + 1);



  mols->n_impropers = 0;
  mols->ip = (int *) realloc (mols->ip, sizeof (int) * mols->n_bonds * 6);
  mols->jp = (int *) realloc (mols->jp, sizeof (int) * mols->n_bonds * 6);
  mols->kp = (int *) realloc (mols->kp, sizeof (int) * mols->n_bonds * 6);
  mols->lp = (int *) realloc (mols->lp, sizeof (int) * mols->n_bonds * 6);
  mols->tp = (int *) realloc (mols->tp, sizeof (int) * mols->n_bonds * 6);


  for (i = 0; i < mols->n_atoms; ++i)
    {

      for (j = 0; j < 40; ++j)
	{
	  if (mols->gaff_types[i] == impropers[j][2])
	    {			/* Look for central atom */
	      for (l = 0; l < 8; ++l)
		{
		  vecinos[l] = 0;
		  vecinos2[l] = 0;
		}

	      k = 0;
	      for (l = 0; l < mols->n_bonds; ++l)
		{
		  if (mols->bond_a1[l] == (i + 1))
		    {
		      vecinos[k] = mols->gaff_types[mols->bond_a2[l] - 1];
		      vecinos2[k] = mols->bond_a2[l] - 1;
		      k++;
		    }
		  else if (mols->bond_a2[l] == (i + 1))
		    {
		      vecinos[k] = mols->gaff_types[mols->bond_a1[l] - 1];
		      vecinos2[k] = mols->bond_a1[l] - 1;
		      k++;
		    }
		}

	      nmatchs = 0;
	      matchs[0] = 0;
	      matchs[1] = 0;
	      matchs[2] = 0;

	      for (l = 0; l < k; ++l)
		{
		  /* 71 = wildcard */
		  if (vecinos[l] == impropers[j][0])
		    {
		      matchs[0] = vecinos2[l];
		      nmatchs++;
		      break;
		    }
		}

	      for (l = 0; l < k; ++l)
		{
		  /* 71 = wildcard */
		  if ((vecinos[l] == impropers[j][1])
		      && vecinos2[l] != matchs[0])
		    {
		      matchs[1] = vecinos2[l];
		      nmatchs++;
		      break;
		    }
		}

	      for (l = 0; l < k; ++l)
		{
		  /* 71 = wildcard */
		  if ((vecinos[l] == impropers[j][3])
		      && vecinos2[l] != matchs[1] && vecinos2[l] != matchs[1])
		    {
		      matchs[2] = vecinos2[l];
		      nmatchs++;
		      break;
		    }
		}


	      if (impropers[j][0] == 71)
		{
		  if (vecinos2[0] != matchs[1] && vecinos2[0] != matchs[2])
		    {
		      matchs[0] = vecinos2[0];
		      nmatchs++;
		    }
		  else if (vecinos2[1] != matchs[1]
			   && vecinos2[1] != matchs[2])
		    {
		      matchs[0] = vecinos2[1];
		      nmatchs++;
		    }
		  else if (vecinos2[2] != matchs[1]
			   && vecinos2[2] != matchs[2])
		    {
		      matchs[0] = vecinos2[2];
		      nmatchs++;
		    }
		  else if (vecinos2[3] != matchs[1]
			   && vecinos2[3] != matchs[2])
		    {
		      matchs[0] = vecinos2[3];
		      nmatchs++;
		    }

		}


	      if (impropers[j][1] == 71)
		{
		  if (vecinos2[0] != matchs[0] && vecinos2[0] != matchs[2])
		    {
		      matchs[1] = vecinos2[0];
		      nmatchs++;
		    }
		  else if (vecinos2[1] != matchs[0]
			   && vecinos2[1] != matchs[2])
		    {
		      matchs[1] = vecinos2[1];
		      nmatchs++;
		    }
		  else if (vecinos2[2] != matchs[0]
			   && vecinos2[2] != matchs[2])
		    {
		      matchs[1] = vecinos2[2];
		      nmatchs++;
		    }
		  else if (vecinos2[3] != matchs[0]
			   && vecinos2[3] != matchs[2])
		    {
		      matchs[1] = vecinos2[3];
		      nmatchs++;
		    }

		}



	      if (impropers[j][3] == 71)
		{
		  if (vecinos2[0] != matchs[0] && vecinos2[0] != matchs[1])
		    {
		      matchs[2] = vecinos2[0];
		      nmatchs++;
		    }
		  else if (vecinos2[1] != matchs[0]
			   && vecinos2[1] != matchs[1])
		    {
		      matchs[2] = vecinos2[1];
		      nmatchs++;
		    }
		  else if (vecinos2[2] != matchs[0]
			   && vecinos2[2] != matchs[1])
		    {
		      matchs[2] = vecinos2[2];
		      nmatchs++;
		    }
		  else if (vecinos2[3] != matchs[0]
			   && vecinos2[3] != matchs[1])
		    {
		      matchs[2] = vecinos2[3];
		      nmatchs++;
		    }

		}



	      if (nmatchs == 3)
		{		/* Bingo! We have an improper */

		  mols->ip[mols->n_impropers] = matchs[0];
		  mols->jp[mols->n_impropers] = matchs[1];
		  mols->kp[mols->n_impropers] = i;
		  mols->lp[mols->n_impropers] = matchs[2];
		  mols->tp[mols->n_impropers] = j;
		  mols->n_impropers = mols->n_impropers + 1;

		}

	    }
	}
    }


  printf ("%i impropers found.\n", mols->n_impropers);

  mols->ik = (int *) realloc (mols->ik, sizeof (int) * mols->n_bonds * 6);
  mols->jk = (int *) realloc (mols->jk, sizeof (int) * mols->n_bonds * 6);
  mols->kk = (int *) realloc (mols->kk, sizeof (int) * mols->n_bonds * 6);
  mols->lk = (int *) realloc (mols->lk, sizeof (int) * mols->n_bonds * 6);

  ttor = 0;
  bonds1 = bonds2 = bonds3 = 0;
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

	  for (k = 0; k < bonds2; ++k)
	    {

	      for (step = 0; step < 4; step++)
		vecinos4[step] = 0;

	      bonds3 = 0;
	      for (step = 0; step < mols->n_bonds; ++step)
		{
		  if (mols->bond_a1[step] == (vecinos2[k] + 1))
		    {
		      vecinos4[bonds3] = mols->bond_a2[step] - 1;
		      bonds3++;
		    }
		  else if (mols->bond_a2[step] == (vecinos2[k] + 1))
		    {
		      vecinos4[bonds3] = mols->bond_a1[step] - 1;
		      bonds3++;
		    }
		}

	      for (l = 0; l < bonds3; l++)
		{


		  if (i != vecinos[j] && i != vecinos2[k]
		      && i != vecinos4[l] && vecinos[j] != vecinos2[k]
		      && vecinos[j] != vecinos4[l]
		      && vecinos2[k] != vecinos4[l]
		      && check_dihedral (&mols, ttor, i, vecinos[j],
					 vecinos2[k], vecinos4[l]) == 0)
		    {

		      mols->ik[ttor] = i;
		      mols->jk[ttor] = vecinos[j];
		      mols->kk[ttor] = vecinos2[k];
		      mols->lk[ttor] = vecinos4[l];
		      ttor++;
		    }
		}
	    }
	}
    }


  printf ("%i dihedrals found.\n", ttor);

  mols->n_torsionals = ttor;
  *mymol = mols;
  free (vecinos);
  free (vecinos2);

}
