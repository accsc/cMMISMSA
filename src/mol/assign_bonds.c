/**
 *
 *      @file assign_bonds.c
 *      @brief Based on distances, assign bond orders
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


void
assign_bond_types (MOL2 ** mymol)
{
  int i = 0, j = 0, k = 0, order = 0;
  float rx = 0.0f, ry = 0.0f, rz = 0.0f, dist = 0.0f;
  MOL2 *mols = NULL;


  mols = *mymol;

  mols->bonds =
    (int *) realloc (mols->bonds, sizeof (int) * mols->n_atoms * 8);
  mols->bond_a1 =
    (int *) realloc (mols->bond_a1, sizeof (int) * mols->n_atoms * 8);
  mols->bond_a2 =
    (int *) realloc (mols->bond_a2, sizeof (int) * mols->n_atoms * 8);


/*
 * These numbers are not totally made-up. Extracted from CCCBDB. http://cccbdb.nist.gov/
 * Some common sense corrections are applied though.
 */
  for (i = 0; i < mols->n_atoms; ++i)
    {
      for (k = i + 1; k < (mols->n_atoms); k++)
	{
	  rx = mols->x[i] - mols->x[k];
	  ry = mols->y[i] - mols->y[k];
	  rz = mols->z[i] - mols->z[k];

	  dist = sqrt (rx * rx + ry * ry + rz * rz);
	  order = 0;


	  /* C - C or C = C or aromatic rings */
	  if (mols->atoms[i] == 1 && mols->atoms[k] == 1)
	    {
	      if (dist <= 1.26 && dist >= 0.5)
		order = 3;
	      else if (dist > 1.26 && dist <= 1.35)
		order = 2;
	      else if (dist > 1.35 && dist <= 1.9)
		order = 1;
	      /* C  - N or C=N */
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 3)
		   || (mols->atoms[i] == 3 && mols->atoms[k] == 1))
	    {
	      if (dist <= 1.85 && dist >= 1.33)
		order = 1;
	      else if (dist < 1.33 && dist >= 1.19)
		order = 2;
	      else if (dist < 1.19)
		order = 3;

	      /* C - O or C=O */
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 2)
		   || (mols->atoms[i] == 2 && mols->atoms[k] == 1))
	    {
	      if (dist > 1.30 && dist <= 1.69)
		order = 1;
	      else if (dist <= 1.30)
		order = 2;
	      /* C - H */
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 4)
		   || (mols->atoms[i] == 4 && mols->atoms[k] == 1))
	    {
	      if (dist <= 1.20)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 3 && mols->atoms[k] == 4)
		   || (mols->atoms[i] == 4 && mols->atoms[k] == 3))
	    {
	      if (dist <= 1.20)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 4 && mols->atoms[k] == 2)
		   || (mols->atoms[i] == 2 && mols->atoms[k] == 4))
	    {
	      if (dist <= 1.20)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 4 && mols->atoms[k] == 6)
		   || (mols->atoms[i] == 6 && mols->atoms[k] == 4))
	    {
	      if (dist <= 1.48)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 6)
		   || (mols->atoms[i] == 6 && mols->atoms[k] == 1))
	    {
	      if (dist <= 1.65)
		order = 2;
	      else if (dist > 1.65 && dist <= 1.91)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 3 && mols->atoms[k] == 6)
		   || (mols->atoms[i] == 6 && mols->atoms[k] == 3))
	    {
	      if (dist <= 2.0)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 2 && mols->atoms[k] == 6)
		   || (mols->atoms[i] == 6 && mols->atoms[k] == 2))
	    {
	      if (dist <= 1.55)
		order = 2;
	      else if (dist > 1.55 && dist <= 1.8)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 2 && mols->atoms[k] == 3)
		   || (mols->atoms[i] == 3 && mols->atoms[k] == 2))
	    {
	      if (dist <= 1.26)
		order = 2;
	      else if (dist <= 1.51 && dist > 1.26)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 10)
		   || (mols->atoms[i] == 10 && mols->atoms[k] == 1))
	    {
	      if (dist <= 1.6)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 8)
		   || (mols->atoms[i] == 8 && mols->atoms[k] == 1))
	    {
	      if (dist <= 2.0)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 9)
		   || (mols->atoms[i] == 9 && mols->atoms[k] == 1))
	    {
	      if (dist <= 1.9)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 1 && mols->atoms[k] == 7)
		   || (mols->atoms[i] == 7 && mols->atoms[k] == 1))
	    {
	      if (dist <= 2.2)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 2 && mols->atoms[k] == 2))
	    {
	      if (dist >= 1.28 && dist <= 2.1)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 3 && mols->atoms[k] == 3))
	    {
	      if (dist >= 1.4 && dist <= 1.8)
		order = 1;
	      else if (dist < 1.4 && dist >= 1.20)
		order = 2;
	      else if (dist < 1.20 && dist >= 0.9)
		order = 3;
	    }
	  else if ((mols->atoms[i] == 5 && mols->atoms[k] == 2)
		   || (mols->atoms[i] == 2 && mols->atoms[k] == 5))
	    {
	      if (dist <= 1.9)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 5 && mols->atoms[k] == 3)
		   || (mols->atoms[i] == 3 && mols->atoms[k] == 5))
	    {
	      if (dist <= 1.99)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 5 && mols->atoms[k] == 1)
		   || (mols->atoms[i] == 1 && mols->atoms[k] == 5))
	    {
	      if (dist <= 2.0)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 5 && mols->atoms[k] == 6)
		   || (mols->atoms[i] == 6 && mols->atoms[k] == 5))
	    {
	      if (dist <= 2.0)
		order = 1;
	    }
	  else if ((mols->atoms[i] == 6 && mols->atoms[k] == 6))
	    {
	      if (dist <= 1.88)
		order = 2;
	      else if (dist > 1.88 && dist < 2.2)
		order = 1;
	    }

	  if (order > 0)
	    {
	      mols->bonds[j] = order;
	      mols->bond_a1[j] = i + 1;
	      mols->bond_a2[j] = k + 1;
	      j++;
	    }


	}
    }

  mols->n_bonds = j;
  *mymol = mols;
}
