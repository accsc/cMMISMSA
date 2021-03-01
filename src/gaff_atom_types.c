/**
 *
 *      @file gaff_atom_types.c
 *      @brief Assign GAFF atom types to a molecule
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
GAFF_atom_typing (MOL2 ** mymols)
{
  int i = 0, j = 0, k = 0;
  MOL2 *mols = NULL;
  int *vecinos = NULL;
  int *vecinos2 = NULL;
  int verboseflag = 0;
  int ewg = 0;
  int ewg2 = 0;
  int bondis = 0;
  int flag = 0;
  int flag2 = 0;
  int flag3 = 0;



  mols = *mymols;

  for (i = 0; i < mols->n_atoms; ++i)
    mols->gaff_types[i] = -1;


  vecinos = (int *) calloc (sizeof (int *), mols->n_atoms + 1);
  vecinos2 = (int *) calloc (sizeof (int *), mols->n_atoms + 1);

  for (i = 0; i < mols->n_atoms; ++i)
    {
/* Get vecinos */
      for (j = 0; j < mols->n_atoms; ++j)
	{
	  vecinos[j] = 0;
	  vecinos2[j] = 0;
	}
      k = 0;
      for (j = 0; j < mols->n_bonds; ++j)
	{
	  if (mols->bond_a1[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a2[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	  else if (mols->bond_a2[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a1[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	}

/* Assign types */

      if (mols->atoms[i] == 1)
	{

	  if (mols->aromatic[i] == 1)
	    {
	      mols->gaff_types[i] = CA;
	      if (verboseflag == 1)
		printf ("%i tipo CA.\n", i);
	    }
	  else
	    {
	      if (k == 4)
		{

		  if (mols->ringer[i] == 3)
		    {
		      mols->gaff_types[i] = CX;
		      if (verboseflag == 1)
			printf ("%i tipo CX.\n", i);
		    }
		  else if (mols->ringer[i] == 4)
		    {
		      mols->gaff_types[i] = CY;
		      if (verboseflag == 1)
			printf ("%i tipo CY.\n", i);
		    }
		  else
		    {
		      mols->gaff_types[i] = C3;
		      if (verboseflag == 1)
			printf ("%i tipo C3.\n", i);
		    }
		}
	      else if (k == 3)
		{
		  if (mols->ringer[i] == 3)
		    {
		      mols->gaff_types[i] = CU;
		      if (verboseflag == 1)
			printf ("%i tipo CU.\n", i);
		    }
		  else if (mols->ringer[i] == 4)
		    {
		      mols->gaff_types[i] = CV;
		      if (verboseflag == 1)
			printf ("%i tipo CV.\n", i);
		    }
		  else
		    {
		      ewg = 0;
		      ewg += get_number_bond_by_atom (mols[0], i + 1, 2, 2);
		      ewg += get_number_bond_by_atom (mols[0], i + 1, 2, 6);
		      ewg2 = get_number_bond_by_atom (mols[0], i + 1, 2, 3);
		      if (ewg != 0)
			{
			  mols->gaff_types[i] = C;
			  if (verboseflag == 1)
			    printf ("%i tipo C.\n", i);
			}
		      else if (ewg2 == 3)
			{
			  mols->gaff_types[i] = CZ;
			  if (verboseflag == 1)
			    printf ("%i tipo CZ.\n", i);
			}
		      else
			{
			  ewg = get_bonds (mols[0], i + 1, 2);
			  ewg2 = get_bonds (mols[0], i + 1, 1);
			  if (ewg > 0 && ewg2 > 0)
			    {
			      ewg = 0;
			      for (k = 0; k < 4; ++k)
				if (vecinos2[k] == 1)
				  ewg +=
				    get_bonds (mols[0], vecinos[k] + 1, 2);
			      if (ewg > 0)
				{

				  if (mols->ringer[i] != 0)
				    {
				      mols->gaff_types[i] = CC;
				      if (verboseflag == 1)
					printf ("%i tipo CC.\n", i);
				    }
				  else
				    {
				      mols->gaff_types[i] = CE;
				      if (verboseflag == 1)
					printf ("%i tipo CE.\n", i);
				    }
				}
			      else
				{
				  mols->gaff_types[i] = C2;
				  if (verboseflag == 1)
				    printf ("%i tipo C2.\n", i);
				}
			    }
			  else
			    {
			      mols->gaff_types[i] = C2;
			      if (verboseflag == 1)
				printf ("%i tipo C2.\n", i);
			    }

			}
		    }
		}
	      else if (k == 2)
		{
		  ewg = get_bonds (mols[0], i + 1, 2);
		  ewg2 = get_bonds (mols[0], i + 1, 1);
		  if (ewg > 0 && ewg2 > 0)
		    {
		      ewg = 0;
		      for (k = 0; k < 4; ++k)
			if (vecinos2[k] == 1)
			  ewg += get_bonds (mols[0], vecinos[k] + 1, 2);

		      if (ewg > 0)
			{
			  mols->gaff_types[i] = CG;
			  if (verboseflag == 1)
			    printf ("%i tipo CG.\n", i);
			}
		      else
			{
			  mols->gaff_types[i] = C1;
			  if (verboseflag == 1)
			    printf ("%i tipo C1.\n", i);
			}
		    }
		  else
		    {
		      mols->gaff_types[i] = C1;
		      if (verboseflag == 1)
			printf ("%i tipo C1.\n", i);
		    }
		}
	      else if (k == 1)
		{
		  mols->gaff_types[i] = C1;
		  if (verboseflag == 1)
		    printf ("%i tipo C1.\n", i);
		}
	      else if (verboseflag == 1)
		printf ("Carbon Unknow type!!!\n");
	    }

	}
      else if (mols->atoms[i] == 2)
	{
	  if (k == 1)
	    {
	      mols->gaff_types[i] = O;
	      if (verboseflag == 1)
		printf ("%i tipo O.\n", i);
	    }
	  else
	    {
	      flag = 0;
	      for (k = 0; k < 4; ++k)
		if (mols->atoms[vecinos[k]] == 4)
		  flag = 1;

	      if (flag == 1)
		{
		  mols->gaff_types[i] = OH;
		  if (verboseflag == 1)
		    printf ("%i tipo OH.\n", i);
		}
	      else
		{
		  mols->gaff_types[i] = OS;
		  if (verboseflag == 1)
		    printf ("%i tipo OS.\n", i);
		}


	    }

	}
      else if (mols->atoms[i] == 4)
	{
	  if (k > 1)
	    fprintf (stderr, "Warning: H bonded to many atoms.\n");

	  for (k = 0; k < 1; ++k)
	    {
	      ewg = 0;
	      bondis = get_number_any_bond (mols[0], vecinos[k] + 1, 0);
	      if (mols->atoms[vecinos[k]] == 1)
		{
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 3);
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 2);
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 6);
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 9);
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 8);
		  ewg +=
		    get_number_bond_by_atom (mols[0], vecinos[k] + 1, 0, 7);
		  if (bondis == 4 && ewg == 0)
		    {
		      mols->gaff_types[i] = HC;
		      if (verboseflag == 1)
			printf ("%i tipo HC.\n", i);
		    }
		  else if (bondis == 4 && ewg == 1)
		    {
		      mols->gaff_types[i] = H1;
		      if (verboseflag == 1)
			printf ("%i tipo H1.\n", i);
		    }
		  else if (bondis == 4 && ewg == 2)
		    {
		      mols->gaff_types[i] = H2;
		      if (verboseflag == 1)
			printf ("%i tipo H2.\n", i);
		    }
		  else if (bondis == 4 && ewg == 3)
		    {
		      mols->gaff_types[i] = H3;
		      if (verboseflag == 1)
			printf ("%i tipo H3.\n", i);
		    }
		  else if (bondis == 3 && ewg == 1)
		    {
		      mols->gaff_types[i] = H4;
		      if (verboseflag == 1)
			printf ("%i tipo H4.\n", i);
		    }
		  else if (bondis == 3 && ewg == 2)
		    {
		      mols->gaff_types[i] = H5;
		      if (verboseflag == 1)
			printf ("%i tipo H5.\n", i);
		    }
		  else if (bondis == 5 && ewg == 4)
		    {
		      mols->gaff_types[i] = HX;
		      if (verboseflag == 1)
			printf ("%i tipo HX.\n", i);
		    }
		  else
		    {
		      mols->gaff_types[i] = HA;
		      if (verboseflag == 1)
			printf ("%i tipo HA.\n", i);
		    }
		}
	      else if (mols->atoms[vecinos[k]] == 2)
		{
		  mols->gaff_types[i] = HO;
		  if (verboseflag == 1)
		    printf ("%i tipo HO.\n", i);
		}
	      else if (mols->atoms[vecinos[k]] == 3)
		{
		  mols->gaff_types[i] = HN;
		  if (verboseflag == 1)
		    printf ("%i tipo HN.\n", i);
		}
	      else if (mols->atoms[vecinos[k]] == 5)
		{
		  mols->gaff_types[i] = HP;
		  if (verboseflag == 1)
		    printf ("%i tipo HP.\n", i);
		}
	      else if (mols->atoms[vecinos[k]] == 6)
		{
		  mols->gaff_types[i] = HS;
		  if (verboseflag == 1)
		    printf ("%i tipo HS.\n", i);
		}
	    }

	}
      else if (mols->atoms[i] == 6)
	{

	  if (k == 1)
	    {
	      mols->gaff_types[i] = S;
	      if (verboseflag == 1)
		printf ("%i tipo S.\n", i);

	    }
	  else if (k == 5 || k == 6)
	    {
	      mols->gaff_types[i] = S6;
	      if (verboseflag == 1)
		printf ("%i tipo S6.\n", i);
	    }
	  else if (k == 2)
	    {
	      flag = 0;
	      for (j = 0; j < 4; ++j)
		if (mols->atoms[vecinos[j]] == 4)
		  flag = 1;
	      if (flag == 1)
		{
		  mols->gaff_types[i] = SH;
		  if (verboseflag == 1)
		    printf ("%i tipo SH.\n", i);
		}
	      else if (get_bonds (mols[0], i + 1, 2) > 0
		       || get_bonds (mols[0], i + 1, 3) > 0)
		{
		  mols->gaff_types[i] = S2;
		  if (verboseflag == 1)
		    printf ("%i tipo S2.\n", i);
		}
	      else
		{
		  mols->gaff_types[i] = SS;
		  if (verboseflag == 1)
		    printf ("%i tipo SS.\n", i);
		}
	    }
	  else if (k == 3)
	    {
	      ewg = get_bonds (mols[0], i + 1, 2);
	      ewg2 = get_bonds (mols[0], i + 1, 1);
	      if (ewg > 0 && ewg2 > 0)
		{
		  ewg = 0;
		  for (k = 0; k < 4; ++k)
		    if (vecinos2[k] == 1)
		      ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		  if (ewg > 0)
		    {
		      mols->gaff_types[i] = SX;
		      if (verboseflag == 1)
			printf ("%i tipo SX.\n", i);
		    }
		  else
		    {
		      mols->gaff_types[i] = S4;
		      if (verboseflag == 1)
			printf ("%i tipo S4.\n", i);
		    }
		}
	      else
		{
		  mols->gaff_types[i] = S4;
		  if (verboseflag == 1)
		    printf ("%i tipo S4.\n", i);
		}
	    }
	  else if (k == 4)
	    {
	      ewg = get_bonds (mols[0], i + 1, 2);
	      ewg2 = get_bonds (mols[0], i + 1, 1);
	      if (ewg > 0 && ewg2 > 0)
		{
		  ewg = 0;
		  for (k = 0; k < 4; ++k)
		    if (vecinos2[k] == 1)
		      ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		  if (ewg > 0)
		    {
		      mols->gaff_types[i] = SY;
		      if (verboseflag == 1)
			printf ("%i tipo SY.\n", i);
		    }
		  else
		    {
		      mols->gaff_types[i] = S6;
		      if (verboseflag == 1)
			printf ("%i tipo S6.\n", i);
		    }
		}
	      else
		{
		  mols->gaff_types[i] = S6;
		  if (verboseflag == 1)
		    printf ("%i tipo S6.\n", i);
		}
	    }
	  else
	    {

	      if (verboseflag == 1)
		printf ("Unkown type of S.\n");
	    }


	}
      else if (mols->atoms[i] == 5)
	{
	  if (k == 1)
	    {
	      mols->gaff_types[i] = P2;
	      if (verboseflag == 1)
		printf ("%i tipo P2.\n", i);
	    }
	  else if (k == 2)
	    {
	      ewg = get_bonds (mols[0], i + 1, 2);
	      ewg2 = get_bonds (mols[0], i + 1, 1);
	      if (ewg > 0 && ewg2 > 0)
		{
		  ewg = 0;
		  for (k = 0; k < 4; ++k)
		    if (vecinos2[k] == 1)
		      ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		  if (ewg > 0)
		    {

		      if (mols->ringer[i] != 0)
			{
			  mols->gaff_types[i] = PC;
			  if (verboseflag == 1)
			    printf ("%i tipo PC.\n", i);
			}
		      else
			{
			  mols->gaff_types[i] = PE;
			  if (verboseflag == 1)
			    printf ("%i tipo PE.\n", i);
			}
		    }
		  else
		    {
		      mols->gaff_types[i] = P2;
		      if (verboseflag == 1)
			printf ("%i tipo P2.\n", i);
		    }
		}
	      else
		{
		  mols->gaff_types[i] = P2;
		  if (verboseflag == 1)
		    printf ("%i tipo P2.\n", i);
		}
	    }
	  else if (k == 3)
	    {
	      ewg = 0;
	      ewg += get_number_bond_by_atom (mols[0], i + 1, 2, 2);
	      ewg += get_number_bond_by_atom (mols[0], i + 1, 2, 6);
	      if (ewg != 0)
		{
		  mols->gaff_types[i] = P4;
		  if (verboseflag == 1)
		    printf ("%i tipo P4.\n", i);
		}
	      else
		{
		  ewg = get_bonds (mols[0], i + 1, 2);
		  ewg2 = get_bonds (mols[0], i + 1, 1);
		  if (ewg > 0 && ewg2 > 0)
		    {
		      ewg = 0;
		      for (k = 0; k < 4; ++k)
			if (vecinos2[k] == 1)
			  ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		      if (ewg > 0)
			{
			  mols->gaff_types[i] = PX;
			  if (verboseflag == 1)
			    printf ("%i tipo PX.\n", i);
			}
		      else
			{
			  mols->gaff_types[i] = P3;
			  if (verboseflag == 1)
			    printf ("%i tipo P3.\n", i);
			}
		    }
		  else
		    {
		      mols->gaff_types[i] = P3;
		      if (verboseflag == 1)
			printf ("%i tipo P3.\n", i);
		    }
		}
	    }
	  else if (k == 4)
	    {
	      ewg = get_bonds (mols[0], i + 1, 2);
	      ewg2 = get_bonds (mols[0], i + 1, 1);
	      if (ewg > 0 && ewg2 > 0)
		{
		  ewg = 0;
		  for (k = 0; k < 4; ++k)
		    if (vecinos2[k] == 1)
		      ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		  if (ewg > 0)
		    {
		      mols->gaff_types[i] = PY;
		      if (verboseflag == 1)
			printf ("%i tipo PY.\n", i);
		    }
		  else
		    {
		      mols->gaff_types[i] = P5;
		      if (verboseflag == 1)
			printf ("%i tipo P5.\n", i);
		    }
		}
	      else
		{
		  mols->gaff_types[i] = P5;
		  if (verboseflag == 1)
		    printf ("%i tipo P5.\n", i);
		}
	    }
	  else if (k == 5 || k == 6)
	    {
	      mols->gaff_types[i] = P5;
	      if (verboseflag == 1)
		printf ("%i tipo P5.\n", i);
	    }
	  else if (verboseflag == 1)
	    printf ("P Unknow type.\n");
	}
      else if (mols->atoms[i] == 10)
	{
	  mols->gaff_types[i] = F;
	  if (verboseflag == 1)
	    printf ("%i tipo F.\n", i);
	}
      else if (mols->atoms[i] == 9)
	{
	  mols->gaff_types[i] = CL;
	  if (verboseflag == 1)
	    printf ("%i tipo CL.\n", i);
	}
      else if (mols->atoms[i] == 8)
	{
	  mols->gaff_types[i] = BR;
	  if (verboseflag == 1)
	    printf ("%i tipo BR.\n", i);
	}
      else if (mols->atoms[i] == 7)
	{
	  mols->gaff_types[i] = I;
	  if (verboseflag == 1)
	    printf ("%i tipo I.\n", i);
	}
    }


/* Another loop for N. I need carbonyl checks before typing amides*/
/* Lazzy-crap(tm) technique */

  for (i = 0; i < mols->n_atoms; ++i)
    {
/* Get vecinos */
      for (j = 0; j < mols->n_atoms; ++j)
	{
	  vecinos[j] = 0;
	  vecinos2[j] = 0;
	}
      k = 0;
      for (j = 0; j < mols->n_bonds; ++j)
	{
	  if (mols->bond_a1[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a2[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	  else if (mols->bond_a2[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a1[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	}

      if (mols->atoms[i] == 3)
	{
	  if (k == 1)
	    {
	      mols->gaff_types[i] = N1;
	      if (verboseflag == 1)
		printf ("%i tipo N1.\n", i);
	    }
	  else if (k == 2)
	    {
	      ewg = get_bonds (mols[0], i + 1, 2);
	      if (mols->aromatic[i] == 1)
		{
		  mols->gaff_types[i] = NB;
		  if (verboseflag == 1)
		    printf ("%i tipo NB.\n", i);
		}
	      else if (get_bonds (mols[0], i + 1, 3) > 0 || ewg == 2)
		{
		  mols->gaff_types[i] = N1;
		  if (verboseflag == 1)
		    printf ("%i tipo N1.\n", i);
		}
	      else
		{
		  ewg2 = get_bonds (mols[0], i + 1, 1);
		  if (ewg > 0 && ewg2 > 0)
		    {
		      ewg = 0;
		      for (k = 0; k < 4; ++k)
			if (vecinos2[k] == 1)
			  ewg += get_bonds (mols[0], vecinos[k] + 1, 2);
		      if (ewg > 0)
			{
			  if (mols->ringer[i] != 0)
			    {
			      mols->gaff_types[i] = NC;
			      if (verboseflag == 1)
				printf ("%i tipo NC.\n", i);
			    }
			  else
			    {
			      mols->gaff_types[i] = NE;
			      if (verboseflag == 1)
				printf ("%i tipo NE.\n", i);
			    }
			}
		      else
			{
			  mols->gaff_types[i] = N2;
			  if (verboseflag == 1)
			    printf ("%i tipo N2.\n", i);
			}
		    }
		  else
		    {
		      mols->gaff_types[i] = N2;
		      if (verboseflag == 1)
			printf ("%i tipo N2.\n", i);
		    }
		}
	    }
	  else if (k == 3)
	    {
	      if (mols->aromatic[i] == 1)
		{
		  mols->gaff_types[i] = NA;
		  if (verboseflag == 1)
		    printf ("%i tipo NA.\n", i);
		}
	      else
		{

		  flag2 = 0;
		  flag = 0;
		  flag3 = 0;
		  ewg = 0;
		  for (j = 0; j < k; ++j)
		    {
	              if( mols->aromatic[vecinos[j]] != 0)
			     ewg = 1;
		      if (mols->atoms[vecinos[j]] == 2)
			flag2++;
		      else if (mols->gaff_types[vecinos[j]] == C)
			flag = 1;
		      else if (vecinos2[j] == 2
			       && (mols->gaff_types[vecinos[j]] == C
				   || mols->atoms[vecinos[j]] == 3
				   || mols->atoms[vecinos[j]] == 1))
			flag3++;
		    }

		  if (flag2 == 2 && k == 3)
		    {
		      mols->gaff_types[i] = NO;
		      if (verboseflag == 1)
			printf ("%i tipo NO.\n", i);
		    }
		  else
		    {
		      if (flag == 1)
			{
			  mols->gaff_types[i] = N;
			  if (verboseflag == 1)
			    printf ("%i tipo N.\n", i);
			}
		      else if (flag3 == 1 || ewg == 1)
			{

			  mols->gaff_types[i] = NH;
			  if (verboseflag == 1)
			    printf ("%i tipo NH.\n", i);

			}
		      else
			{
			  mols->gaff_types[i] = N3;
			  if (verboseflag == 1)
			    printf ("%i tipo N3.\n", i);
			}
		    }
		}

	    }
	  else if (k == 4)
	    {
	      mols->gaff_types[i] = N4;
	      if (verboseflag == 1)
		printf ("%i tipo N4.\n", i);
	    }
	  else if (verboseflag == 1)
	    printf ("N  Unknow type.\n");
	}
    }


#ifdef DEBUG
  for (j = 0; j < mols->n_atoms; ++j)
    printf ("%i - %i.\n", j, mols->gaff_types[j]);

#endif
/* Fixes */


  for (i = 0; i < mols->n_atoms; ++i)
    {
      flag = 0;
      flag2 = 0;
      for (j = 0; j < mols->n_atoms; ++j)
	{
	  vecinos[j] = 0;
	  vecinos2[j] = 0;
	}
      k = 0;
      for (j = 0; j < mols->n_bonds; ++j)
	{
	  if (mols->bond_a1[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a2[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	  else if (mols->bond_a2[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a1[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	}

      flag3 = 0;
      flag2 = 0;
      flag = 0;

      if (mols->gaff_types[i] == C || mols->gaff_types[i] == C2)
	{
	  for (j = 0; j < k; j++)
	    {
	      if (mols->gaff_types[vecinos[j]] == O)
		flag++;
	      if (mols->gaff_types[vecinos[j]] == S)
		flag3++;

	      if (mols->atoms[vecinos[j]] == 2
		  && mols->ringer[vecinos[j]] <= 0)
		flag2++;
	    }

	  if (mols->gaff_types[i] == C
	      && (flag >= 2 || (mols->ringer[i] > 0 && flag2 == 0)))
	    mols->gaff_types[i] = C2;
	  if (mols->gaff_types[i] == C2 && flag3 > 0)
	    mols->gaff_types[i] = C;

	}
      flag = flag2 = flag3 = 0;

      if (mols->gaff_types[i] == N3 || mols->gaff_types[i] == S2
	  || mols->gaff_types[i] == N || mols->gaff_types[i] == CC)
	{
	  for (j = 0; j < k; j++)
	    {
	      if (mols->gaff_types[vecinos[j]] == CA)
		flag++;

	      if (mols->gaff_types[vecinos[j]] == C2
		  || mols->gaff_types[vecinos[j]] == C
		  || mols->gaff_types[vecinos[j]] == N2
		  || mols->gaff_types[vecinos[j]] == S2
		  || mols->gaff_types[vecinos[j]] == SS)
		flag2++;
	      if (mols->aromatic[vecinos[j]] > 0)
		flag3++;

	    }

	  if (mols->gaff_types[i] == N3 && flag == 1 && flag2 == 1)
	    mols->gaff_types[i] = NA;
	  if (mols->gaff_types[i] == N3 && flag2 >= 2)
	    mols->gaff_types[i] = NA;
	  if (mols->gaff_types[i] == S2 && (flag2 == 2 || flag == 2))
	    mols->gaff_types[i] = SS;
	  if (mols->gaff_types[i] == S2 && flag2 == 1 && flag == 1)
	    mols->gaff_types[i] = SS;
	  if (mols->gaff_types[i] == CC && flag2 + flag >= 2)
	    mols->gaff_types[i] = C2;
	  /*if (mols->gaff_types[i] == N && flag3 == 0)
	   *            mols->gaff_types[i] = NH; */
	}

    }

  for (i = 0; i < mols->n_atoms; ++i)
    {

      if (mols->atoms[i] != 4)
	continue;

      flag = 0;
      flag2 = 0;
      for (j = 0; j < mols->n_atoms; ++j)
	{
	  vecinos[j] = 0;
	  vecinos2[j] = 0;
	}

      k = 0;
      for (j = 0; j < mols->n_bonds; ++j)
	{
	  if (mols->bond_a1[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a2[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	  else if (mols->bond_a2[j] == (i + 1))
	    {
	      vecinos[k] = mols->bond_a1[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	}
      flag = vecinos[0];

      k = 0;
      for (j = 0; j < mols->n_bonds; ++j)
	{
	  if (mols->bond_a1[j] == (flag + 1))
	    {
	      vecinos[k] = mols->bond_a2[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	  else if (mols->bond_a2[j] == (flag + 1))
	    {
	      vecinos[k] = mols->bond_a1[j] - 1;
	      vecinos2[k] = mols->bonds[j];
	      k++;
	    }
	}

      for (j = 0; j < k; j++)
	{
	  if (mols->gaff_types[vecinos[j]] == N4)
	    flag2++;
	}

      if (k == 4 && flag2 > 0 && mols->atoms[flag] == 1)
	mols->gaff_types[i] = HX;

    }


  free (vecinos);
  free (vecinos2);

  *mymols = mols;
}


void assign_gaff_vdw_types(MOL2 **mymol)
{
    MOL2 *mols = NULL;
    int i = 0;

    mols = *mymol;

    for ( i = 0; i < mols->n_atoms; i++)
    {
    /* Assign epsilon (well depth) parameter from GAFF */
    mols->vdw_parm1[i] = gaff_vdw_epsilon[ mols->gaff_types[i] -1 ];
    /* Assign r0 (equilibrium distance) parameter from GAFF */
    mols->vdw_parm2[i] = gaff_vdw_r0[ mols->gaff_types[i] -1 ];
    }

    *mymol = mols;

}
