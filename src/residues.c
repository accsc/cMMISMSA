/**
 *
 *      @file residues.c
 *      @brief Handle residues topology and structures
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


int get_number_of_residues(MOL2 *mols)
{
    int i = 0, max_res = -1;
    for( i = 0; i < mols->n_atoms; i++)
    {
        if( mols->internal_res_num[i] > max_res)
            max_res = mols->internal_res_num[i];
    }

    return max_res+1;
}

int get_atoms_for_residue(MOL2 *mols, int res_num)
{
    int i = 0, counter = 0;

    for( i = 0; i < mols->n_atoms; i++)
    {
        if (mols->internal_res_num[i] == res_num)
            ++counter;
    }

    return counter;

}


void process_rings_residue(MOL2 **myparent, int res_num)
{
    int natoms_res = 0, i = 0, j = 0;
    MOL2 *res_mol = NULL, *parent = NULL;
    int *ringer = NULL, *aro = NULL;

    parent = *myparent;
    if ((res_mol = (MOL2 *) calloc (sizeof (MOL2), 1)) == NULL)
    {
      fprintf (stderr, "Error. Cannot allocate big memory chunk.\n");
      fflush (stderr);
      exit (-2);
    }

    natoms_res = get_atoms_for_residue(parent, res_num);
    if(natoms_res == 0)
    return;
    init_molecule (&res_mol, natoms_res, 1);
    j = 0;
    for( i = 0; i< parent->n_atoms; i++)
    {
        if (parent->internal_res_num[i] == res_num)
        {
            res_mol->atoms[j] = parent->atoms[i];
            res_mol->x[j] = parent->x[i];
            res_mol->y[j] = parent->y[i];
            res_mol->z[j] = parent->z[i];
            ++j;
        }
    }

  assign_bond_types (&res_mol);

  ringer = res_mol->ringer;
  aro = res_mol->aromatic;

  get_number_of_rings2 (res_mol[0], &ringer, &aro);
  j = 0;

  for( i = 0; i < parent->n_atoms; i++)
  {
        if (parent->internal_res_num[i] == res_num)
        {
            parent->aromatic[i] = aro[j];
            parent->ringer[i] = ringer[j];
            ++j;
        }
  }
  *myparent = parent;
  cleanup(&res_mol);

}


