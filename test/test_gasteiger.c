/**
 *
 *	@file test_gaff_typing.c
 *	@brief Test GAFF typing module
 *
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 2021/02
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

#include <mol/mol.c>
#include <amber/amber.c>


int
main (argc, argv)
     int argc;
     char *argv[];
{

  int i = 0;
  float s = 0.0;
  MOL2 *mol = NULL;
  MultiPDB_reader (&mol, argv[1], 0);
  mol_percieve(&mol);
  assign_gasteiger(&mol);
  for( i = 0; i < mol->n_atoms; i++)
  {
	  	printf("%f\n", mol->pcharges[i]);
		s = s + mol->pcharges[i];
  }
  printf("%f\n",s);

}
