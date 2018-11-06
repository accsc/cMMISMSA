/**
 *
 *      @file masks.c
 *      @brief Atom selector masks
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 09/04/2014
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *
 */

#include <string.h>

void select_on_mask(MOL2 **mymol, char *mask, int verbose)
{
	MOL2 *mol = NULL;
	char *token = NULL, *token2_1 = NULL, *token2_2 = NULL;
	char *current_mask = NULL;
	char *global_mask;

	char *local_mask = NULL;

	int first = -1, last = -1, i = 0;
	int total = 0;

	mol = *mymol;
	current_mask = mask;

   	token = strtok_r(mask, ",",&global_mask);
   
	while( token != NULL ) 
   	{
		first = -1; last = -1;
/*      		printf( "%s\n", token );*/
		local_mask = NULL;
		token2_2 = NULL;

		token2_1 = strtok_r(token,"-",&local_mask);
		if( local_mask != NULL)
		 token2_2 = strtok_r(local_mask,"-",&local_mask);

		if( token2_2 != NULL)
		  last = atoi(token2_2);
		first = atoi(token2_1);

		/*
	        if( verbose )
	        {
			fprintf(stderr,"From: %i to: %i\n",first,last);
		}*/

		if( last == -1)
		 last = first;
		for( i = 0; i < mol->n_atoms; i++)
		{
			if( mol->res_num[i] >= first && mol->res_num[i] <= last)
			{
				mol->backbone[i] = 1;
				++total;
			}
		}

      		token = strtok_r(global_mask, ",",&global_mask);
   	}

	if( verbose )
	{
		fprintf(stderr,"Mask: Selected %i atoms\n",total);
		fflush(stderr);
	}

	return;
}

void select_on_mask_atoms(MOL2 **mymol, char *mask, int verbose)
{
	MOL2 *mol = NULL;
	char *token = NULL, *token2_1 = NULL, *token2_2 = NULL;
	char *current_mask = NULL;
	char *global_mask;

	char *local_mask = NULL;

	int first = -1, last = -1, i = 0;
	int total = 0;

	mol = *mymol;
	current_mask = mask;

   	token = strtok_r(mask, ",",&global_mask);
   
	while( token != NULL ) 
   	{
		first = -1; last = -1;
/*      		printf( "%s\n", token );*/
		local_mask = NULL;
		token2_2 = NULL;

		token2_1 = strtok_r(token,"-",&local_mask);
		if( local_mask != NULL)
		 token2_2 = strtok_r(local_mask,"-",&local_mask);

		if( token2_2 != NULL)
		  last = atoi(token2_2);
		first = atoi(token2_1);

		/*
	        if( verbose )
	        {
			fprintf(stderr,"From: %i to: %i\n",first,last);
		}*/

		if( last == -1)
		 last = first;

		if( last >= mol->n_atoms)
			last = mol->n_atoms-1;

		for( i = first-1; i <= last-1; i++)
		{
			
				mol->backbone[i] = 1;
				++total;
		}

      		token = strtok_r(global_mask, ",",&global_mask);
   	}

	if( verbose )
	{
		fprintf(stderr,"Mask: Selected %i atoms\n",total);
		fflush(stderr);
	}

	return;
}

