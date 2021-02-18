/**
 *
 *
 *
 *      @file amber.c
 *      @brief Read a AMBER top file to a MOL2 structure
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@date 01/10/2010
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

#include "amber.h"


int TOP_reader(MOL2 **mymol, char *finput_name, TOP **mytop, int verbose)
{

	int i = 0, j = 0, k = 0, l = 0;
	FILE *input = NULL;
	char *line = NULL;
	MOL2 *mol = NULL;
	TOP *topo = NULL;

	int natoms = 0, nres = 0, conformers = 1, ntypes = 0;
	char tresnum[9], tcharge[17], tresname[5];
	int current_res = 0;

	int *res_pointers = NULL;
	char **res_labels = NULL;
	char **atom_names = NULL;
	float *pcharges = NULL;
	int real_atoms = 0;
	int *atom_type_index = NULL;

	char keyword[100];
	int *nb_index = NULL;
	float *lj_A = NULL, *lj_B = NULL;

	int init_flag = 0, resnum_flag = 0, resname_flag = 0;
	int atom_name_flag = 0, pcharges_flag = 0;

        if ( (input = fopen(finput_name, "r")) == NULL) {
                fprintf(stderr, "TOP_reader: Error. Cant open file %s.\n", finput_name);
                fflush(stderr);
                return -1;
        }

        if ( (line = (char*)calloc(1, 1024 + 2)) == NULL) {
                fprintf(stderr, "TOP_reader: Error. Cant allocate memory.\n");
                fflush(stderr);
             	return -2;
        }

        if( (topo = (TOP *) calloc(sizeof(TOP),1)) == NULL)
        {
                fprintf(stderr,"TOP_reader: Cannot allocate memory for the structure\n");
                return -2;
        }


        if( (mol = (MOL2 *) calloc(sizeof(MOL2),1)) == NULL)
        {
                fprintf(stderr,"TOP_reader: Cannot allocate memory for the structure\n");
                return -2;
        }


	if( verbose )
		fprintf(stderr,"TOP_reader: Scanning file...");


        while ( fgets(line, 1024, input) ) 
	{
		if( line[0] == '%' && line[1] == 'F' && line[2] == 'L' && line[3] == 'A' && line[4] == 'G')
		{
			sscanf(line,"%*s %99s",&keyword);
			fgets(line, 1024, input); /* This is the format line */

			if( strstr(keyword,"POINTERS") != NULL && strspn(keyword,"POINTERS") == 8)
			{
	                        fgets(line, 1024, input); 
				strncpy(tresnum,line,8);
				tresnum[8] = '\0';
				natoms = atoi(tresnum);
                                strncpy(tresnum,&line[8],8);
                                tresnum[8] = '\0';
                                ntypes = atoi(tresnum);
                                fgets(line, 1024, input);
                                strncpy(tresnum,&line[8],8);
                                tresnum[8] = '\0';
                                nres = atoi(tresnum);

				if( natoms <= 0 || ntypes <= 0 || nres <= 0)
				{
					fprintf(stderr,"Corrupted topology file\n");
					fflush(stderr);
					return -3;
				}

				topo->n_atoms = natoms;
				topo->ntypes = ntypes;
				topo->nb_index = (int *) calloc(sizeof(int), ntypes*ntypes);
				topo->lj_A = (float *) calloc(sizeof(float), ntypes*(ntypes+1)/2);
				topo->lj_B = (float *) calloc(sizeof(float), ntypes*(ntypes+1)/2);
				atom_type_index = (int *) calloc(sizeof(int), natoms);
				topo->excluded_list = (int *) calloc(sizeof(int), natoms);

				/* Initialize the structure for the following sections */
				/* If we find another structure before we throw an error */

				pcharges = (float*)calloc(sizeof(float), natoms + 1);
				atom_names = (char**)calloc(sizeof(char*), natoms + 1);
                                for( i = 0; i < natoms; i++)
                                {
                                        atom_names[i] = (char *) calloc(sizeof(char),5);
				}


			}else if( strstr(keyword,"ATOM_NAME") != 0 && strspn(keyword,"ATOM_NAME") == 9){
				if( natoms <= 0)
				{
					fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
					fflush(stderr);
					return -3;
				}

				j = (int) ceilf((float) natoms / 20.0f);
				k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
				{
                                        fgets(line, 1024, input); 
					for( k = 0; k < 20; k++)
					{
					 strncpy(atom_names[(i*20)+k],&line[k*4],4);
					 atom_names[(i*20)+k][4] = '\0';
					}

				}
                                fgets(line, 1024, input);  /* Last line */
				if( (j = natoms % 20) == 0){
					j = 20; /* Complete then */
				}
                                for( k = 0; k < j; k++)
				{
                                        strncpy(atom_names[(i*20)+k],&line[k*4],4);
                                        atom_names[(i*20)+k][4] = '\0';
				}
				atom_name_flag = 1;
			/* End of atom name reader */
			}else if( strstr(keyword,"RESIDUE_LABEL") != NULL && strspn(keyword,"RESIDUE_LABEL") == 13){
                                if( nres <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }


                                res_labels = (char **) calloc(sizeof(char*), nres+1);
                                for( i = 0; i < nres; i++)
                                {
                                        res_labels[i] = (char *) calloc(sizeof(char),5);
                                }



                                j = (int) ceilf((float) nres / 20.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 20; k++)
                                        {
                                         strncpy(res_labels[(i*20)+k],&line[k*4],4);
					 res_labels[(i*20)+k][4] = '\0';
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = nres % 20) == 0){
                                        j = 20; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(res_labels[(i*20)+k],&line[k*4],4);
                                         res_labels[(i*20)+k][4] = '\0';
                                }

				resname_flag = 1;
                        /* End of residue name reader */
			}else if( strstr(keyword,"RESIDUE_POINTER") != NULL && strspn(keyword,"RESIDUE_POINTER") == 15){
                                if( nres <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                res_pointers = (int *) calloc(sizeof(int), nres+1);

				current_res = 0;
                                j = (int) ceilf((float) nres / 10.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 10; k++)
                                        {
                                         strncpy(tresnum,&line[k*8],8);
					 tresnum[8] = '\0';
                                         res_pointers[current_res] = atoi(tresnum);
					 ++current_res;
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = nres % 10) == 0){
                                        j = 10; /* Complete then */
                                }
                                        for( k = 0; k < j; k++)
                                        {
                                         strncpy(tresnum,&line[k*8],8);
                                         tresnum[8] = '\0';
                                         res_pointers[current_res] = atoi(tresnum);
                                         ++current_res;
                                        }

				res_pointers[current_res] = natoms+1;

				resnum_flag = 1;

                        /* End of res name reader */
			}else if( strstr(keyword,"CHARGE") != NULL && strspn(keyword,"CHARGE") == 6){
                                if( natoms <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                j = (int) ceilf((float) natoms / 5.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 5; k++)
                                        {
                                         strncpy(tcharge,&line[k*16],16);
                                         tcharge[16] = '\0';
					 pcharges[(i*5)+k] = atof(tcharge) / 18.2223f; 
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = natoms % 5) == 0){
                                        j = 5; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(tcharge,&line[k*16],16);
                                        tcharge[16] = '\0';
                                        pcharges[(i*5)+k] = atof(tcharge) / 18.2223f;
                                }
                        /* End of charges reader */
			}else if( strstr(keyword,"ATOM_TYPE_INDEX") != NULL && strspn(keyword,"ATOM_TYPE_INDEX") == 15){
                                if( natoms <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                j = (int) ceilf((float) natoms / 10.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 10; k++)
                                        {
                                         strncpy(tresnum,&line[k*8],8);
                                         tresnum[8] = '\0';
                                         atom_type_index[(i*10)+k] = atoi(tresnum);
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = natoms % 10) == 0){
                                        j = 10; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(tresnum,&line[k*8],8);
                                        tresnum[8] = '\0';
                                        atom_type_index[(i*10)+k] = atoi(tresnum);
                                }
                        /* End of atomtype reader */
                        }else if( strstr(keyword,"NONBONDED_PARM_INDEX") != NULL && strspn(keyword,"NONBONDED_PARM_INDEX") == 20){
                                if( ntypes <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                j = (int) ceilf((float) ntypes*ntypes / 10.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 10; k++)
                                        {
                                         strncpy(tresnum,&line[k*8],8);
                                         tresnum[8] = '\0';
					 topo->nb_index[(i*10)+k] = atoi(tresnum);
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = (ntypes*ntypes) % 10) == 0){
                                        j = 10; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(tresnum,&line[k*8],8);
                                        tresnum[8] = '\0';
                                        topo->nb_index[(i*10)+k] = atoi(tresnum);
                                }
                        /* End of nonbonded parm index reader */
			}else if( strstr(keyword,"LENNARD_JONES_ACOEF") != NULL && strspn(keyword,"LENNARD_JONES_ACOEF") == 19){
                                if( ntypes <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                j = (int) ceilf((float) (ntypes*(ntypes+1)/2.0f) / 5.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 5; k++)
                                        {
                                         strncpy(tcharge,&line[k*16],16);
                                         tcharge[16] = '\0';
                                         topo->lj_A[(i*5)+k] = atof(tcharge);
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = (ntypes*(ntypes+1)/2) % 5) == 0){
                                        j = 5; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(tcharge,&line[k*16],16);
                                        tcharge[16] = '\0';
                                        topo->lj_A[(i*5)+k] = atof(tcharge);
                                }
                        /* End of LJ_A parm index reader */
                        }else if( strstr(keyword,"LENNARD_JONES_BCOEF") != NULL && strspn(keyword,"LENNARD_JONES_BCOEF") == 19){
                                if( ntypes <= 0)
                                {
                                        fprintf(stderr,"Error. POINTERS section must be placed before any other atom section\n");
                                        fflush(stderr);
                                        return -3;
                                }

                                j = (int) ceilf((float) (ntypes*(ntypes+1)/2.0f) / 5.0f);
                                k = 0;
                                for( i = 0; i < j-1; i++) /* Till the last line */
                                {
                                        fgets(line, 1024, input);
                                        for( k = 0; k < 5; k++)
                                        {
                                         strncpy(tcharge,&line[k*16],16);
                                         tcharge[16] = '\0';
                                         topo->lj_B[(i*5)+k] = atof(tcharge);
                                        }

                                }
                                fgets(line, 1024, input);  /* Last line */
                                if( (j = (ntypes*(ntypes+1)/2) % 5) == 0){
                                        j = 5; /* Complete then */
                                }
                                for( k = 0; k < j; k++)
                                {
                                        strncpy(tcharge,&line[k*16],16);
                                        tcharge[16] = '\0';
                                        topo->lj_B[(i*5)+k] = atof(tcharge);
                                }
                        /* End of LJ_B parm index reader */
                        }
		}
        }

	if( verbose )
		fprintf(stderr," done.\nTOP_reader: Topology contains Natoms: %i. Nres: %i. Ntypes: %i\n",natoms,nres,ntypes);

	if( resnum_flag )
	{
		for( i = 0 ; i < current_res; i++)
		{
		        for( j = res_pointers[i]-1; j < res_pointers[i+1]-1; j++)
		        {
/*		                mol->res_num[j] = i+1;*/
				if( resname_flag )
				{
				  if( (res_labels[i][0] == 'W' && res_labels[i][1] == 'A' && res_labels[i][2] == 'T') || (res_labels[i][0] == 'N' && res_labels[i][1] == 'a' && res_labels[i][2] == '+') || (res_labels[i][0] == 'C' && res_labels[i][1] == 'l' && res_labels[i][2] == '-') )
					  topo->excluded_list[j] = 1;
				  else
					  real_atoms++;
/*				  strncpy(mol->res_names[j],res_labels[i],3);*/
				}
		        }
		}

/*		free(res_pointers);*/
	}

	if( verbose )
	{
		fprintf(stderr,"TOP_reader: Number of atoms after water/ion exclusion: %i\n",real_atoms);
	}

        mol->n_atoms = real_atoms;
        mol->x = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->y = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->z = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->pcharges = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->radius = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->atoms = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->ringer = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->aromatic = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->bond_dist = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->grads_X = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->grads_Y = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->grads_Z = (float*)calloc(sizeof(float), mol->n_atoms);
        mol->backbone = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->exclude = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->selection = (int*)calloc(sizeof(int), mol->n_atoms + 10);
        mol->fragment_flag = (int*)calloc(sizeof(int), mol->n_atoms);
        mol->conformers = (CONFORMER*)calloc(sizeof(CONFORMER), conformers);
        mol->n_fragments = 0;

        mol->res_num = (int*)calloc(sizeof(int), mol->n_atoms + 1);
        mol->internal_res_num = (int*)calloc(sizeof(int), mol->n_atoms + 1);
        mol->res_names = (char**)calloc(sizeof(char*), mol->n_atoms + 1);
        mol->atom_names = (char**)calloc(sizeof(char*), mol->n_atoms + 1);
        mol->res_type = (int*)calloc(sizeof(int), mol->n_atoms + 1);


        for( i = 0; i < mol->n_atoms; i++)
        {
                mol->res_names[i] = (char *) calloc(sizeof(char),5);
                mol->atom_names[i] = (char *) calloc(sizeof(char),5);
        }

        mol->gaff_types = (int*)calloc(sizeof(int), mol->n_atoms + 1);
        mol->ism_types = (int *)calloc(sizeof(int), mol->n_atoms + 1);
        mol->hyde_types = (int *)calloc(sizeof(int), mol->n_atoms + 1);

        mol->ism_selection = (int *)calloc(sizeof(int), mol->n_atoms + 1);
        mol->hyde_selection = (int *)calloc(sizeof(int), mol->n_atoms + 1);

        mol->vdw_selection = (int *)calloc(sizeof(int), mol->n_atoms + 1);

        topo->atom_type_index = (int *) calloc(sizeof(int), mol->n_atoms);

	/* TODO. this does not work when there are other residues after a WAT block */
        if( resnum_flag )
        {
		k = 0;
                for( i = 0 ; i < current_res; i++)
                {
                     if( (res_labels[i][0] == 'W' && res_labels[i][1] == 'A' && res_labels[i][2] == 'T') || (res_labels[i][0] == 'N' && res_labels[i][1] == 'a' && res_labels[i][2] == '+') || (res_labels[i][0] == 'C' && res_labels[i][1] == 'l' && res_labels[i][2] == '-'))
			continue;
			++k;
                        for( j = res_pointers[i]-1; j < res_pointers[i+1]-1; j++)
                        {
                              mol->res_num[j] = k;
                              mol->internal_res_num[j] = k;
                              if( resname_flag )
                              {
                                strncpy(mol->res_names[j],res_labels[i],3);
                              }
                        }
                }

              free(res_pointers);
        }

	j = 0;
	for( i = 0; i < natoms; i++)
	{
		if( !topo->excluded_list[i] )
		{
			strncpy(mol->atom_names[j],atom_names[i],4);
			topo->atom_type_index[j] = atom_type_index[i];
			mol->pcharges[j] = pcharges[i];
			++j;
		}
		free(atom_names[i]);
	}
	free(atom_type_index);
	free(pcharges);
	free(atom_names);

	if( resname_flag )
	{
		for( i = 0; i < nres; i++)
		{
			free(res_labels[i]);
		}
		free(res_labels);
	}

	if( atom_name_flag )
	{
		for( i = 0; i < mol->n_atoms; i++)
		{
			if ( mol->atom_names[i][0] == 'C' && ( mol->atom_names[i][1] != 'l' && mol->atom_names[i][1] != 'L'))
                		mol->atoms[i] = 1;
	                else if ( mol->atom_names[i][0] == 'O')
	                        mol->atoms[i] = 2;
	                else if ( mol->atom_names[i][0] == 'N')
                                mol->atoms[i] = 3;
                        else if ( mol->atom_names[i][0] == 'H')
                                mol->atoms[i] = 4;
                        else if ( mol->atom_names[i][0] == 'P')
                                mol->atoms[i] = 5;
                        else if ( mol->atom_names[i][0] == 'S')
                                mol->atoms[i] = 6;
                        else if ( mol->atom_names[i][0] == 'I')
                                mol->atoms[i] = 7;
                        else if ( mol->atom_names[i][0] == 'B' && (mol->atom_names[i][1] == 'r' || mol->atom_names[i][1] == 'R'))
                                mol->atoms[i] = 8;
                        else if ( mol->atom_names[i][0] == 'C' && (mol->atom_names[i][1] == 'l' || mol->atom_names[i][1] == 'L'))
                                mol->atoms[i] = 9;
                        else if( mol->atom_names[i][0] == 'F' &&  mol->atom_names[i][1] != 'E' && mol->atom_names[i][1] != 'e')
                                mol->atoms[i] = 10;
                        else if(mol->atom_names[i][0] == 'Z' && ( mol->atom_names[i][1] == 'N' || mol->atom_names[i][1] == 'n' ))
                                mol->atoms[i] = 100;
                        else if(mol->atom_names[i][0] == 'M' && (mol->atom_names[i][1] == 'N' || mol->atom_names[i][1] == 'n') )
                                mol->atoms[i] = 101;
                        else if(mol->atom_names[i][0] == 'M' && (mol->atom_names[i][1] == 'G' || mol->atom_names[i][1] == 'g' ))
                                mol->atoms[i] = 102;
                        else if(mol->atom_names[i][0] == 'C' && mol->atom_names[i][1] == '0')
                                mol->atoms[i] = 105;
                        else if(mol->atom_names[i][0] == 'F' && (mol->atom_names[i][1] == 'E' || mol->atom_names[i][1] == 'e'))
                                mol->atoms[i] = 106;
                        else if(mol->atom_names[i][0] == 'K')
                                mol->atoms[i] = 103;
                        else if(mol->atom_names[i][0] == 'N' && (mol->atom_names[i][1] == 'A' || mol->atom_names[i][1] == 'a'))
                                mol->atoms[i] = 104;
                        else{

                               if ( mol->atom_names[i][1] == 'H')
                                       mol->atoms[i] = 4;
                               else{
                                       mol->atoms[i] = 11;
                                       mol->gaff_types[i] = C3;
/*                                       fprintf(stderr, "Reader %i: Error. Do not know the atom type $%s$\n", __LINE__,mol->atom_names[i]);*/
                               }
                        }
		}
	}

/*        mol_percieve(&mol, 0);*/  /* I only trust my own connectivity */

        *mymol = mol;
	*mytop = topo;

        for ( i = 0; i < conformers; ++i) {
                mol->conformers[i].x = (float*)calloc(sizeof(float), mol->n_atoms);
                mol->conformers[i].y = (float*)calloc(sizeof(float), mol->n_atoms);
                mol->conformers[i].z = (float*)calloc(sizeof(float), mol->n_atoms);
                mol->conformers[i].pcharges = (float*)calloc(sizeof(float), mol->n_atoms);
                mol->conformers[i].radius = (float*)calloc(sizeof(float), mol->n_atoms);
        }

	fclose(input);
	free(line);
	return 0;
}


int load_crd(MOL2 **mymol, char *fname, TOP *topo)
{

	MOL2 *mol = NULL;
	FILE *input = NULL;
	char *line = NULL;
	int natoms = 0;
	int i = 0, j = 0, k = 0;
	char tcoor[13];
	int index = 0, r = 0;

	mol = *mymol;

	if( (input = fopen(fname,"rb")) == NULL)
	{
		fprintf(stderr,"Cannot open the crd file\n");
		fflush(stderr);
		return -1;
	}

	if( (line = (char *) calloc(sizeof(char),1024)) == NULL)
	{
		fprintf(stderr,"Cannot allocate buffer\n");
		fflush(stderr);
		return 2;
	}

        fgets(line, 1024, input); 
        fgets(line, 1024, input); 

	strncpy(tcoor,line,6);
	tcoor[6] = '\0';
	natoms = atoi(tcoor);

	fprintf(stderr,"Topology atoms: %i | Coordinates in file for atoms: %i\n",topo->n_atoms,natoms);
	fflush(stderr);
	
	if( natoms != topo->n_atoms)
	{
		fprintf(stderr,"Different number of atoms!");
		fflush(stderr);
		free(line);
		return -3;
	}

	j = (int) ceil((float) natoms*3.0f/6.0f);

        k = 0;
        for( i = 0; i < j-1; i++) 
        {
                fgets(line, 1024, input);
                if( line == NULL || strlen(line) < (6*12))
                {
                  fprintf(stderr,"CRD Reader: File corrupted\n");
                  fflush(stderr);
                  return -1;
                }
/*
                for( k = 0; k < 6; k++)
                {
                 strncpy(tcoor,&line[k*12],12);
                 tcoor[12] = '\0';
		 if( k == 0 || k == 3)
                 mol->x[(i*2)+k/3] = atof(tcoor);
	 	 else if( k == 1 || k == 4)
                 mol->y[(i*2)+k/3] = atof(tcoor);
                 else if( k == 2 || k == 5)
                 mol->z[(i*2)+k/3] = atof(tcoor);
                }
*/

                for( k = 0; k < 6; k++)
                {
                 index = floor((float) ((i*6)+k)/3.0f);
                 if( topo->excluded_list[index])
                     continue;
                 strncpy(tcoor,&line[k*12],12);
                 tcoor[12] = '\0';
                 r = ((i*6)+k) %  3;
                 if( r == 0)
                   mol->x[index] = atof(tcoor);
                 else if( r == 1)
                   mol->y[index] = atof(tcoor);
                 else if( r == 2)
                   mol->z[index] = atof(tcoor);
                }




        }
        fgets(line, 1024, input);  
        if( (j = (natoms*3) % 6) == 0){
                j = 6; 
        }

        if( line == NULL || strlen(line) < (j*12))
        {
           fprintf(stderr,"CRD Reader: File corrupted\n");
           fflush(stderr);
           return -1;
        }

 /*       for( k = 0; k < j; k++)
        {
                strncpy(tcoor,&line[k*12],12);
                tcoor[12] = '\0';
                 if( k == 0 || k == 3)
                 mol->x[(i*2)+k/3] = atof(tcoor);
                 else if( k == 1 || k == 4)
                 mol->y[(i*2)+k/3] = atof(tcoor);
                 else if( k == 2 || k == 5)
                 mol->z[(i*2)+k/3] = atof(tcoor);
        }*/

        for( k = 0; k < j; k++)
        {
                 index = floor((float) ((i*6)+k)/3.0f);
                 if( topo->excluded_list[index])
                     continue;

                 strncpy(tcoor,&line[k*12],12);
                 tcoor[8] = '\0';
                 r = ((i*6)+k) %  3;
                 if( r == 0)
                   mol->x[index] = atof(tcoor);
                 else if( r == 1)
                   mol->y[index] = atof(tcoor);
                 else if( r == 2)
                   mol->z[index] = atof(tcoor);

        }

	

	fclose(input);
	free(line);
	return 0;
}

float get_pdb_vdw_selected(MOL2 *mol, float ***byres, float *vdw_result, float *qq_result)
{
        int i = 0, j = 0, index = 0;
        float dx = 0.0, dy = 0.0f, dz = 0.0f;
        float dist = 0.0f, dist2 = 0.0f, distr = 0.0f;
        float r12 = 0.0f, r6 =0.0f, acoeff = 0.0, bcoeff = 0.0;
        float vdwE = 0.0f, qqE = 0.0f, tmpE = 0.0f, dScreening = 0.0f;
        float cte = 0.0f, lambda = 0.0f, epsilon = 0.0f, req = 0.0f;
        float alpha = 1.0367f;
        int first = 0, second = 0;

        float **res_energies = NULL;
        int current_res = 0;
        int res_index = 0;
        int last_atom = -1;

        cte = (78.39f - 1.0) / 2.0;
        lambda = alpha / (78.39 + 1.0);

        res_energies = *byres;

        for( i = 0; i < mol->n_atoms; i++)
        {
                res_index = 0;
                last_atom = -1;
                if( mol->backbone[i] == 0 || mol->exclude[i] == 1)
                        continue;

                for ( j = 0; j < mol->n_atoms; ++j)
                {

                        if( mol->vdw_selection[j] != 1 || mol->backbone[j] != 0 || mol->exclude[j] == 1)
                                continue;

                        res_index = mol->internal_res_num[j];
                                dx = (mol->x[i] - mol->x[j]);
                                dy = (mol->y[i] - mol->y[j]);
                                dz = (mol->z[i] - mol->z[j]);

                                dist = dx * dx + dy * dy + dz * dz;
                                dist2 = dist;
                                dist = sqrt(dist);
                                distr = 1 / dist;

                                dScreening = ((78.39 + 1.0) / (1.0 + cte * exp(-lambda * (78.39 + 1.0) * dist))) - 1;
                                tmpE = 332.0 * mol->pcharges[i]*mol->pcharges[j]*distr/dScreening;

                                qqE += tmpE;

                                res_energies[res_index][0] += tmpE;

                                epsilon = sqrt(mol->vdw_parm1[i] * mol->vdw_parm1[j]);
                                req = mol->vdw_parm2[i] + mol->vdw_parm2[j];
                                req = req * req * req * req * req * req; /* r0 ^6 */
                                acoeff = epsilon * (req*req);
                                bcoeff = 2.0*epsilon*req;

                                r6 = (dist2 * dist2 * dist2);
                                r12 = r6 * r6;

                                tmpE = (acoeff / r12) - (bcoeff/r6);

                                vdwE = vdwE + tmpE;

                                res_energies[res_index][1] += tmpE;

                }
        }

        *vdw_result = vdwE;
        *qq_result = qqE;
        return (vdwE+qqE);

}

float get_amber_vdw_selected( MOL2 *mol, TOP *topo, float ***byres, float *vdw_result, float *qq_result)
{
	int i = 0, j = 0, index = 0;
	float dx = 0.0, dy = 0.0f, dz = 0.0f;
	float dist = 0.0f, dist2 = 0.0f, distr = 0.0f;
	float r12 = 0.0f, r6 =0.0f;
	float vdwE = 0.0f, qqE = 0.0f, tmpE = 0.0f, dScreening = 0.0f;
	float cte = 0.0f, lambda = 0.0f;
        float alpha = 1.0367f;
	int first = 0.0f, second = 0.0f;

	float **res_energies = NULL;
	int current_res = 0;
	int res_index = 0;
	int last_atom = -1;

        cte = (78.39f - 1.0) / 2.0;
        lambda = alpha / (78.39 + 1.0);

	res_energies = *byres;

        for( i = 0; i < mol->n_atoms; i++)
        {
		res_index = 0;
		last_atom = -1;
		if( mol->backbone[i] == 0 || mol->exclude[i] == 1) 
			continue;

                for ( j = 0; j < mol->n_atoms; ++j) 
		{

                        if( mol->vdw_selection[j] != 1 || mol->backbone[j] != 0 || mol->exclude[j] == 1)
                                continue;
				
			if( last_atom == -1)
			{
				last_atom = j;
				current_res = mol->internal_res_num[j];
			}else{
				if( current_res != mol->internal_res_num[j] )
				{
/*				   ++res_index;*/
				   res_index = mol->internal_res_num[j]-1;
				   current_res = mol->internal_res_num[j];
				}
			}
                                dx = (mol->x[i] - mol->x[j]);
                                dy = (mol->y[i] - mol->y[j]);
                                dz = (mol->z[i] - mol->z[j]);

                                dist = dx * dx + dy * dy + dz * dz;
                                dist2 = dist;
                                dist = sqrt(dist);
                                distr = 1 / dist;

/*                                tmpE = (332.0 * lig->pcharges[i] * prot->pcharges[j]) *  distr * distr ;
 *                                                                vdwE += tmpE;*/

                                dScreening = ((78.39 + 1.0) / (1.0 + cte * exp(-lambda * (78.39 + 1.0) * dist))) - 1;
                                tmpE = 332.0 * mol->pcharges[i]*mol->pcharges[j]*distr/dScreening;

                                qqE += tmpE;
				/*
				fprintf(stderr,"Res Idx: %i\n",res_index);
				fflush(stderr);*/
                                res_energies[res_index][0] += tmpE;


				first = topo->atom_type_index[i];
				second = topo->atom_type_index[j];
/*				if( topo->atom_type_index[i] > topo->atom_type_index[j])
				{
	                                first = topo->atom_type_index[j];
	                                second = topo->atom_type_index[i];
				}*/

				index = topo->nb_index[(topo->ntypes * (first-1) + second)-1]-1;
				r6 = dist2 * dist2 * dist2;
				r12 = r6 * r6;
				tmpE = (topo->lj_A[index]/r12) - (topo->lj_B[index]/r6);

                                vdwE = vdwE + tmpE;

				res_energies[res_index][1] += tmpE;

                }
        }

/*	fprintf(stderr,"vdw: %f. qqE: %f\n",vdwE,qqE);*/
	*vdw_result = vdwE;
	*qq_result = qqE;
	return (vdwE+qqE);
}

float get_amber_vdw( MOL2 *mol, TOP *topo)
{
	int i = 0, j = 0, index = 0;
	float dx = 0.0, dy = 0.0f, dz = 0.0f;
	float dist = 0.0f, dist2 = 0.0f, distr = 0.0f;
	float r12 = 0.0f, r6 =0.0f;
	float vdwE = 0.0f, qqE = 0.0f, tmpE = 0.0f, dScreening = 0.0f;
	float cte = 0.0f, lambda = 0.0f;
        float alpha = 1.0367f;
	int first = 0.0f, second = 0.0f;

        cte = (78.39f - 1.0) / 2.0;
        lambda = alpha / (78.39 + 1.0);


        for( i = 0; i < mol->n_atoms; i++)
        {
                for ( j = i+1; j < mol->n_atoms; ++j) 
		{

                                dx = (mol->x[i] - mol->x[j]);
                                dy = (mol->y[i] - mol->y[j]);
                                dz = (mol->z[i] - mol->z[j]);

                                dist = dx * dx + dy * dy + dz * dz;
                                dist2 = dist;
                                dist = sqrt(dist);
                                distr = 1 / dist;

                                tmpE = (332.0 * mol->pcharges[i] * mol->pcharges[j]) *  distr * distr;
                                qqE += tmpE;

				first = topo->atom_type_index[i];
				second = topo->atom_type_index[j];

				index = topo->nb_index[(topo->ntypes * (first-1) + second)-1]-1;
				r6 = dist2 * dist2 * dist2;
				r12 = r6 * r6;
				tmpE = (topo->lj_A[index]/r12) - (topo->lj_B[index]/r6);

                                vdwE = vdwE + tmpE;

                }
        }

	fprintf(stderr,"vdw: %f. qqE: %f\n",vdwE,qqE);

	return (vdwE+qqE);
}

void get_next_mdcrd_name(char *base_mdcrd, char **mymdcrd_file, int current_mdcrd)
{
	char *mdcrd_file = NULL;
		mdcrd_file = (char *) calloc(sizeof(char), (strlen(base_mdcrd)+8));
	sprintf(mdcrd_file,"%s%i.mdcrd",base_mdcrd,current_mdcrd);
        *mymdcrd_file = mdcrd_file;


	return;
}

/* We have to call this routine before reading the next MDCRD file of a series */
void reset_mdcrd(TOP **mytopo)
{
	TOP *topo = NULL;

        topo = *mytopo;

	topo->frames = 0;
	topo->box = 0;
        if( topo->mdcrd_index != NULL)
        {
                free(topo->mdcrd_index);
        }

	topo->mdcrd_index = NULL;
	*mytopo = topo;
	return;
}

int load_mdcrd(MOL2 **mymol, char *fname, TOP **mytopo, int frame, int verbose)
{

        MOL2 *mol = NULL;
	TOP *topo = NULL;
        FILE *input = NULL;
        char *line = NULL;
        int natoms = 0;
        int i = 0, j = 0, k = 0, r = 0;
        char tcoor[13];
        int index = 0;
	long position = 0L;
	int total_lines = 0;
	int lines_per_frame = 0;

        mol = *mymol;
	topo = *mytopo;

        if( (input = fopen(fname,"rb")) == NULL)
        {
                fprintf(stderr,"Cannot open the multi-crd file\n");
                fflush(stderr);
                return -1;
        }

        if( (line = (char *) calloc(sizeof(char),1024)) == NULL)
        {
                fprintf(stderr,"Cannot allocate buffer\n");
                fflush(stderr);
                return 2;
        }

        fgets(line, 1023, input); /* title */

	natoms = topo->n_atoms;
	position = ftell(input);

        j = (int) ceil((float) natoms*3.0f/10.0f);

	/* Set frames to 0 to read another mdcrd in the same structure */
	if( topo->frames == 0 && topo->mdcrd_index != NULL)
	{
		free(topo->mdcrd_index);
		topo->box = 0;
	}

        if( topo->mdcrd_index == NULL)
	{

		if( verbose )
		{
			fprintf(stderr,"Analyzing mdcrd file...");
			fflush(stderr);
		}

	        while( fgets(line, 1024, input) )
		{
			++total_lines;
		}

		if( verbose )
		{
			fprintf(stderr," done\n");
	                fflush(stderr);
		}

		if( total_lines % j == 0)
		{
			topo->frames = total_lines/j;

	                if( verbose )
				fprintf(stderr,"MDCRD contains: %i frames\n",topo->frames);

			topo->box = 0;
		}else if( total_lines % (j+1) == 0){
			topo->frames = total_lines/(j+1);

			if( verbose)
		                fprintf(stderr,"MDCRD-boxed contains: %i frames\n",topo->frames);

			topo->box = 1;
		}else{
			fprintf(stderr,"MDCRD is corrupted\n");
			fflush(stderr);
                        fclose(input);
                        free(line);
			return -1;
		}
	}
	

        if( topo->mdcrd_index == NULL && topo->frames != 0)
        {

		if( verbose )
		{
	                fprintf(stderr,"No index - Please wait - Indexing...");
	                fflush(stderr);
		}

                topo->mdcrd_index = (long *) calloc(sizeof(long),topo->frames);
		if( topo->box )
			lines_per_frame = (j+1);
		else
	                lines_per_frame = (j);

		topo->mdcrd_index[0] = position;
		rewind(input);
                fgets(line, 1024, input);

/*		topo->mdcrd_index[0] = ftell(input);*/
		for( i = 0; i < topo->frames-1; i++)
		{
			for( k = 0; k < lines_per_frame; k++)
			{
				fgets(line, 1024, input);
			}
			topo->mdcrd_index[i+1] = ftell(input);
		}

		if( verbose )
		{
	                fprintf(stderr," done\n");
	                fflush(stderr);
		}

	}
        k = 0;
        rewind(input);
        fseek(input,topo->mdcrd_index[frame],SEEK_SET);

        for( i = 0; i < j-1; i++) /* Till the last line */
        {
                fgets(line, 1024, input);
		if( line == NULL || strlen(line) < (10*8))
		{
		  fprintf(stderr,"CRD Reader: Initial Block file corrupted: %s\n",line);
		  fflush(stderr);	
		  fclose(input);
		  free(line);
		  return -1;
		}
                for( k = 0; k < 10; k++)
                {
		 index = floor((float) ((i*10)+k)/3.0f);
		 if( topo->excluded_list[index])
                     continue;
                 strncpy(tcoor,&line[k*8],8);
                 tcoor[8] = '\0';
		 r = ((i*10)+k) %  3;
		 if( r == 0)
                   mol->x[index] = atof(tcoor);
		 else if( r == 1)
                   mol->y[index] = atof(tcoor);
		 else if( r == 2)
                   mol->z[index] = atof(tcoor);
                }

        }
        fgets(line, 1024, input);  /* Last line */
        if( (j = (natoms*3) % 10) == 0){
                j = 10; /* Complete then */
        }

        if( line == NULL || strlen(line) < (j*8))
        {
           fprintf(stderr,"CRD Reader: Last block line corrupted\n");
           fflush(stderr);
           fclose(input);
           free(line);
           return -1;
        }


        for( k = 0; k < j; k++)
        {
                 index = floor((float) ((i*10)+k)/3.0f);
                 if( topo->excluded_list[index])
                     continue;

                 strncpy(tcoor,&line[k*8],8);
                 tcoor[8] = '\0';
                 r = ((i*10)+k) %  3;
                 if( r == 0)
                   mol->x[index] = atof(tcoor);
                 else if( r == 1)
                   mol->y[index] = atof(tcoor);
                 else if( r == 2)
                   mol->z[index] = atof(tcoor);

        }


        free(line);
	fclose(input);
        return 0;
}


int split_amber_mol(MOL2 *mol, MOL2 **mylig, MOL2 **myprot)
{
	int i = 0, j = 0, k = 0;
	int lig_atoms = 0, prot_atoms = 0;

	MOL2 *lig = NULL, *prot = NULL;

	if( (lig = (MOL2 *) calloc(sizeof(MOL2),1)) == NULL)
	{
		fprintf(stderr,"Error. Cannot allocate memory.\n");
		fflush(stderr);
		return -2;
	}

        if( (prot = (MOL2 *) calloc(sizeof(MOL2),1)) == NULL)
        {
                fprintf(stderr,"Error. Cannot allocate memory.\n");
                fflush(stderr);
		free(lig);
                return -2;
        }

	for( i = 0; i < mol->n_atoms; i++)
	{
		if( mol->backbone[i] != 0)
			++lig_atoms;
		else
			++prot_atoms;
	}


        lig->n_atoms = lig_atoms;
        lig->x = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->y = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->z = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->pcharges = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->atoms = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->ringer = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->aromatic = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->bond_dist = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->grads_X = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->grads_Y = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->grads_Z = (float*)calloc(sizeof(float), lig->n_atoms);
        lig->backbone = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->exclude = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->selection = (int*)calloc(sizeof(int), lig->n_atoms + 10);
        lig->fragment_flag = (int*)calloc(sizeof(int), lig->n_atoms);
        lig->conformers = (CONFORMER*)calloc(sizeof(CONFORMER), 0);
        lig->n_fragments = 0;
        lig->gaff_types = (int*)calloc(sizeof(int), lig->n_atoms + 1);
        lig->vdw_parm1 = (float*)calloc(sizeof(float), lig->n_atoms + 1);
        lig->vdw_parm2 = (float*)calloc(sizeof(float), lig->n_atoms + 1);
        lig->ism_types = (int *)calloc(sizeof(int), lig->n_atoms + 1);
        lig->ism_selection = (int *)calloc(sizeof(int), lig->n_atoms + 1);
        lig->vdw_selection = (int *)calloc(sizeof(int), lig->n_atoms + 1);
        lig->res_num = (int *)calloc(sizeof(int), lig->n_atoms + 1);
        lig->internal_res_num = (int *)calloc(sizeof(int), lig->n_atoms + 1);


        prot->n_atoms = prot_atoms;
        prot->x = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->y = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->z = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->pcharges = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->atoms = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->ringer = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->aromatic = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->bond_dist = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->grads_X = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->grads_Y = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->grads_Z = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->backbone = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->exclude = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->selection = (int*)calloc(sizeof(int), prot->n_atoms + 10);
        prot->fragment_flag = (int*)calloc(sizeof(int), prot->n_atoms);
        prot->conformers = (CONFORMER*)calloc(sizeof(CONFORMER), 0);
        prot->n_fragments = 0;
        prot->gaff_types = (int*)calloc(sizeof(int), prot->n_atoms + 1);
        prot->vdw_parm1 = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->vdw_parm2 = (float*)calloc(sizeof(float), prot->n_atoms);
        prot->ism_types = (int *)calloc(sizeof(int), prot->n_atoms + 1);
        prot->ism_selection = (int *)calloc(sizeof(int), prot->n_atoms + 1);
        prot->vdw_selection = (int *)calloc(sizeof(int), prot->n_atoms + 1);
        prot->res_num = (int *)calloc(sizeof(int), prot->n_atoms + 1);
        prot->internal_res_num = (int *)calloc(sizeof(int), prot->n_atoms + 1);


	/* Initial coordinates */

	j = k = 0;
	for( i = 0; i < mol->n_atoms; i++)
	{
		if( mol->backbone[i] != 0 )
		{
			lig->x[j] = mol->x[i];
			lig->y[j] = mol->y[i];
			lig->z[j] = mol->z[i];
			lig->atoms[j] = mol->atoms[i];
			lig->pcharges[j] = mol->pcharges[i];
            lig->res_num[j] = mol->res_num[i];
            lig->internal_res_num[j] = mol->internal_res_num[i];
            lig->vdw_selection[k] = mol->vdw_selection[i];

			++j;
		}else{
            prot->x[k] = mol->x[i];
            prot->y[k] = mol->y[i];
            prot->z[k] = mol->z[i];
            prot->vdw_selection[k] = mol->vdw_selection[i];
			prot->atoms[k] = mol->atoms[i];
			prot->pcharges[k] = mol->pcharges[i];
			prot->res_num[k] = mol->res_num[i];
			prot->internal_res_num[k] = mol->internal_res_num[i];
			++k;
		}
	}

	mol_percieve(&lig);
	mol_percieve(&prot);

	*mylig = lig;
	*myprot = prot;

	return 0;
}

/* Only copy crds - do not initialize */
int split_amber_mol_crd(MOL2 *mol, MOL2 **mylig, MOL2 **myprot) 
{
	MOL2 *lig = NULL, *prot = NULL;
	int i = 0, j = 0, k = 0;

	lig = *mylig;
	prot = *myprot;

        j = k = 0;
        for( i = 0; i < mol->n_atoms; i++)
        {
                if( mol->backbone[i] != 0 )
                {
                        lig->x[j] = mol->x[i];
                        lig->y[j] = mol->y[i];
                        lig->z[j] = mol->z[i];
                        lig->vdw_selection[k] = mol->vdw_selection[i];
                        ++j;
                }else{
                        prot->x[k] = mol->x[i];
                        prot->y[k] = mol->y[i];
                        prot->z[k] = mol->z[i];
                        prot->vdw_selection[k] = mol->vdw_selection[i];
                        ++k;
                }
        }

	return 0;
}

void clean_topology(TOP **mytop)
{
	TOP *top = NULL;

	top = *mytop;


        free(top->nb_index );
        free(top->lj_A );
        free(top->lj_B );
        free(top->excluded_list );
	free(top->atom_type_index);

	if( top->mdcrd_index != NULL)
	  free(top->mdcrd_index);


	return;
}

void clean_amber_split_mol(MOL2 **mymol)
{

	MOL2 *lig = NULL;
	
	lig = *mymol;

        free(lig->x );
        free(lig->y );
        free(lig->z );
        free(lig->pcharges );
        free(lig->atoms );
        free(lig->ringer );
        free(lig->aromatic );
        free(lig->bond_dist );
        free(lig->grads_X );
        free(lig->grads_Y );
        free(lig->grads_Z );
        free(lig->backbone );
        free(lig->exclude );
        free(lig->selection );
        free(lig->fragment_flag );
        free(lig->conformers );
        free(lig->n_fragments );
        free(lig->gaff_types );
        free(lig->ism_types );
        free(lig->ism_selection );
        free(lig->vdw_selection );
        free(lig->res_num );
        free(lig->internal_res_num );

        free(lig->bonds);
        free(lig->bond_a1 );
        free(lig->bond_a2);

	return;

}

void clean_amber_mol(MOL2 **mymol)
{
	int i = 0;
	MOL2 *mol = NULL;

	mol = *mymol;


        free(mol->x );
        free(mol->y );
        free(mol->z );
        free(mol->pcharges );
        free(mol->radius );
        free(mol->atoms );
        free(mol->ringer );
        free(mol->aromatic );
        free(mol->bond_dist );
        free(mol->grads_X );
        free(mol->grads_Y );
        free(mol->grads_Z );
        free(mol->backbone );
        free(mol->exclude );
        free(mol->selection );
        free(mol->fragment_flag );
        free(mol->n_fragments );

        free(mol->res_num );
        free(mol->internal_res_num );
        free(mol->res_type );

        for( i = 0; i < mol->n_atoms; i++)
        {
                free(mol->res_names[i] );
                free(mol->atom_names[i] );
        }
	free(mol->res_names);
	free(mol->atom_names);
        free(mol->gaff_types );
        free(mol->ism_types );
        free(mol->hyde_types );

        free(mol->ism_selection );
        free(mol->hyde_selection );

        free(mol->vdw_selection );
	mol->nconformers = 1;
        for ( i = 0; i < mol->nconformers; ++i) {
                free(mol->conformers[i].x );
                free(mol->conformers[i].y );
                free(mol->conformers[i].z );
                free(mol->conformers[i].pcharges );
                free(mol->conformers[i].radius );
        }

	free(mol->conformers);

	return;
}

double quasiharmonic_entropy(double *eigenvalues, int n)
{
	int i = 0;
        double hwkT, w, dS, entropy = 0;
        double hbar, lambda;
	float temp = 298.0f;

    hbar = 6.62606957e-34/(2*3.14159265358979323846);
    for (i = 6; i < n; i++)
    {
        if (eigenvalues[i] > 0)
        {
            lambda = eigenvalues[i]*1.660538921e-27;
            w      = sqrt(1.3806488e-23*temp/lambda)/1e-10;
            hwkT   = (hbar*w)/(1.3806488e-23*temp);
            dS     = (hwkT/(exp(hwkT)-1.0) - log(1.0-exp(-hwkT)));
            entropy     += dS;
                fprintf(stdout, "i = %5d eig = %g w = %10g lam = %10g hwkT = %10g dS = %10g\n",
                        i, eigenvalues[i], w, lambda, hwkT, dS);
		fflush(stdout);
        }
    }

	return entropy*1.3806488e-23*6.02214129e23;
}

int update_split_coordinates(MOL2 *mol, MOL2 **mylig, MOL2 **myprot)
{
	int i = 0, j = 0, k = 0;
	int lig_atoms = 0, prot_atoms = 0;

	MOL2 *lig = NULL, *prot = NULL;

	lig = *mylig;
	prot = *myprot;

	for( i = 0; i < mol->n_atoms; i++)
	{
		if( mol->backbone[i] != 0)
			++lig_atoms;
		else
			++prot_atoms;
	}

	/* Initial coordinates */

	j = k = 0;
	for( i = 0; i < mol->n_atoms; i++)
	{
		if( mol->backbone[i] != 0 )
		{
			lig->x[j] = mol->x[i];
			lig->y[j] = mol->y[i];
			lig->z[j] = mol->z[i];
			lig->atoms[j] = mol->atoms[i];
			lig->pcharges[j] = mol->pcharges[i];
                        lig->res_num[j] = mol->res_num[i];
                        lig->internal_res_num[j] = mol->internal_res_num[i];
                        lig->vdw_selection[k] = mol->vdw_selection[i];

			++j;
		}else{
                        prot->x[k] = mol->x[i];
                        prot->y[k] = mol->y[i];
                        prot->z[k] = mol->z[i];
                        prot->vdw_selection[k] = mol->vdw_selection[i];
			prot->atoms[k] = mol->atoms[i];
			prot->pcharges[k] = mol->pcharges[i];
			prot->res_num[k] = mol->res_num[i];
			prot->internal_res_num[k] = mol->internal_res_num[i];
			++k;
		}
	}

	*mylig = lig;
	*myprot = prot;

	return 0;
}
