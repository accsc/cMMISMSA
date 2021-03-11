/**
 *
 *      @file amber_topology.c
 *      @brief Read a AMBER top file to a MOL2 structure
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 2021/02
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

