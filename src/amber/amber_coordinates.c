/**
 *
 *
 *
 *      @file amber_coordinates.c
 *      @brief Read a AMBER coordinate files to a MOL2 structure
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


