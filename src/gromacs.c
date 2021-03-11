/**
 *
 *
 *
 *      @file gromacs.c
 *      @brief Read XTC files
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 11/03/2018
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

int load_xtc(MOL2 **mymol, char *fname, TOP **mytopo, int verbose_flag)
{
	XDRFILE *xdf = NULL;

        MOL2 *mol = NULL;
	TOP *topo = NULL;
        int natoms = 0;
	int step = 0;
	int result = 0;
	float time = 0.0f, prec = 0.0f;
	rvec *x = NULL;
	matrix box;
        int index = 0, i = 0;

        mol = *mymol;
	topo = *mytopo;

	xdf = topo->xdr;
	
	if( xdf == NULL)
	{
        	if( (xdf = xdrfile_open(fname,"r")) == NULL){
               	 fprintf(stderr,"Cannot open XTC file\n");
               	 fflush(stderr);
               	 return -1;
        	}
		
		result = xtc_header(xdf,&natoms,&step,&time,0);
		topo->n_atoms = natoms;
		fprintf(stderr,"XTC reader: Header reported %i atoms\n",natoms);
		fflush(stderr);
		x = calloc(topo->n_atoms,sizeof(*x));
		/*result = xtc_coord(xdf,&natoms,box,x,&prec,0);
		if( result != exdrOK)
		{
			fprintf(stderr,"XTC reader. Error reading first step of file\n");
			fflush(stderr);
			return -1;
		}	*/
		xdrfile_close(xdf);
		fprintf(stderr,"XTC reader: Re-opening xtc file after header extraction\n");
		fflush(stderr);
		xdf = xdrfile_open(fname,"r");
		topo->xdr = xdf;
		result = read_xtc(xdf,natoms,&step,&time,box,x,&prec);
		fprintf(stderr,"XTC reader. Atoms: %i Step: %i Time: %f. Prec: %f\n",topo->n_atoms,step,time,prec);
		fflush(stderr);
                xdrfile_close(xdf);
                xdf = xdrfile_open(fname,"r");
                topo->xdr = xdf;
		topo->frames = -1; /* No indexing support for XTC */
	}else{
		natoms = topo->n_atoms;
        	x = calloc(topo->n_atoms,sizeof(*x));
		result = read_xtc(xdf,natoms,&step,&time,box,x,&prec);
		if (exdrENDOFFILE == result)
		{
			free(x);
			if( verbose_flag)
			{
                        	fprintf(stderr,"XTC reader: End of file\n");
                        	fflush(stderr);
			}
                        return -1;
		}else{
                        if (exdrOK != result)
                        {
				free(x);
                                fprintf(stderr,"XTC reader. Error reading step of file\n");
                                fflush(stderr);
                                return -1;
                        }else{
				if( verbose_flag)
				{
					fprintf(stderr,"XTC reader. Atoms: %i Step: %i. Time: %f, Precision: %f\n",topo->n_atoms,step,time,prec);
					fflush(stderr);
				}
			}
		}
	}

        for( index = 0; index < topo->n_atoms; index++)
        {
                 if( topo->excluded_list[index])
                     continue;
                   mol->x[index] = x[index][0]*10.0; /* nm to A */
                   mol->y[index] = x[index][1]*10.0; /* nm to A */
                   mol->z[index] = x[index][2]*10.0; /* nm to A */
        }

	free(x);
        return 0;
}


