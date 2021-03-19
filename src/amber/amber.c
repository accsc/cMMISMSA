/**
 *
 *
 *
 *      @file amber.c
 *      @brief Read a AMBER top file to a MOL2 structure
 *
 *      @author Alvaro Cortes Cabrera <alvarocortes@gmail.com>
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

/** XDRFILE lib for XTC file support */
#include <gromacs/xdrfile_xtc.h>
#include <gromacs/xdrfile.c>
#include <gromacs/xdrfile_xtc.c>

#include <amber/amber.h>
#include <amber/amber_topology.c>
#include <amber/amber_coordinates.c>
#include <amber/amber_energy.c>



int split_amber_mol(MOL2 **mymol, MOL2 **mylig, MOL2 **myprot)
{
	int i = 0, j = 0, k = 0;
	int lig_atoms = 0, prot_atoms = 0;

	MOL2 *lig = NULL, *prot = NULL, *mol = NULL;

    mol = *mymol;

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

    init_molecule (&lig, lig_atoms, 1);
    init_molecule (&prot, prot_atoms, 1);


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
            lig->vdw_selection[j] = mol->vdw_selection[i];

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
    check_assign_types(lig);
    check_assign_types(prot);

    assign_gasteiger(&lig);

    j = k = 0;
    for( i = 0; i < mol->n_atoms; i++)
    {
        if( mol->backbone[i] != 0 )
        {
            mol->pcharges[i] = lig->pcharges[j];
            ++j;
        }
    }

	*mylig = lig;
	*myprot = prot;
    *mymol = mol;

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
        free(lig->vdw_parm1 );
        free(lig->vdw_parm2 );
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

        free(mol->vdw_parm1 );
        free(mol->vdw_parm2 );


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
			/*lig->pcharges[j] = mol->pcharges[i];*/
                        lig->res_num[j] = mol->res_num[i];
                        lig->internal_res_num[j] = mol->internal_res_num[i];
                        lig->vdw_selection[j] = mol->vdw_selection[i];

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

char *ltrim(char *s)
{

    while(isspace(*s)) s++;
    return s;
}

char *rtrim(char *s)
{

    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}


char *trim(char *s)
{
    return rtrim(ltrim(s)); 
}


int detect_AMBER_HIS_type(MOL2 *mol, int res_number)
{
    int i = 0, type1 = 0, type2 = 0, res = 0;
    char *atom_name = NULL, *res_name = NULL;

    for(i = 0; i < mol->n_atoms; i++)
    {
        atom_name = mol->atom_names[i];
        res_name = mol->res_names[i];

        if (mol->internal_res_num[i] == res_number)
        {
            if ( (strcmp(res_name, "HIS") == 0) && (strcmp(trim(atom_name), "HD1") == 0))
                type1 = 1; /* HID */
            else if( (strcmp(res_name, "HIS") == 0) && (strcmp(trim(atom_name), "HE2") == 0))
                type2 = 1; /* HIE */
        }
    
    }

    if( type1 == 1 && type2 == 1)
        res = 3; /* HIP */
    else if(type1 == 1) 
        res = 1; /* HID */
    else if(type2 == 1)
        res = 2; /* HIE */

    /* res = 0 means no protons found */

    return res;
}


void AMBER_atom_typing (MOL2 ** mymol)
{
    int i = 0, j = 0;
    MOL2 *mols = NULL;
    char *atom_name = NULL, *res_name = NULL;

    mols = *mymol;

    for(i = 0; i < mols->n_atoms; i++)
    {
        atom_name = mols->atom_names[i];
        res_name = mols->res_names[i];

        for(j = 0; j < 2257; j++)
        {
           if ( (strcmp(res_name, amber_res_names[j]) == 0) && (strcmp(trim(atom_name), amber_atom_names[j]) == 0))
           {
                mols->vdw_parm1[i] = amber_vdw_epsilon[j];
                mols->vdw_parm2[i] = amber_vdw_r0[j];
                mols->pcharges[i] = amber_pcharges[j];
           }
        }
    }

    *mymol = mols;
}


