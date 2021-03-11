/**
 *
 *      @file waste.c
 *      @brief Memory clean up for ligand structures.
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @date 01/10/2011
 *
 * 	This program is free software; you can redistribute it and/or modify
 * 	it under the terms of the GNU General Public License as published by
 * 	the Free Software Foundation version 2 of the License.
 *
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 	GNU General Public License for more details.
 *
 */


/**
 *
 *	@brief Clean a MOL2 structure allocated by a reader
 *	@param mymols MOL2 structure to clean
 *	@return 0 if success
 *
 */
int cleanup(MOL2 **mymols)
{
	MOL2 *mols = NULL;
	int i = 0, j = 0;

	mols = *mymols;
	free(mols->x);
	free(mols->y);
	free(mols->z);
	free(mols->pcharges);
	free(mols->radius);
	free(mols->atoms);
	free(mols->ringer);
	free(mols->aromatic);
	free(mols->exclude);
	free(mols->bonds);
	free(mols->bond_a1);
	free(mols->bond_a2);
	free(mols->bond_dist);
	free(mols->grads_X);
	free(mols->grads_Y);
	free(mols->grads_Z);
	free(mols->ik);
	free(mols->jk);
	free(mols->kk);
	free(mols->lk);
	free(mols->ip);
	free(mols->jp);
	free(mols->kp);
	free(mols->lp);
	free(mols->tp);
	free(mols->ia);
	free(mols->ja);
	free(mols->ka);
	free(mols->gaff_types);
	free(mols->vdw_parm1);
	free(mols->vdw_parm2);
        free(mols->ism_types);
        free(mols->hyde_types);
        free(mols->ism_selection);
        free(mols->hyde_selection);
	free(mols->backbone);
	free(mols->fragment_flag);
	free(mols->selection);
	free(mols->vdw_type);
	free(mols->vdwpairs_a);
	free(mols->vdwpairs_b);
	free(mols->internal_res_num);
	free(mols->res_num);
	free(mols->res_type);
        for( j = 0; j < mols->n_atoms; j++)
        {
		free(mols->res_names[j]);
		free(mols->atom_names[j]);
        }
	free(mols->res_names);
        free(mols->atom_names);

        for ( i = 0; i < mols->nconformers; ++i) 
	{
                free(mols->conformers[i].x);
                free(mols->conformers[i].y);
                free(mols->conformers[i].z);
                free(mols->conformers[i].pcharges);
                free(mols->conformers[i].radius);
        }

	free(mols->conformers);
        free(mols->vdw_selection);


	free(mols);
	return 0;
}
