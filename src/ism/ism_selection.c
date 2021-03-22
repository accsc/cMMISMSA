/**
 *
 *      @file ism.c
 *      @brief Implicit solvation model atom selection outines
 *
 *      @author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *      @date 2021/03
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
 */


/**
 *	@brief Select only atoms in the grid for Implicit solvent calculation
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@param mymol Pointer to protein MOL2 structure
 *	@param max_grid Max grid coordinates
 *	@param min_grid Min grid coordinates
 *	@return 0 on success
 */
int ism_select_atoms(MOL2 **mymol, float *max_grid, float *min_grid)
{
	MOL2 *mol = NULL;
	int i = 0;
	mol = *mymol;


	for( i = 0; i < mol->n_atoms; i++)
	{
		mol->ism_selection[i] = 0;

	        if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
		{
			mol->ism_selection[i] = 1;
		}
	}


	*mymol = mol;
	return 0;
}

/**
 *	@brief Calculate the distance matrix for a molecule
 *	@param mymol Molecule MOL2 structure
 *	@param mymatrix Matrix
 *	@return 0 on success
 */
int ism_get_distance_matrix(MOL2 *mymol, float ***mymatrix)
{
	float **matrix = NULL;
	int i = 0, j = 0;
	float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;

	matrix = *mymatrix;


	for( i = 0; i < mymol->n_atoms; i++)
	{
		matrix[i][i] = 0.0f;
		for( j = i + 1; j < mymol->n_atoms; j++)
		{
			dx = mymol->x[i] - mymol->x[j];
                        dy = mymol->y[i] - mymol->y[j];
                        dz = mymol->z[i] - mymol->z[j];
			dist = ( dx * dx) + ( dy * dy ) + ( dz * dz);
			dist = sqrtf(dist);
			matrix[i][j] = dist;
			matrix[j][i] = dist;
		}
	}

	*mymatrix = matrix;
	return 0;
}


int ism_select_atoms_by_ligand_noconformer(MOL2 **mymol, MOL2 *lig)
{
        MOL2 *mol = NULL;
        int i = 0, j = 0;
        float max_grid[3], min_grid[3];
        mol = *mymol;

        max_grid[0] = max_grid[1] = max_grid[2] = -9999.9f;
        min_grid[0] = min_grid[1] = min_grid[2] = 9999.9f;

        for( j = 0; j < lig->n_atoms; j++)
        {
                if( lig->x[j] > max_grid[0])
                 max_grid[0] = lig->x[j];
                if( lig->y[j] > max_grid[1])
                 max_grid[1] = lig->y[j];
                if( lig->z[j] > max_grid[2])
                 max_grid[2] = lig->z[j];

                if( lig->x[j] < min_grid[0])
                 min_grid[0] = lig->x[j];
                if( lig->y[j] < min_grid[1])
                 min_grid[1] = lig->y[j];
                if( lig->z[j] < min_grid[2])
                 min_grid[2] = lig->z[j];
        }
        for( i = 0; i < 3 ; i++)
        {
                max_grid[i] += 10.0f;
                min_grid[i] -= 10.0f;
        }


        for( i = 0; i < mol->n_atoms; i++)
        {
		if( mol->vdw_selection[i] == 2)
		   continue;

                mol->ism_selection[i] = 0;
                mol->vdw_selection[i] = 0;

                if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
                {
                        mol->ism_selection[i] = 1;
                        mol->vdw_selection[i] = 1;

                }
        }


        *mymol = mol;
        return 0;

}


int ism_select_atoms_by_ligand(MOL2 **mymol, MOL2 *lig)
{
        MOL2 *mol = NULL;
        int i = 0, j = 0;
        float max_grid[3], min_grid[3];
        mol = *mymol;

        max_grid[0] = max_grid[1] = max_grid[2] = -9999.9f;
        min_grid[0] = min_grid[1] = min_grid[2] = 9999.9f;


        for( i = 0; i < lig->nconformers; i++)
        {
                set_default_conformer(&lig, i);
                for( j = 0; j < lig->n_atoms; j++)
                {
                        if( lig->x[j] > max_grid[0])
                         max_grid[0] = lig->x[j];
                        if( lig->y[j] > max_grid[1])
                         max_grid[1] = lig->y[j];
                        if( lig->z[j] > max_grid[2])
                         max_grid[2] = lig->z[j];

                        if( lig->x[j] < min_grid[0])
                         min_grid[0] = lig->x[j];
                        if( lig->y[j] < min_grid[1])
                         min_grid[1] = lig->y[j];
                        if( lig->z[j] < min_grid[2])
                         min_grid[2] = lig->z[j];
                }
        }

        for( i = 0; i < 3 ; i++)
        {
                max_grid[i] += 10.0f;
                min_grid[i] -= 10.0f;
        }


        for( i = 0; i < mol->n_atoms; i++)
        {
                mol->ism_selection[i] = 0;
                mol->vdw_selection[i] = 0;

                if( mol->x[i] >= min_grid[0] && mol->x[i] <= max_grid[0] &&
                    mol->y[i] >= min_grid[1] && mol->y[i] <= max_grid[1] &&
                    mol->z[i] >= min_grid[2] && mol->z[i] <= max_grid[2])
                {
                        mol->ism_selection[i] = 1;
                        mol->vdw_selection[i] = 1;

                }
        }


        *mymol = mol;
        return 0;

}

