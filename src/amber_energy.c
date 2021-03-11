/**
 *
 *
 *
 *      @file amber_energy.c
 *      @brief Calculate several energy terms from the AMBER force field
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
                fprintf(stderr,"Current charge: %i %f\n",i,mol->pcharges[i]);

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
                                /*fprintf(stderr,"Prot charge: %i %f\n",j, mol->pcharges[j]);*/

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


