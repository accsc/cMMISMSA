/**
 *
 *      @file assign_gasteiger.c
 *      @brief Calculate PEOE charges using the Gasteiger method
 *
 *      @author Alvaro Cortes Cabrera <alvarocortes@gmail.com>
 *      @date 2021/02
 *     
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation version 2 of the License.
 *     
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 */



float polynomial_terms[20][4] = {{7.17, 6.24, -0.56, 12.85}, /* H */
                                {7.98, 9.18, 1.88, 19.04}, /* C.3 C.cat */
                                {8.79 + 0.5, 9.32, 1.51, 19.62}, /* C.2 */
                                {7.98 + 0.55, 9.18, 1.88, 19.04}, /* C.ar */
                                {10.39, 9.45, 0.73, 20.57}, /* C.1 */
                                {11.54 + 6.0, 10.28, 1.36, 28.00}, /* N.3 N.4 */
                                {12.87 - 1.29, 11.15, 0.85, 24.87}, /* N.ar */
                                {12.87, 11.15, 0.85, 24.87}, /* N.2 */
                                {12.87 + 0.5, 11.15, 0.85, 24.87}, /* N.pl3 */
                                {12.87 + 3.5, 11.15, 0.85, 24.87}, /* N.am */
                                {15.68, 11.70, -0.27, 27.11}, /* N.1 */
                                {14.18 + 0.8, 12.92, 1.39, 28.49}, /* O.oh */
                                {14.18 - 3.1, 12.92, 1.39, 28.49}, /* Not valid, used to be O.3 but all O.3 are O.oh */
                                {14.18, 12.92, 1.39, 28.49}, /* O.2 */
                                {15.25, 13.79, 0.47, 31.33}, /* o.co2 */ /* Pending! */
                                {12.36, 13.85, 2.31, 30.82}, /* F */
                                {9.38 + 1.0, 9.69, 1.35, 22.04}, /* Cl */
                                {10.08 + 0.8, 8.47, 1.16, 19.71}, /* Br */
                                {9.90 + 1.0, 7.96, 0.96, 18.82}, /* I */
                                {10.13 + 0.5, 9.13, 1.38, 20.65}, /* S.3, S.2, S.o2, P.3 */
                               };

int get_polynomial_index(MOL2 *mol, int atom)
{
    int idx = -1;

    if( mol->atoms[atom] == 4) /* H */
        idx = 0;
    else if( mol->atoms[atom] == 1) /* C */
    {
            if( mol->aromatic[atom] != 0)
                idx = 3;
            else if( mol->gaff_types[atom] == C1)
                idx = 4;
            else if( mol->gaff_types[atom] == C || mol->gaff_types[atom] == C2 || (mol->gaff_types[atom] >= CC && mol->gaff_types[atom] <= CF))
                idx = 2;
            else
                idx = 1;
    }else if( mol->atoms[atom] == 3) /* N */
    {
            if( mol->aromatic[atom] != 0)
                idx = 6;
            else if( mol->gaff_types[atom] == N1)
                idx = 10;
            else if( mol->gaff_types[atom] == N || mol->gaff_types[atom] == NH)
                idx = 9;
            else if( mol->gaff_types[atom] == N2 || mol->gaff_types[atom] == NC || mol->gaff_types[atom] == NE)
                idx = 7;
            else if( mol->gaff_types[atom] == N3 && get_number_any_bond_p(mol, atom+1, 0) == 3) /* N.pl3 */
                idx = 8;
            else
                idx = 5;
        
    }else if( mol->atoms[atom] == 2) /* O */
    {
            if( mol->gaff_types[atom] == OH || mol->gaff_types[atom] == OW)
                idx = 11;
            else if( mol->gaff_types[atom] == O)
                idx = 13;
            else
                idx = 11;
    }else if( mol->atoms[atom] == 10) /* F */
        idx = 15;
    else if( mol->atoms[atom] == 9) /* Cl */
        idx = 16;
    else if( mol->atoms[atom] == 8) /* Br */
        idx = 17;
    else if( mol->atoms[atom] == 7) /* I */
        idx = 18;
    else if( mol->atoms[atom] == 5 || mol->atoms[atom] == 6)
        idx = 19;

    return idx;

}


float assign_nonbonded(int atom_type, int gaff_type)
{
    float nb_type = 0.0;
    if (atom_type >= 7 && atom_type <=10) /* Halogen */
        nb_type = 6;
    else if( atom_type == 3 && gaff_type != N4 && gaff_type != N) /* Missing N.pl3 in this list */
        nb_type = 2;
    else if (atom_type == 2) /* COO should be 4.5 */
        nb_type = 4;
    else if (atom_type == 6) {
        nb_type = 4;
        if( gaff_type == SX) /* S.o */
            nb_type =2;
        else if( gaff_type == SY) /* S.o2 */
            nb_type = 0;
    }

    return nb_type;

}

float assign_valence(int atom_type)
{

    float valence = 0.0;
    if( atom_type == 1)
        valence = 4;
    else if (atom_type == 2)
        valence = 6;
    else if (atom_type == 3)
        valence = 5;
    else if (atom_type == 4)
        valence = 1;
    else if (atom_type == 5)
        valence = 5;
    else if (atom_type == 6)
        valence = 6;
    else if (atom_type == 7)
        valence = 7;
    else if (atom_type == 8)
        valence = 7;
    else if (atom_type == 9)
        valence = 7;
    else if (atom_type == 10)
        valence = 7;
    return valence;
}


float assign_bond_order(MOL2 *mol, int at_idx)
{
    int i = 0, a1 = 0, a2 = 0;
    float order = 0, aro = 0;

    for(i = 0; i < mol->n_bonds; i++)
    {
        a1 = mol->bond_a1[i]-1;
        a2 = mol->bond_a2[i]-1;
        if( a1 == at_idx || a2 == at_idx)
        {
            if( mol->aromatic[a1] != 0 && mol->aromatic[a2] != 0)
            {
                aro = aro + 1.0;
                /*fprintf(stderr,"Aromatic flag: %s %s %i\n", let[mol->gaff_types[a1]-1],let[mol->gaff_types[a2]-1],mol->bonds[i]);*/
            }else{
                order += mol->bonds[i];
                /*fprintf(stderr,"%s %s %i\n", let[mol->gaff_types[a1]-1],let[mol->gaff_types[a2]-1],mol->bonds[i]);*/
            }
        }
        
    }

   if (aro > 0)
   {
            order = order + aro + 1.0;
   }

   return order;

}

float assign_formal_charge(MOL2 *mol, int i, int atom_type, int gaff_type)
{   
    float valence = 0.0f, nb_term = 0.0f, bond_order = 0.0f, formal_charge = 0.0f;
    
    valence = assign_valence(atom_type);
    nb_term = assign_nonbonded(atom_type, gaff_type);
    bond_order = assign_bond_order(mol, i);
    formal_charge = valence - nb_term - bond_order;

    if( gaff_type == N && bond_order == 3)
        formal_charge = 0;
    else if( gaff_type == NA && bond_order == 4)
        formal_charge = 0;
    else if( gaff_type == CA && bond_order == 5)
        formal_charge = 0;
    else if( gaff_type == C2 && bond_order == 5 && formal_charge == -1.0)
        formal_charge = 0;
    else if( gaff_type == N3 && bond_order == 4 && formal_charge == -1.0)
        formal_charge = 1.0;
    else if( gaff_type == O && bond_order == 1)
    {
        /*fprintf(stderr,"Prev %f\n",formal_charge);*/
        formal_charge = -0.5;
    }

    /*fprintf(stderr,"%f %f %f %f\n",valence,nb_term,bond_order,formal_charge);*/
    
    return formal_charge;
}

float electronegativity(float charge, float poly_terms[4], int atom_type)
{
    float chi = 0.0;

    if( charge > 1.1)
    {
        charge = 1.1;
    }else if( charge < -1.1){
        charge = -1.1;
    }

    if( atom_type == 4 && fabs(charge-1.0) < 1E-5) 
        chi = 20.02;
    else{
        chi = poly_terms[0] + (poly_terms[1]*charge) + (poly_terms[2] * charge * charge) + (poly_terms[3]*charge*charge*charge);
    }

   /*fprintf(stderr,"Electro %f %f %f %f %f %f\n",charge, poly_terms[0], poly_terms[1], poly_terms[2], poly_terms[3], chi);*/

    return chi;
}




void assign_gasteiger(MOL2 **mymol)
{
    int i = 0, j = 0, k = 0, *poly = NULL, neig = 0;
    MOL2 *mol = NULL;
    float *fcharges = NULL, *deltas = NULL, *ccharges = NULL, scaling_factor = 1.56, abs_chg = 0.0f ;
    float chi1 = 0.0f, chi2 = 0.0, diff_chi = 0.0, chi_norm = 0.0f;

    mol = *mymol;

    if( (fcharges = (float *) calloc(sizeof(float), mol->n_atoms+1)) == NULL)   
    {
        fprintf(stderr,"Memory error on PEOE charges assisgment\n");
        return;
    }

    if( (deltas = (float *) calloc(sizeof(float), mol->n_atoms+1)) == NULL)
    {
        fprintf(stderr,"Memory error on PEOE charges assisgment\n");
        free(fcharges);
        return;
    }

    if( (poly = (int *) calloc(sizeof(int), mol->n_atoms+1)) == NULL)
    {
        fprintf(stderr,"Memory error on PEOE charges assisgment\n");
        free(fcharges);
        free(deltas);
        return;
    }

    if( (ccharges = (float *) calloc(sizeof(float), mol->n_atoms+1)) == NULL)
    {
        fprintf(stderr,"Memory error on PEOE charges assisgment\n");
        free(fcharges);
        free(deltas);
        free(poly);
        return;
    }



/*    fprintf(stderr,"Charges:\n");*/
    for( i = 0; i < mol->n_atoms; i++)
    {
        /*fprintf(stderr,"%i %s %f\n",i, let[mol->gaff_types[i]-1], assign_formal_charge(mol, i, mol->atoms[i], mol->gaff_types[i]));*/
        fcharges[i] = assign_formal_charge(mol, i, mol->atoms[i], mol->gaff_types[i]) / scaling_factor;
        ccharges[i] = 0.0f;
        deltas[i] = 0.0f;
        abs_chg += fabs(fcharges[i]);
        poly[i] = get_polynomial_index(mol, i);
/*        if (poly[i] < 0)
            fprintf(stderr,"Error for %i\n",i);
        else
            fprintf(stderr,"FC: %f\n",fcharges[i]);*/

    }

/*    fprintf(stderr,"Abs %f\n",abs_chg);*/

    for( i = 0; i < 6; i++)
    {
        for( j = 0; j < mol->n_atoms; j++)
        {
            chi1 = electronegativity(ccharges[j], polynomial_terms[poly[j]], mol->atoms[j]);
            deltas[j] = 0.0f;
            for( k = 0; k < mol->n_bonds; k++)
            {
                neig = -1;
                if( mol->bond_a1[k]-1 == j)
                    neig = mol->bond_a2[k]-1;
                else if( mol->bond_a2[k]-1 == j)
                    neig = mol->bond_a1[k]-1;

                if( neig >= 0)
                {
                   /*fprintf(stderr,"Atom %i. Neig %i. %s %s \n",j,neig,let[mol->gaff_types[j]-1],let[mol->gaff_types[neig]-1]);*/
                   chi2 = electronegativity(ccharges[neig], polynomial_terms[poly[neig]], mol->atoms[neig]);
                   diff_chi = chi2 - chi1;
                   
                   if (chi2 > chi1)
                       chi_norm = electronegativity(1.0, polynomial_terms[poly[j]], mol->atoms[j]);
                   else
                       chi_norm = electronegativity(1.0, polynomial_terms[poly[neig]], mol->atoms[neig]);

                   deltas[j] += (diff_chi / chi_norm) * (powf(0.778, (i + 1)));
                   /*fprintf(stderr,"%i chi1: %f chi2: %f chi_norm: %f Adding: %f\n",j,chi1,chi2, chi_norm, (diff_chi / chi_norm) * (powf(0.778, (i + 1))));*/
                }
            }
            /*fprintf(stderr,"%i %f %f\n", j, chi1, deltas[j]);*/
        }


        for( j = 0; j < mol->n_atoms; j++)
        {
            if (abs_chg < 0.001)
                ccharges[j] += deltas[j];
            else{
                ccharges[j] += deltas[j] + (1.0/6.0) * fcharges[j];
            }

            /*fprintf(stderr,"Cycle %i: %i %f\n",i,j,deltas[j]);*/

        }


    }

    for( j = 0; j < mol->n_atoms; j++)
    {
        ccharges[j] *= scaling_factor;
        /*fprintf(stderr,"%i %f\n", j, ccharges[j]);*/
        mol->pcharges[j] = ccharges[j];
    }


    *mymol = mol;

    free(deltas);
    free(fcharges);
    free(ccharges);
    free(poly);
}



