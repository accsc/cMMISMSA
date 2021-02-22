


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
                                {4.18 - 3.1, 12.92, 1.39, 28.49}, /* O.3 */
                                {14.18, 12.92, 1.39, 28.49}, /* O.2 */
                                {15.25, 13.79, 0.47, 31.33}, /* o.co2 */
                                {12.36, 13.85, 2.31, 30.82}, /* F */
                                {9.38 + 1.0, 9.69, 1.35, 22.04}, /* Cl */
                                {10.08 + 0.8, 8.47, 1.16, 19.71}, /* Br */
                                {9.90 + 1.0, 7.96, 0.96, 18.82}, /* I */
                                {10.13 + 0.5, 9.13, 1.38, 20.65}, /* S.3, S.2, S.o2, P.3 */
                               };


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

float assign_formal_charge(MOL2 *mol, int atom_type, int gaff_type)
{
    float valence = 0.0f, nb_term = 0.0f, bond_order = 0.0f, formal_charge = 0.0f;

    valence = assign_valence(atom_type);
    nb_term = assign_nonbonded(atom_type, gaff_type);
    bond_order = assign_nonbonded(mol, gaff_type);

    formal_charge = valence - nb_term - bond_order;

    return formal_charge;
}


void assign_gasteiger(MOL2 **mymol)
{
    int i = 0;
    MOL2 *mol = NULL;

    mol = *mymol;
    fprintf(stderr,"Charges:\n");
    for( i = 0; i < mol->n_atoms; i++)
    {
        fprintf(stderr,"%i %s %f\n",i, let[mol->gaff_types[i]-1], assign_formal_charge(mol, mol->atoms[i], mol->gaff_types[i]));
    }

    *mymol = mol;
}


float assign_bond_order(MOL2 *mol, int at_idx)
{
    int i = 0, order = 0, aro = 0;

    for(i = 0; i < mol->n_bonds; i++)
    {
        if( mol->bond_a1[i] == (at_idx+1) || mol->bond_a2[i] == (at_idx+1))
        {
            if( mol->aromatic[mol->bond_a1[i]-1] == 1 || mol->aromatic[mol->bond_a2[i]-1] == 1)
            {
                aro++;
                order++;
            }else{
                order += mol->bonds[i];
            }
        }
        
    }

   if (aro > 0)
            order = order + aro + 1;

   return (float) order;

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

    if( atom_type == 1 && fabs(charge-1.0) < 1E-5) /* H */
        chi = 1.0;
    else{
        chi = poly_terms[0] + poly_terms[1]*charge + poly_terms[2] * charge * charge + poly_terms[3]*charge*charge*charge;
    }

    return chi;
}

