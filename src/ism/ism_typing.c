/**
 *
 *      @file ism_typing.c
 *      @brief Implicit solvation model atom typing routines
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
 *	@brief Atom typing for implicit solvation model
 *	@param mymol Pointer to molecule MOL2 structure
 *	@return 0 on success
 */
int ism_typing(MOL2 **mymol)
{
	MOL2 *mols = NULL;
        int vecinos[8];
        int vecinos2[8];
	int k = -1, k_heavy = -1, i = 0, j = 0;

	mols = *mymol;

	
        for (i = 0; i < mols->n_atoms; ++i)
               mols->ism_types[i] = -1;


        for (i = 0; i < mols->n_atoms; ++i) {

                k = 0;
		k_heavy = 0;
                for (j = 0; j < mols->n_bonds; ++j) {
			if ( mols->bond_a1[j] == (i + 1)) {
                                k++;
				if( mols->atoms[mols->bond_a2[j] - 1] != 4 )
					k_heavy++;
                        }else if ( mols->bond_a2[j] == (i + 1)) {
                                k++;
                                if( mols->atoms[mols->bond_a1[j] - 1] != 4 )
                                        k_heavy++;
                        }
                }


                if ( mols->atoms[i] == 4 ) {  /* H */
                	mols->ism_types[i] = 1;
                }else if ( mols->atoms[i] == 9 ) {  /* Cl */
                        mols->ism_types[i] = 2;
                }else if ( mols->atoms[i] == 105 ) {  /* Ca or C0 */
                        mols->ism_types[i] = 3;
                }else if ( mols->atoms[i] == 5 ) {  /* P */
                        mols->ism_types[i] = 4;
                }else if ( mols->atoms[i] == 3 ) {  /* N */
			if ( k >= 4)
			{
				if( k_heavy == 1)
				  mols->ism_types[i] = 7;
				else if ( k_heavy == 2)
                                  mols->ism_types[i] = 6;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 5;
				else
                                  mols->ism_types[i] = 12;

			}else if( k == 3){

                               if( k_heavy == 1)
                                  mols->ism_types[i] = 10;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 9;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 8;
				else
                                  mols->ism_types[i] = 10;

			}else if( k == 2){
                              if( k_heavy == 1)
                                  mols->ism_types[i] = 12;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 11;
                                else
                                  mols->ism_types[i] = 12;


			}else{
				mols->ism_types[i] = 12;
			}

		}else if ( mols->atoms[i] == 2 ) {  /* O */
                        if ( k >= 2)
                        {
                                if( k_heavy == 1)
                                  mols->ism_types[i] = 14;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 13;
                                else
                                  mols->ism_types[i] = 14;

                        }else if( k == 1){
                                mols->ism_types[i] = 15;


                        }else{
                                mols->ism_types[i] = 15;
                        }


		}else if ( mols->atoms[i] == 6 ) {  /* S */
	                if( k_heavy == 2)
	                        mols->ism_types[i] = 18;
			else if( k_heavy == 3 )
                                mols->ism_types[i] = 17;
                        else if( k_heavy == 4 )
                                mols->ism_types[i] = 16;
			else
                                mols->ism_types[i] = 16;
		}else if ( mols->atoms[i] == 1 ) {  /* C */

                        if ( k >= 4)
                        {
                                if( k_heavy == 1)
                                  mols->ism_types[i] = 22;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 21;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 20;
                                else if ( k_heavy >= 4)
                                  mols->ism_types[i] = 19;
				else
                                  mols->ism_types[i] = 24;
			}else if ( k == 3){
                               if( k_heavy == 1)
                                  mols->ism_types[i] = 25;
                                else if ( k_heavy == 2)
                                  mols->ism_types[i] = 24;
                                else if ( k_heavy == 3)
                                  mols->ism_types[i] = 23;
				else
                                  mols->ism_types[i] = 24;
			}else if ( k == 2){
                               if( k_heavy == 1)
                                  mols->ism_types[i] = 26;
                                else
                                  mols->ism_types[i] = 24;
			}else{
                               mols->ism_types[i] = 24;
			}

		}else if ( mols->atoms[i] == 10) {
                        mols->ism_types[i] = 27;
                }else if ( mols->atoms[i] == 8) {
                        mols->ism_types[i] = 28;
                }else if ( mols->atoms[i] == 7) {
                        mols->ism_types[i] = 30;
                }else if ( mols->atoms[i] ==102 ) {
                        mols->ism_types[i] = 31;
                }else if ( mols->atoms[i] == 100) {
                        mols->ism_types[i] = 29;
                }else{
			fprintf(stderr,"LCPO model has no parameters for atom no %i. Assuming Csp2. Atom type: %i. GAFF type: %i\n",i+1,mols->atoms[i],mols->gaff_types[i]);
                        mols->ism_types[i] = 24;
		}

	}

	*mymol = mols;
	return 0;
}


