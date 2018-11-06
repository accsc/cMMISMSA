/**
 *
 *      @file amber.h
 *      @brief AMBER files reader headers
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 28/01/2014
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */



/**
 * Structure to AMBER vdw params
 */
typedef struct{
/** Total number of atoms */
int n_atoms;
/** Number of different atom types */
int ntypes;
/** Atom type index array */
int *atom_type_index;  
/** Atoms that are not considered for anything like WATs */
int *excluded_list;
/** Non bonded parm index */
int *nb_index;
/** Lennard Jones A coeff */
float *lj_A;     
/** Lennard Jones B coeff */
float *lj_B;  
/** mdcrd frames index */
long *mdcrd_index;
/** Box */
int box;
/** frames */
int frames;
/** XTC file handle if needed */
XDRFILE *xdr;
} TOP;


