#ifndef MAGNETISM_H
#define MAGNETISM_H

#include "tools.h"

class Magnetism
{
public:
    // constructor and deconstructor
    Magnetism();
    ~Magnetism();

    // notice : becast is done in unitcell
    double *start_magnetization;

    // tot_magnetization : majority spin - minority spin (nelup - neldw).
    double tot_magnetization;
    double tot_magnetization_nc[3];
    double abs_magnetization;

	// nelup and neldw is different only the magnetism is fixed,
	// otherwise they are half of the total electrons.
    double get_nelup();
    double get_neldw();

    double nelup_percent(const int &it) {
        return 0.5 * (1.0 + start_magnetization[it]) ;
    }
    double neldw_percent(const int &it) {
        return 0.5 * (1.0 - start_magnetization[it]) ;
    }

    void compute_magnetization();

    //LiuXh add 2023.03.02
    void cal_local_mag();

private:

};

/*
 A comment about variables nelup, neldw, multiplicity and tot_magnetization:
 All these variables contain the same information and must be kept harmonized.
 Variables nelup and neldw will be removed in future versions of the code.
 Variables multiplicity and tot_magnetization, though redundent will probably
 coexist since multiplicity is the more natural way (?)for defining the spin
 configuratio in the quantum-chemistry community while tot_magnetization is
 more natural (?) when dealing with extended systems.
*/

#endif
