#include "magnetism.h"
#include "global.h"

Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = new double[10];
}

Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

void Magnetism::compute_magnetization()
{
    if (NSPIN==2)
    {
        this->tot_magnetization = 0.00;
        this->abs_magnetization = 0.00;

		//chr.check_ne(chr.rho[0]);
		//chr.check_ne(chr.rho[1]);
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            double diff = chr.rho[0][ir] - chr.rho[1][ir];
            this->tot_magnetization += diff;
            this->abs_magnetization += abs(diff);
        }

        Parallel_Reduce::reduce_double_pool( this->tot_magnetization );
        Parallel_Reduce::reduce_double_pool( this->abs_magnetization );
        this->tot_magnetization *= ucell.omega / pw.ncxyz;
        this->abs_magnetization *= ucell.omega / pw.ncxyz;

		OUT(ofs_running,"total magnetism (Bohr mag/cell)",this->tot_magnetization);
		OUT(ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
		
		if(TWO_EFERMI)
		{
			OUT(ofs_running,"nelup",get_nelup());
			OUT(ofs_running,"neldw",get_neldw());
		}
		else
		{
			OUT(ofs_running,"nelec",ucell.nelec);
		}

//        cout << "\n tot_mag = " << setprecision(6) << this->tot_magnetization << " Bohr mag/cell" << endl;
  //      cout << " abs_mag = " << setprecision(6) << this->abs_magnetization << " Bohr mag/cell" << endl;
    }
	// noncolliear :
	else if(NSPIN==4)
	{
		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] = 0.00;
		this->abs_magnetization = 0.00;
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			double diff = sqrt(pow(chr.rho[1][ir], 2) + pow(chr.rho[2][ir], 2) +pow(chr.rho[3][ir], 2));
 
			for(int i=0;i<3;i++)this->tot_magnetization_nc[i] += chr.rho[i+1][ir];
			this->abs_magnetization += abs(diff);
		}
		Parallel_Reduce::reduce_double_pool( this->tot_magnetization_nc, 3 );
		Parallel_Reduce::reduce_double_pool( this->abs_magnetization );

		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] *= ucell.omega / pw.ncxyz;
		this->abs_magnetization *= ucell.omega / pw.ncxyz;
		ofs_running<<"total magnetism (Bohr mag/cell)"<<'\t'<<this->tot_magnetization_nc[0]<<'\t'<<this->tot_magnetization_nc[1]<<'\t'<<this->tot_magnetization_nc[2]<<'\n';
		OUT(ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
	}

    return;
}


double Magnetism::get_nelup(void)
{
	double nelup = 0.0;
	if(TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		nelup = 0.5 * ucell.nelec + 0.5 * tot_magnetization;
	}
	else
	{
		nelup = 0.5 * ucell.nelec;
	}
    return nelup;

    // for constrained magnetism calculation : not used now
    // nelup = ( nelec + mcons(3,1) ) * 0.5D0

//	double nelup = 0.0;
//	for(int i=0; i<pw.ntype; i++)
//	{
//		nelup += ucell.nelec * (1.0+start_magnetization[i])/2.0/pw.ntype;
//	}
//	return nelup;
}


double Magnetism::get_neldw(void)
{
	double neldw = 0.0;
	if(TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		neldw = 0.5 * ucell.nelec - 0.5 * tot_magnetization;
	}
	else
	{
		neldw = 0.5 * ucell.nelec;
	}
    return neldw ;

//	double neldw = 0.0;
//	for(int i=0; i<pw.ntype; i++)
//	{
//		neldw += pw.nelec * (1.0-start_magnetization[i])/2.0/pw.ntype;
//	}
//	return neldw;
}

// LiuXh add 2023.03.02
void Magnetism::cal_local_mag()
{
    int natom = ucell.nat;
    //double* local_mag;
    //double* local_mag_nc;

    //if(NSPIN==2)
    //{
    //    local_mag = new double[natom];
    //    ZEROS(local_mag, natom);
    //}
    //else if(NSPIN==4)
    //{
    //    local_mag_nc = new double[natom];
    //    ZEROS(local_mag_nc, natom);
    //}

    double* local_mag;
    local_mag = new double[natom];
    ZEROS(local_mag, natom);

    int ncx = pw.ncx;
    int ncy = pw.ncy;
    int nczp = pw.nczp;
    int nczp_start = pw.nczp_start;

    //double mcell_pos[3] = {0.0};
    Vector3<double> mcell_pos;
    Vector3<double> atom_pos;
    Vector3<double> dis_mcell_atom;
    double distance = 0.0;
    double atom_rcut = 5.0; // Unit: Bohr
    for(int T0=0; T0<ucell.ntype; T0++)
    {
        for(int I0=0; I0<ucell.atoms[T0].na; I0++)
        {
            int iat = ucell.itia2iat(T0,I0);
            atom_pos = ucell.atoms[T0].tau[I0];

            for(int i=0; i<ncx; i++)
            {
                for(int j=0; j<ncy; j++)
                {
                    for(int k=nczp_start; k<nczp_start+nczp; k++)
                    {
                        //mcell_pos[0] = i*GridT.meshcell_vec1[0] + j*GridT.meshcell_vec2[0] + k*GridT.meshcell_vec3[0];
                        //mcell_pos[1] = i*GridT.meshcell_vec1[1] + j*GridT.meshcell_vec2[1] + k*GridT.meshcell_vec3[1];
                        //mcell_pos[2] = i*GridT.meshcell_vec1[2] + j*GridT.meshcell_vec2[2] + k*GridT.meshcell_vec3[2];
                        int mcell_index = i*ncy*nczp + j*nczp + (k-nczp_start);

                        mcell_pos.x = i*GridT.meshcell_vec1[0] + j*GridT.meshcell_vec2[0] + k*GridT.meshcell_vec3[0];
                        mcell_pos.y = i*GridT.meshcell_vec1[1] + j*GridT.meshcell_vec2[1] + k*GridT.meshcell_vec3[1];
                        mcell_pos.z = i*GridT.meshcell_vec1[2] + j*GridT.meshcell_vec2[2] + k*GridT.meshcell_vec3[2];

                        dis_mcell_atom = atom_pos - mcell_pos;
                        distance = dis_mcell_atom.norm();
                        if(distance < atom_rcut)
                        {
                            if(NSPIN==2)
                            {
                                double diff = chr.rho[0][mcell_index] - chr.rho[1][mcell_index];
                                local_mag[iat] += diff;
                            }
                            if(NSPIN==4)
                            {
                                double diff = sqrt(pow(chr.rho[1][mcell_index], 2) + pow(chr.rho[2][mcell_index], 2) +pow(chr.rho[3][mcell_index], 2));
                                //local_mag_nc[iat] += diff;
                                local_mag[iat] += diff;
                            }
                        }
                    }
                }
            }
        }
    }

    //if(NSPIN==2)
    //{
    //    Parallel_Reduce::reduce_double_pool(local_mag, natom);
    //    for(int iat=0; iat<natom; iat++)
    //    {
    //        local_mag[iat] *= ucell.omega / pw.ncxyz;
    //        OUT(ofs_running, "atom label", iat);
    //        OUT(ofs_running, "local magnetic moment", local_mag[iat]);
    //    }

    //    delete[] local_mag;
    //}
    //else if(NSPIN==4)
    //{
    //    Parallel_Reduce::reduce_double_pool(local_mag_nc, natom);
    //    for(int iat=0; iat<natom; iat++)
    //    {
    //        local_mag_nc[iat] *= ucell.omega / pw.ncxyz;
    //        OUT(ofs_running, "atom label", iat);
    //        OUT(ofs_running, "local magnetic moment", local_mag_nc[iat]);
    //    }

    //    delete[] local_mag_nc;
    //}

    Parallel_Reduce::reduce_double_pool(local_mag, natom);
    for(int iat=0; iat<natom; iat++)
    {
        local_mag[iat] *= ucell.omega / pw.ncxyz;
        OUT(ofs_running, "atom label", iat);
        OUT(ofs_running, "local magnetic moment", local_mag[iat]);
    }

    delete[] local_mag;
}
