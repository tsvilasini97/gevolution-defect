#ifndef IC_DEFECT_PREVOLUTION_HEADER
#define IC_DEFECT_PREVOLUTION_HEADER

#include "parser.hpp"
#include "global_defect.hpp"
#include "metadata.hpp"
#include "background.hpp"


using namespace std;
using namespace LATfield2;

void generateIC_defects_prevolution(cosmology & cosmo, defects_metadata & defects_sim, const double fourpiG, DefectBase *defects, GlobalDefect  & defects_, double z_ic_defect,double z_in ,double dx, string h5filename)
{

	defects = &defects_; 

	defects_.generate_init_cond(h5filename);

	COUT << "Starting the prevolution of defects!" << endl;

	double endz = z_in;
	int i =0;
	double tmp = 1 / (z_ic_defect + 1);
	double sum = 0;
	double sum_;
	int N = 10;
	int zdefectscount = 0;
	
	COUT << "initial redshift is: " << tmp << endl;
	defects->compute_Tuv_defect(tmp);
	defects->write_Tuv_defect(h5filename,999);
	
	ofstream phifile;
	if(parallel.rank() == 0)
	{
		phifile.open (h5filename + "average_rho_phi_defect.txt",ios::trunc);
		phifile << "#i" << " " << "#z" << " " << "#a" << " " << "#adotovera" << " " << "#average phi" << " " << "#average rho" << endl;
		phifile.close();
	}


	do 
	{
		double vari = Hconf(tmp, fourpiG, cosmo);
		double DT = 0.1;
		double T;
		defects->update_phi(&DT);
		rungekutta4bg(tmp, fourpiG, cosmo, DT/2.0);
		defects->update_pi(&DT,&tmp,&vari);
		rungekutta4bg(tmp, fourpiG, cosmo, DT/2.0);

		defects->compute_Tuv_defect(tmp);
		

		double averagephi = defects_.averagephi();
		double averagerho = defects_.averagerhodefect();

		T += DT;
		i++;
		COUT << "Current redshift is at: " << tmp << endl;
		if(parallel.rank() == 0)
		{
			phifile.open (h5filename + "average_rho_phi_defect.txt",std::ios_base::app);
			phifile << i << " " << (1/tmp) - 1 << " " << tmp << " " << vari << " " << averagephi << " " << averagerho << endl;
			phifile.close();
		}

//		sum += averagephi;
//		if(i%N == 0)
//		{
//			sum = sum/N;
//			sum_ = i;
//			COUT << sum << " " ;
//			if(tmp <= 1/(endz + 1))
//			{
//				sum = 0;
//			}
//		}
		
//		if(i%500 == 0)
//		{
//			defects->write_Tuv_defect(h5filename,i);
//			defects_.test_global_defect(h5filename, i, (1/tmp) - 1, DT, T);
//		}
		
		if (zdefectscount < defects_sim.num_defect_output && 1. / tmp < defects_sim.z_defects[zdefectscount] + 1.)
			{	
				COUT << " Number of outputs are = " << defects_sim.num_defect_output  << endl;
				
				COUT << COLORTEXT_BLUE << " writing defect output " << COLORTEXT_RESET << " at z = " << ((1./tmp) - 1.) << endl;
				defects->defects_stat_output();

				defects->write_Tuv_defect(h5filename,zdefectscount);
				defects->writedefectSnapshots(h5filename, zdefectscount);

				zdefectscount++;
			}

	}
	while (tmp <= 1/(endz + 1) );

	sum = sum /(i-sum_);
	double z = (1-tmp)/tmp;
	COUT << endl << endl << "The total steps for defects prevolution is: " << i << " and the redshift at the end of defect prevolution now is: " << z << endl;

#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif

	if(sum < 0.988 || sum > 1.02)
	{
	    defects_.write_Tuv_defect(h5filename, i);
		COUT << endl << COLORTEXT_RED << " error " << COLORTEXT_RESET << ": The defect has not reached stable condition !!!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	else if(sum >= 0.988 || sum <= 1.02)
	{
		defects_.write_Tuv_defect(h5filename, i);
		COUT << endl << COLORTEXT_BLUE << "The defect has reached stable condition in prevolution " << COLORTEXT_RESET << endl;
	}


}

#endif
