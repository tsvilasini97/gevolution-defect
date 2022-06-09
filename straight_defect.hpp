#ifndef STRAIGHT_DEFECT_HPP
#define STRAIGHT_DEFECT_HPP

#include "defect_base.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>

#include "LATfield2.hpp"
using namespace LATfield2;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"
#include "TUV_comm.hpp"

//using namespace std;


class StraightDefect:public DefectBase
{
private:
	double v;
	double x_;
	double x0;
	double y0;

	vector<vector<double>> Tuv_straight_;

public:
	void initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim);
	void update_phi(double *dt);
	void compute_Tuv_defect(double a);
	void write_Tuv_defect(string h5filename, const int snapcount);

	void projectionCIC_straight_defect(double a);
	void projection_Ti0_straight_defect(double a);
	void write_straight_defect_position(string h5filename, const int snapcount, double z, double dtau, double tau);
};

void StraightDefect::initialize(Lattice * lat, Lattice * klat, double *dx, metadata * sim, defects_metadata * defects_sim)
{
	dx_ = dx;
	lat_ = lat;
	klat_ = klat;
	defects_sim_ = defects_sim;
	sim_ = sim;

	x0 = defects_sim_->string_x_pos;
	y0 = defects_sim_->string_y_pos;
	v = defects_sim_->string_velocity;
	
	
	x0 /= sim_->boxsize;
	y0 /= sim_->boxsize;
	
	x_ = x0;

	double gamma;
	double mu;

	mu = 75.3154 / sim_->boxsize / sim_->boxsize;

	COUT << "the val of mu is: " << mu << endl;
	gamma = 1/sqrt(1 - v*v);

	Tuv_defect_.initialize(*lat_,4,4,symmetric);
	Tuv_defect_.alloc();

	double t00 = mu * gamma / *dx_ / *dx_;
	double t11 = mu * gamma * v*v/ *dx_ / *dx_ ;
	double t33 = - mu / gamma / *dx_ / *dx_;
	double t10 = mu * gamma * v / *dx_ / *dx_;

	int m = 4;
	int n = 4;
	Tuv_straight_.resize(m, vector<double>(n,0));
	Tuv_straight_[0][0] = t00;
	Tuv_straight_[0][1] = t10;
	Tuv_straight_[1][1] = t11;
	Tuv_straight_[3][3] = t33;
	Tuv_straight_[1][0] = t10;
	
//	projectionCIC_straight_defect();

	COUT << "initialization of straight defect done!" << endl;
}

void StraightDefect::update_phi(double *dt)
{ 
	x_ = x_ + v * *dt;
	if(x_>1)
	{
#ifdef LATFIELD2_HPP
		COUT << COLORTEXT_RED << "The string has reached the edge of the box, so aborting further computation now!!" << COLORTEXT_RESET << endl;
		parallel.abortForce();
#endif
	}
}



void StraightDefect::compute_Tuv_defect(double a)
{
	/// Emptying Halo
	
	long sizeLocalGross[3];
	long sizeLocal[3];
	int comp = 10;
	int halo = Tuv_defect_.lattice().halo();
	
	for(int i=0;i<3;i++)
	{
		sizeLocal[i]=Tuv_defect_.lattice().sizeLocal(i);
		sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
	}

	for(int k=0;k<sizeLocalGross[2];k++)
	{
		for(int j=0;j<sizeLocalGross[1];j++)
		{
			for(int i=0;i<sizeLocalGross[0];i++)
			{
				for(int c=0;c<comp;c++)Tuv_defect_(setIndex(sizeLocalGross,i,j,k),c) = 0;
			}
		}
	}
	projectionCIC_straight_defect(a);
	TUV_comm(&Tuv_defect_);
	
}



void StraightDefect::write_Tuv_defect(string h5filename, const int snapcount)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%05d", snapcount);

#ifdef EXTERNAL_IO
	COUT << "Currently defect snapshot does not work with external IO" << endl;
#else
	Tuv_defect_.saveHDF5(h5filename + filename_def + "_Tuv_straight_defect_.h5" );
#endif
}



void StraightDefect::write_straight_defect_position(string h5filename, const int snapcount, double z, double dtau, double tau)
{
	char filename_def[2*PARAM_MAX_LENGTH+24];
	sprintf(filename_def, "%05d", snapcount);
	if(snapcount ==0)
	{
		ofstream straightdefectfile;
		if(parallel.rank() == 0)
		{
			straightdefectfile.open (h5filename + "straight_defect.txt",ios::trunc);
			straightdefectfile << "#z" << " " << "#dt" << " " << "#tau" << " " << "#x" << " " << "#snapcount"<< endl;
			straightdefectfile.close();
			
			straightdefectfile.open (h5filename + "straight_defect.txt",std::ios_base::app);
			straightdefectfile << z << " " << dtau << " " << tau << " " << x_ << " " << snapcount << endl;
			straightdefectfile.close();
		}
	}
	else
	{
		ofstream straightdefectfile;
		if(parallel.rank() == 0)
		{
			straightdefectfile.open (h5filename + "straight_defect.txt",std::ios_base::app);
			straightdefectfile << z << " " << dtau << " " << tau << " " << x_ << " " << snapcount << endl;
			straightdefectfile.close();
		}
	}
}

void StraightDefect::projectionCIC_straight_defect(double a)
{	
	Site xTuv(Tuv_defect_.lattice());

	double cord[3];
	double rescalPos[3];
	double rescalPosDown[3];
	double string_coord[2];

//	std::cout << std::setprecision(12) << " String x = " << x_ << " and y position = " << y0 << endl;
	COUT << " String x = " << x_ << " and y position = " << y0 << endl;
	cord[0] = (x_ - fmod(x_,*dx_)) / *dx_;
	cord[1] = (y0 - fmod(y0,*dx_)) / *dx_;

	string_coord[0] = x_;
	string_coord[1] = y0;

	for(int zcord = 0; zcord < sim_->numpts; zcord++)
	{
//		COUT << "coordinates === " << cord[0] << " " << cord[1] << " " << zcord << endl;
		if(xTuv.setCoord(cord[0], cord[1], zcord))
		{
			for(int i =0;i<2;i++)
			{
				rescalPos[i]= string_coord[i] - cord[i]* *dx_;
				rescalPosDown[i]= *dx_ - rescalPos[i];
			}
			for (int i=0; i<3; i++)
			{
				for (int j=0; j<=i; j++)
				{
					Tuv_defect_(xTuv,i,j)        = rescalPosDown[0] * rescalPosDown[1] * Tuv_straight_[i][j] / *dx_ / a / a / *dx_;
					Tuv_defect_(xTuv+1,i,j)      = rescalPosDown[0] * rescalPos[1] * Tuv_straight_[i][j] / *dx_ / a / a / *dx_;
					Tuv_defect_(xTuv+0,i,j)      = rescalPos[0] * rescalPosDown[1] * Tuv_straight_[i][j] / *dx_  / a / a / *dx_;
					Tuv_defect_(xTuv+0+1,i,j)    = rescalPos[0] * rescalPos[1] * Tuv_straight_[i][j] / *dx_ / a / a / *dx_; 
					
//					Tuv_defect_(xTuv+0+1+2,i,j)  = rescalPos[0] * rescalPos[1] * Tuv_straight_[i][j];
//					Tuv_defect_(xTuv+2,i,j)      = rescalPosDown[0] * rescalPosDown[1] * Tuv_straight_[i][j];
//					Tuv_defect_(xTuv+1+2,i,j)    = rescalPosDown[0] * rescalPos[1] * Tuv_straight_[i][j];
//					Tuv_defect_(xTuv+0+2,i,j)    = rescalPos[0] * rescalPosDown[1] * Tuv_straight_[i][j];
				}
			}
		}
	}

}

void StraightDefect::projection_Ti0_straight_defect(double a)
{	
	Site xTuv(Tuv_defect_.lattice());
	
    double rescalPos[3];
    double rescalPosDown[3];
    
    double string_pos[2];
    double string_ref_coord[2];
    
    string_pos[0] = x_;
    string_pos[1] = y0; 
    
    for (int i=0; i<2; i++) string_ref_coord[i] = (string_pos[i]-fmod(string_pos[i],*dx_))/ *dx_;
	
	for(int i=0;i<sim_->numpts;i++)
	{
 		
		if(xTuv.setCoord(string_ref_coord[0], string_ref_coord[1], i))
		{
 			for(int k =0;k<2;k++)
		    {
		        rescalPos[k]= (string_pos[k] - string_ref_coord[k] * *dx_) / *dx_;
		        rescalPosDown[k]= 1.0l - rescalPos[k];
		    }

			
			for (int m=1; m<4; m++)
			{
				Tuv_defect_(xTuv,m,0)        = rescalPos[1] * Tuv_straight_[m][0]  / a / a;
				Tuv_defect_(xTuv+1,m,0)      = rescalPosDown[1] * Tuv_straight_[m][0]  / a / a ;
			}
		}
	}
}

#endif
