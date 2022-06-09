#ifndef IC_DEFECT_READ_HEADER
#define IC_DEFECT_READ_HEADER

#include "parser.hpp"
#include "global_defect.hpp"
#include "metadata.hpp"
#include "background.hpp"


using namespace std;
using namespace LATfield2;

void generateIC_defects_read(DefectBase *defects, GlobalDefect  & defects_, string h5filename, string phifile, string pifile)
{
	defects = &defects_; 
	defects_.generate_init_cond(h5filename, phifile, pifile);
}

#endif
