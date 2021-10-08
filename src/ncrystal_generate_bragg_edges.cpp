//////////////////////////////////////////////////////////////////////////////////
//  This is a C++ wrapper to get the coherent elastic cross section from        //
//  NCrystal and pass it back to LEAPR module of NJOY as pairs of energy        //
//  and the cross section at each bragg edge. We also pass back the value for   //
//  for bound coherent and incoherent cross section, as well as the free atom   //
//  cross section.                                                              //
//  The temperature is being passed in from the current temperature loop        //
//  in LEAPR. We also do some checks to make sure the user is using correct     //
//  physics and proper element is being treated. User has to make sure, in case //
//  of own .NCMAT files to user v4 format.                                      //
//////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/NCrystal.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCMatCfg.hh"
#include <iostream>
#include <cstdlib>
#include <cstdio>

namespace NC = NCrystal;

extern "C" 
{
  void generate_bragg_edges(char* s, int* nbragg, double* data, double* current_temp, int* maxb, unsigned int* atomic_z, unsigned int* atomic_a, double* c_incoherent_bound_xs, double* c_coherent_bound_xs, double* c_spr, double* c_ncrystal_msd, double* c_incoherent_fraction, double* c_fraction, int* c_redistribute_this);
}

void generate_bragg_edges(char* s, int* nbragg, double* data, double* current_temp, int* maxb, unsigned int* atomic_z, unsigned int* atomic_a, double* c_incoherent_bound_xs, double* c_coherent_bound_xs, double* c_spr, double* c_ncrystal_msd, double* c_incoherent_fraction, double* c_fraction, int* c_redistribute_this)
{
  NC::libClashDetect();
  
  //create an Info object from the passed in string with current_temp
  NC::MatCfg cfg( s );
  cfg.set_temp(NC::Temperature{*current_temp});
  auto info = NCrystal::createInfo(cfg);

  //Assert that the material is not a single crystal, this only for works for polycrystalline materials
  if ( cfg.isSingleCrystal())
    NCRYSTAL_THROW(MissingInfo,"This NJOY plugin is not to be used with single crystals!!");

  //get the bound incoherent and coherent cross sections, as well as the free atom cross section
  unsigned nfound = 0;
  *c_incoherent_fraction = -1;
  unsigned int minimum_incoherent_contribution_Z = -1;
  unsigned int minimum_incoherent_contribution_A = -1;
  *c_redistribute_this = 0;
  for( auto& di : info->getDynamicInfoList() ) {
    auto di_vdos_ptr = dynamic_cast<const NC::DI_VDOS*>(di.get());
    if (di_vdos_ptr){
      const NC::DI_VDOS& di_vdos = *di_vdos_ptr;
      const NC::AtomData& atomData = di_vdos.atomData();
      
      if ((di_vdos.fraction()*atomData.incoherentXS().dbl() < *c_incoherent_fraction) ||
        (*c_incoherent_fraction == -1)){
        if ( di_vdos.fraction() == 0.0) {
          *c_incoherent_fraction = 0.0;
        } else {
          *c_incoherent_fraction = di_vdos.fraction()/(1.0 - di_vdos.fraction())*atomData.incoherentXS().dbl();
        }
        minimum_incoherent_contribution_Z = atomData.Z();
        if ( atomData.isNaturalElement() ){
          minimum_incoherent_contribution_A = 0;
        }
        else {
          minimum_incoherent_contribution_A = atomData.A();
        }
      }
      
      if (!atomData.isElement())
        continue;//This code can only handle mixtures with well defined Z
      if (!(atomData.Z() == *atomic_z) )
        continue;//Wrong Z
      if ( ( ( atomData.isSingleIsotope() && atomData.A() == *atomic_a ) ) || ( atomData.isNaturalElement() && *atomic_a == 0 ) ) {
        //Find correct single isotopes or natural element.
        *c_incoherent_bound_xs=atomData.incoherentXS().dbl();
        *c_coherent_bound_xs=atomData.coherentXS().dbl();
        *c_spr=atomData.freeScatteringXS().dbl();
        *c_ncrystal_msd=(di->correspondingAtomInfo())->msd().value();
        *c_fraction=di_vdos.fraction();
        ++nfound;
      }
    }
  }  
  
   //Assert that the material requested in LEAPR is also in .NCMAT file
  if ( nfound == 0 )
      NCRYSTAL_THROW(MissingInfo,"The requested element cannot be found in the .NCMAT file!!");
  if ( nfound > 1 )
      NCRYSTAL_THROW(MissingInfo,"The requested element has multiple roles in the .NCMAT file!!");

  //Mark element for redistribution if iel=99 is used
  if ((minimum_incoherent_contribution_Z == *atomic_z) && 
      (minimum_incoherent_contribution_A == *atomic_a)) *c_redistribute_this = 1;
  
  //Loop over all the planes and send energy of the bragg edge and energy multipled by cross section for that edge to LEAPR
  int n = static_cast<int>(info->nHKL());
  if (!( *maxb > 2*n )){
    NCRYSTAL_THROW(CalcError,"Number of bragg edges is bigger than maxb, increase maxb in leapr.f90");
  }
  *nbragg = n;
  
  const double xsectfact = 0.5/(info->getStructureInfo().volume * info->getStructureInfo().n_atoms);

  for ( auto& hkl : info->hklList()) {
    const double wl = 2.0*hkl.dspacing;
    const double E = NC::wl2ekin(wl);  
    const double fdm = hkl.fsquared * hkl.multiplicity * hkl.dspacing;
    *data++ = E;
    *data++ = E * fdm * xsectfact * wl * wl;
  }
}
