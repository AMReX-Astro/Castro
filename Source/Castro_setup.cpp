#include <winstd.H>
#include <cstdio>

#include "LevelBld.H"
#include <ParmParse.H>
#include "Castro.H"
#include "Castro_F.H"
#include <Derive_F.H>
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif

using std::string;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM >= 2)
static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}
#endif

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

void
Castro::variableSetUp ()
{

  // Castro::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job


  // initialize the start time for our CPU-time tracker
  startCPUTime = ParallelDescriptor::second();


    BL_ASSERT(desc_lst.size() == 0);

    // Initialize the runtime parameters for any of the external
    // microphysics
    extern_init();


    // Initialize the network
    network_init();

    // Get options, set phys_bc
    read_params();
    //
    // Set number of state variables and pointers to components
    //
    int use_sgs = 0;

    int cnt = 0;
    Density = cnt++;
    Xmom = cnt++;
#if (BL_SPACEDIM >= 2)
    Ymom = cnt++;
#endif
#if (BL_SPACEDIM == 3)
    Zmom = cnt++;
#endif
    Eden = cnt++;
    Eint = cnt++;
#ifdef SGS
    Esgs = cnt++;
    use_sgs = 1;
#endif
    Temp = cnt++;

    NumAdv = 0;

    // if we are in 2-d with rotation, carry an advected quantity that
    // will be the velocity through the plane of the domain (i.e. w).
    // This can participate in the Coriolis force
#ifdef ROTATION
#if (BL_SPACEDIM == 2)
    if ( Geometry::IsRZ()) {
      NumAdv = 1;
    }
#endif
#endif
    if (NumAdv > 0)
    {
        FirstAdv = cnt;
        cnt += NumAdv;
    }

    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);

    if (NumSpec > 0)
    {
        FirstSpec = cnt;
        cnt += NumSpec;
    }

    // Get the number of auxiliary quantities from the network model.
    BL_FORT_PROC_CALL(GET_NUM_AUX, get_num_aux)(&NumAux);

    if (NumAux > 0)
    {
        FirstAux = cnt;
        cnt += NumAux;
    }

    NUM_STATE = cnt;

    // Define NUM_GROW from the f90 module.
    BL_FORT_PROC_CALL(GET_METHOD_PARAMS, get_method_params)(&NUM_GROW);

    const Real run_strt = ParallelDescriptor::second() ; 

#ifndef ROTATION
    static Real rotational_period = 0;
#endif

    // we want const_grav in F90, get it here from parmparse, since it
    // it not in the Castro namespace
    ParmParse pp("gravity");
    Real const_grav = 0;
    pp.query("const_grav", const_grav);


    BL_FORT_PROC_CALL(SET_METHOD_PARAMS, set_method_params)
        (dm, Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux, 
         NumAdv, difmag, small_dens, small_temp, small_pres, 
         allow_negative_energy,ppm_type,ppm_reference,
	 ppm_trace_grav,ppm_temp_fix,ppm_tau_in_tracing,ppm_predict_gammae,
	 ppm_reference_edge_limit,
	 ppm_flatten_before_integrals,
	 ppm_reference_eigenvectors,
	 hybrid_riemann, use_colglaz, use_flattening, 
         transverse_use_eos, transverse_reset_density, transverse_reset_rhoe,
         cg_maxiter, cg_tol,
         use_pslope, 
	 grav_source_type, do_sponge,
         normalize_species,fix_mass_flux,use_sgs,rotational_period,
	 const_grav, deterministic);

    Real run_stop = ParallelDescriptor::second() - run_strt;
 
    ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
 
    if (ParallelDescriptor::IOProcessor())
        std::cout << "\nTime in set_method_params: " << run_stop << '\n' ;

    int coord_type = Geometry::Coord();

    Real xmin = Geometry::ProbLo(0);
    Real xmax = Geometry::ProbHi(0);
#if (BL_SPACEDIM >= 2)
    Real ymin = Geometry::ProbLo(1);
    Real ymax = Geometry::ProbHi(1);
#else
    Real ymin = 0.0;
    Real ymax = 0.0;
#endif
#if (BL_SPACEDIM == 3)
    Real zmin = Geometry::ProbLo(2);
    Real zmax = Geometry::ProbHi(2);
#else
    Real zmin = 0.0;
    Real zmax = 0.0;
#endif
    
    BL_FORT_PROC_CALL(SET_PROBLEM_PARAMS, set_problem_params)
         (dm,phys_bc.lo(),phys_bc.hi(),
	  Outflow,Symmetry,SlipWall,NoSlipWall,coord_type,
	  xmin,xmax,ymin,ymax,zmin,zmax);

    Interpolater* interp = &cell_cons_interp;

#ifdef RADIATION
    // cell_cons_interp is not conservative in spherical coordinates.
    // We could do this for other cases too, but I'll confine it to
    // neutrino problems for now so as not to change the results of
    // other people's tests.  Better to fix cell_cons_interp!

    if (Geometry::IsSPHERICAL() && Radiation::nNeutrinoSpecies > 0) {
      interp = &pc_interp;
    }
#endif

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

#ifdef RADIATION
    int ngrow_state = 1;
#else
    int ngrow_state = 0;
#endif

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,ngrow_state,NUM_STATE,
                           interp,state_data_extrap,store_in_checkpoint);

#ifdef GRAVITY
    store_in_checkpoint = true;
    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,BL_SPACEDIM,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

#ifdef LEVELSET
    store_in_checkpoint = true;
    desc_lst.addDescriptor(LS_State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,1,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

#ifdef REACTIONS
    // Components 0:Numspec-1         are      omegadot_i
    // Component    NumSpec            is      enuc =      (eout-ein)
    // Component    NumSpec+1          is  rho_enuc= rho * (eout-ein)
    store_in_checkpoint = true;
    store_in_checkpoint = false;
    desc_lst.addDescriptor(Reactions_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NumSpec+2,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

#ifdef SGS
    // Component 0: prod_sgs
    // Component 1: diss_sgs
    // Component 2: turbulent forcing
    store_in_checkpoint = true;
    desc_lst.addDescriptor(SGS_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,3,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

    Array<BCRec>       bcs(NUM_STATE);
    Array<std::string> name(NUM_STATE);
    
    BCRec bc;
    cnt = 0;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
#if (BL_SPACEDIM >= 2)
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
#endif
#if (BL_SPACEDIM == 3)
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
#endif
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
#ifdef SGS
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_K";
#endif
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";

    for (int i=0; i<NumAdv; ++i)
    {
        char buf[64];
        sprintf(buf, "adv_%d", i);
        cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = string(buf);
    }

    // Get the species names from the network model.
    char* spec_names[NumSpec];
    for (int i = 0; i < NumSpec; i++) {
          int len = 20;
          Array<int> int_spec_names(len);
          // This call return the actual length of each string in "len" 
          BL_FORT_PROC_CALL(GET_SPEC_NAMES, get_spec_names)(int_spec_names.dataPtr(),&i,&len);
          spec_names[i] = new char[len+1];
          for (int j = 0; j < len; j++) 
             spec_names[i][j] = int_spec_names[j];
          spec_names[i][len] = '\0';
    }

    if ( ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Species: " << std::endl;
        for (int i = 0; i < NumSpec; i++)  
           std::cout << spec_names[i] << ' ' << ' ';
        std::cout << std::endl;
    } 

    for (int i=0; i<NumSpec; ++i)
    {
        cnt++; 
        set_scalar_bc(bc,phys_bc); 
        bcs[cnt] = bc; 
        string spec_string(spec_names[i]);
        name[cnt] = "rho_" + spec_string;
    }

    // Get the auxiliary names from the network model.
    char* aux_names[NumAux];
    for (int i = 0; i < NumAux; i++) {
          int len = 20;
          Array<int> int_aux_names(len);
          // This call return the actual length of each string in "len"
          BL_FORT_PROC_CALL(GET_AUX_NAMES, get_aux_names)(int_aux_names.dataPtr(),&i,&len);
          aux_names[i] = new char[len+1];
          for (int j = 0; j < len; j++)
             aux_names[i][j] = int_aux_names[j];
          aux_names[i][len] = '\0';
    }

    if ( ParallelDescriptor::IOProcessor())
    {
        std::cout << NumAux << " Auxiliary Variables: " << std::endl;
        for (int i = 0; i < NumAux; i++)
           std::cout << aux_names[i] << ' ' << ' ';
        std::cout << std::endl;
    }

    for (int i=0; i<NumAux; ++i)
    {
        cnt++;
        set_scalar_bc(bc,phys_bc);
        bcs[cnt] = bc;
        string aux_string(aux_names[i]);
        name[cnt] = "rho_" + aux_string;
    }

    desc_lst.setComponent(State_Type,
                          Density,
                          name,
                          bcs,
                          BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL,ca_denfill),
                                    BL_FORT_PROC_CALL(CA_HYPFILL,ca_hypfill)));

#ifdef GRAVITY
    if (do_grav) {
       set_x_vel_bc(bc,phys_bc);
       desc_lst.setComponent(Gravity_Type,0,"grav_x",bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVXFILL,ca_gravxfill)));
    }
#if (BL_SPACEDIM > 1)
    if (do_grav) {
       set_y_vel_bc(bc,phys_bc);
       desc_lst.setComponent(Gravity_Type,1,"grav_y",bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVYFILL,ca_gravyfill)));
    }
#endif
#if (BL_SPACEDIM > 2)
    if (do_grav) {
       set_z_vel_bc(bc,phys_bc);
       desc_lst.setComponent(Gravity_Type,2,"grav_z",bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVZFILL,ca_gravzfill)));
    }
#endif
#endif

#ifdef LEVELSET
    desc_lst.setComponent(LS_State_Type,0,"LSphi",bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_PHIFILL,ca_phifill)));
#endif

#ifdef REACTIONS
    std::string name_react;
    for (int i=0; i<NumSpec; ++i)
    {
       set_scalar_bc(bc,phys_bc);
       string aux_string(spec_names[i]);
       name_react = "omegadot_" + aux_string;
       desc_lst.setComponent(Reactions_Type, i, name_react, bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
    }
    desc_lst.setComponent(Reactions_Type, NumSpec  , "enuc", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
    desc_lst.setComponent(Reactions_Type, NumSpec+1, "rho_enuc", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
#endif

#ifdef SGS
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(SGS_Type, 0, "prod_sgs", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL,ca_sgsfill)));
    desc_lst.setComponent(SGS_Type, 1, "diss_sgs", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL,ca_sgsfill)));
    desc_lst.setComponent(SGS_Type, 2, "turb_src", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL,ca_sgsfill)));
#endif

#ifdef RADIATION
    int ngrow = 1;
    int ncomp = Radiation::nGroups;
    desc_lst.addDescriptor(Rad_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, ngrow, ncomp,
                           interp);
    set_scalar_bc(bc,phys_bc);

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Radiation::nGroups = " << Radiation::nGroups << std::endl;
      std::cout << "Radiation::nNeutrinoSpecies = "
           << Radiation::nNeutrinoSpecies << std::endl;
      if (Radiation::nNeutrinoSpecies > 0) {
        std::cout << "Radiation::nNeutrinoGroups  = ";
        for (int n = 0; n < Radiation::nNeutrinoSpecies; n++) {
          std::cout << " " << Radiation::nNeutrinoGroups[n];
        }
        std::cout << std::endl;
        if (Temp != BL_SPACEDIM + 3) {
          BoxLib::Error("Neutrino solver assumes Temp == BL_SPACEDIM + 3");
        }
        if (Radiation::nNeutrinoGroups[0] > 0 &&
            NumAdv != 0) {
          BoxLib::Error("Neutrino solver assumes NumAdv == 0");
        }
        if (Radiation::nNeutrinoGroups[0] > 0 &&
            (NumSpec != 1 || NumAux != 1)) {
          BoxLib::Error("Neutrino solver assumes NumSpec == NumAux == 1");
        }
      }
    }

    char rad_name[10];
    if (!Radiation::do_multigroup) {
      desc_lst
        .setComponent(Rad_Type, Rad, "rad", bc,
                      BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
    }
    else {
      if (Radiation::nNeutrinoSpecies == 0 ||
          Radiation::nNeutrinoGroups[0] == 0) {
	for (int i = 0; i < Radiation::nGroups; i++) {
	  sprintf(rad_name, "rad%d", i);
	  desc_lst
            .setComponent(Rad_Type, i, rad_name, bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
	}
      }
      else {
	int indx = 0;
	for (int j = 0; j < Radiation::nNeutrinoSpecies; j++) {
	  for (int i = 0; i < Radiation::nNeutrinoGroups[j]; i++) {
	    sprintf(rad_name, "rads%dg%d", j, i);
	    desc_lst
              .setComponent(Rad_Type, indx, rad_name, bc,
                            BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
	    indx++;
	  }
	}
      }
    }

    // In the following ncomp == 0 would be fine, except that it
    // violates an assertion in StateDescriptor.cpp.  I think we could
    // run with ncomp == 0 just by taking out that assertion.
    ncomp = 2;
    ncomp = (Radiation::Test_Type_lambda) ? 1 : ncomp;
    int nspec = (Radiation::nNeutrinoSpecies>0) ? Radiation::nNeutrinoSpecies : 1;
    if (Radiation::Test_Type_Flux) {
      ncomp = nspec * BL_SPACEDIM;
      if (Radiation::Test_Type_lambda) {
	ncomp++;  // SGFLD only
      }
    }

    ngrow = 0;
    desc_lst.addDescriptor(Test_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, ngrow, ncomp,
                           &cell_cons_interp);
    set_scalar_bc(bc,phys_bc);

    if (Radiation::Test_Type_Flux) {
      Array<std::string> radname(3);
      if (Radiation::nNeutrinoSpecies>0) {
	radname[0] = "nue";
	radname[1] = "nuae";
	radname[2] = "numu";
      }
      else {
	radname[0] = "rad";
      }

      Array<std::string> dimname(BL_SPACEDIM);
      dimname[0] = "x";
#if (BL_SPACEDIM >= 2)
      dimname[1] = "y";
#endif      
#if (BL_SPACEDIM >= 3)
      dimname[2] = "z";
#endif      
      
      int icomp = 0;
      for (int idim=0; idim<BL_SPACEDIM; idim++) {
	for (int ispec=0; ispec<nspec; ispec++) {
	  desc_lst.setComponent(Test_Type, icomp, "F"+dimname[idim]+radname[ispec], bc,
				BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
	  icomp++;
	}
      }
      if (Radiation::Test_Type_lambda) {
	desc_lst.setComponent(Test_Type, icomp, "lambda", bc,
			      BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
      }
    }
    else if (Radiation::Test_Type_lambda) {
      desc_lst.setComponent(Test_Type, 0, "lambda", bc,
			    BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));      
    }
    else {
      char test_name[10];
      for (int i = 0; i < ncomp; i++){
	sprintf(test_name, "test%d", i);
	desc_lst.setComponent(Test_Type, i, test_name, bc,
			      BndryFunc(BL_FORT_PROC_CALL(CA_RADFILL,ca_radfill)));
      }
    }
#endif

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERPRES,ca_derpres),the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERKINENG,ca_derkineng),the_same_box);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Xmom,BL_SPACEDIM);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERSOUNDSPEED,ca_dersoundspeed),the_same_box);
    derive_lst.addComponent("soundspeed",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERMACHNUMBER,ca_dermachnumber),the_same_box);
    derive_lst.addComponent("MachNumber",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM == 1)
    //
    // Wave speed u+c
    //
    derive_lst.add("uplusc",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERUPLUSC,ca_deruplusc),the_same_box);
    derive_lst.addComponent("uplusc",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Wave speed u-c
    //
    derive_lst.add("uminusc",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERUMINUSC,ca_deruminusc),the_same_box);
    derive_lst.addComponent("uminusc",desc_lst,State_Type,Density,NUM_STATE);
#endif

    //
    // Gravitational forcing 
    //
#ifdef GRAVITY
//    derive_lst.add("rhog",IndexType::TheCellType(),1,
//                   BL_FORT_PROC_CALL(CA_RHOG,ca_rhog),the_same_box);
//    derive_lst.addComponent("rhog",desc_lst,State_Type,Density,1);
//    derive_lst.addComponent("rhog",desc_lst,Gravity_Type,0,BL_SPACEDIM);
#endif

    //
    // Entropy (S)
    //
    derive_lst.add("entropy",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERENTROPY,ca_derentropy),the_same_box);
    derive_lst.addComponent("entropy",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM > 1)
    //
    // Vorticity
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERMAGVORT,ca_dermagvort),grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom != Density+1)
       BoxLib::Error("We are assuming Xmom = Density + 1 in Castro_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,BL_SPACEDIM+1);
#endif

    //
    // Div(u)
    //
    derive_lst.add("divu",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERDIVU,ca_derdivu),grow_box_by_one);
    derive_lst.addComponent("divu",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("divu",desc_lst,State_Type,Xmom,BL_SPACEDIM);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DEREINT1,ca_dereint1),the_same_box);
    derive_lst.addComponent("eint_E",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DEREINT2,ca_dereint2),the_same_box);
    derive_lst.addComponent("eint_e",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERLOGDEN,ca_derlogden),the_same_box);
    derive_lst.addComponent("logden",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("StateErr",IndexType::TheCellType(),3,
                   BL_FORT_PROC_CALL(CA_DERSTATE,ca_derstate),grow_box_by_one);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Temp,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,FirstSpec,1);

    //
    // X from rhoX
    //
    for (int i = 0; i < NumSpec; i++){
      string spec_string(spec_names[i]);
      spec_string = "X("+spec_string+")";

      derive_lst.add(spec_string,IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,FirstSpec+i,1);
    }
    //
    // Velocities
    //
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Xmom,1);

#if (BL_SPACEDIM >= 2)
    derive_lst.add("y_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Ymom,1);
#endif

#if (BL_SPACEDIM == 3)
    derive_lst.add("z_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Zmom,1);
#endif

#ifdef SGS
    derive_lst.add("K",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("K",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("K",desc_lst,State_Type,Esgs,1);

    derive_lst.add("forcex",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERFORCEX,ca_derforcex),the_same_box);
    derive_lst.addComponent("forcex",desc_lst,State_Type,Density,1);

    derive_lst.add("forcey",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERFORCEY,ca_derforcey),the_same_box);
    derive_lst.addComponent("forcey",desc_lst,State_Type,Density,1);

    derive_lst.add("forcez",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERFORCEZ,ca_derforcez),the_same_box);
    derive_lst.addComponent("forcez",desc_lst,State_Type,Density,1);
#endif

    derive_lst.add("magvel",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERMAGVEL,ca_dermagvel),the_same_box);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Xmom,BL_SPACEDIM);

#if (BL_SPACEDIM > 1)
    derive_lst.add("radvel",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERRADIALVEL,ca_derradialvel),the_same_box);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Xmom,BL_SPACEDIM);
#endif

    derive_lst.add("magmom",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERMAGMOM,ca_dermagmom),the_same_box);
    derive_lst.addComponent("magmom",desc_lst,State_Type,Xmom,BL_SPACEDIM);

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
    derive_lst.add("maggrav",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERMAGGRAV,ca_dermaggrav),the_same_box);
    derive_lst.addComponent("maggrav",desc_lst,Gravity_Type,0,BL_SPACEDIM);
#endif
#endif

#ifdef PARTICLES
    //
    // We want a derived type that corresponds to the number of particles
    // in each cell.  We only intend to use it in plotfiles for debugging
    // purposes.  We'll just use the DERNULL since don't do anything in
    // fortran for now.  We'll actually set the values in writePlotFile().
    //
    derive_lst.add("particle_count",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERNULL,ca_dernull),the_same_box);
    derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);

    derive_lst.add("particle_mass_density",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERNULL,ca_dernull),grow_box_by_one);
    derive_lst.addComponent("particle_mass_density",desc_lst,State_Type,Density,1);

    derive_lst.add("total_particle_count",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERNULL,ca_dernull),the_same_box);
    derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);
#endif

#ifdef RADIATION
    if (Radiation::do_multigroup) {
      derive_lst.add("Ertot", IndexType::TheCellType(),1,
		     BL_FORT_PROC_CALL(CA_DERERTOT,ca_derertot),the_same_box);
      derive_lst.addComponent("Ertot",desc_lst,Rad_Type,0,Radiation::nGroups);
    }
#endif

#ifdef NEUTRINO
    if (Radiation::nNeutrinoSpecies > 0 &&
        Radiation::plot_neutrino_group_energies_per_MeV) {
      char rad_name[10];
      int indx = 0;
      for (int j = 0; j < Radiation::nNeutrinoSpecies; j++) {
        for (int i = 0; i < Radiation::nNeutrinoGroups[j]; i++) {
          sprintf(rad_name, "Neuts%dg%d", j, i);
          derive_lst.add(rad_name,IndexType::TheCellType(),1,
                         BL_FORT_PROC_CALL(CA_DERNEUT,ca_derneut),
                         the_same_box);
          derive_lst.addComponent(rad_name,desc_lst,Rad_Type,indx,1);
          indx++;
        }
      }
    }

    if (Radiation::nNeutrinoSpecies > 0 &&
        Radiation::nNeutrinoGroups[0] > 0) {
#if 1
      derive_lst.add("Enue", IndexType::TheCellType(),1,
		     BL_FORT_PROC_CALL(CA_DERENUE,ca_derenue),the_same_box);
      derive_lst.addComponent("Enue",desc_lst,Rad_Type,0,Radiation::nGroups);
      derive_lst.add("Enuae", IndexType::TheCellType(),1,
		     BL_FORT_PROC_CALL(CA_DERENUAE,ca_derenuae),the_same_box);
      derive_lst.addComponent("Enuae",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // rho_Yl = rho(Ye + Ynue - Ynuebar)
      //
      derive_lst.add("rho_Yl",IndexType::TheCellType(),1,
                     BL_FORT_PROC_CALL(CA_DERRHOYL,ca_derrhoyl),the_same_box);
      // Don't actually need density for rho * Yl
      derive_lst.addComponent("rho_Yl",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("rho_Yl",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("rho_Yl",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Yl = (Ye + Ynue - Ynuebar)
      //
      derive_lst.add("Yl",IndexType::TheCellType(),1,
                     BL_FORT_PROC_CALL(CA_DERYL,ca_deryl),the_same_box);
      derive_lst.addComponent("Yl",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Yl",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Yl",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Ynue
      //
      derive_lst.add("Ynue",IndexType::TheCellType(),1,
                     BL_FORT_PROC_CALL(CA_DERYNUE,ca_derynue),the_same_box);
      derive_lst.addComponent("Ynue",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Ynue",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Ynue",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Ynuebar
      //
      derive_lst.add("Ynuae",IndexType::TheCellType(),1,
                     BL_FORT_PROC_CALL(CA_DERYNUAE,ca_derynuae),the_same_box);
      derive_lst.addComponent("Ynuae",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Ynuae",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Ynuae",desc_lst,Rad_Type,0,Radiation::nGroups);
#endif
    }
#endif

    for (int i = 0; i < NumSpec; i++)  {
      string spec_string(spec_names[i]);
      derive_lst.add(spec_string,IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,FirstSpec+i,1);
    }

    for (int i = 0; i < NumAux; i++)  {
      string aux_string(aux_names[i]);
      derive_lst.add(aux_string,IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(aux_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(aux_string,desc_lst,State_Type,FirstAux+i,1);
    }

#if 0
    //
    // A derived quantity equal to all the state variables.
    //
    derive_lst.add("FULLSTATE",IndexType::TheCellType(),NUM_STATE,FORT_DERCOPY,the_same_box);
    derive_lst.addComponent("FULLSTATE",desc_lst,State_Type,Density,NUM_STATE);

#endif


    // 
    // Problem-specific adds
#include <Problem_Derives.H>

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    ErrorSetUp();

    for (int i = 0; i < NumSpec; i++) {
      delete[] spec_names[i];
    }

    for (int i = 0; i < NumAux; i++) {
      delete[] aux_names[i];
    }
}
