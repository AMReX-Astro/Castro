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

#include "buildInfo.H"

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

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)    
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

void
Castro::variableSetUp ()
{

  // Castro::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job


  // initialize the start time for our CPU-time tracker
  startCPUTime = ParallelDescriptor::second();


    // Output the git commit hashes used to build the executable.

    if (ParallelDescriptor::IOProcessor()) {

	const char* castro_hash  = buildInfoGetGitHash(1);
	const char* boxlib_hash  = buildInfoGetGitHash(2);
	const char* buildgithash = buildInfoGetBuildGitHash();
	const char* buildgitname = buildInfoGetBuildGitName();

	if (strlen(castro_hash) > 0) {
	  std::cout << "\n" << "Castro git hash: " << castro_hash << "\n";
	}
	if (strlen(boxlib_hash) > 0) {
	  std::cout << "BoxLib git hash: " << boxlib_hash << "\n";
	}
	if (strlen(buildgithash) > 0){
	  std::cout << buildgitname << " git hash: " << buildgithash << "\n";
	}

	std::cout << "\n";
    }

    BL_ASSERT(desc_lst.size() == 0);

    // Initialize the network
    network_init();

    // Get options, set phys_bc
    read_params();

    // Initialize the runtime parameters for any of the external
    // microphysics
    extern_init();

    //
    // Set number of state variables and pointers to components
    //
    int use_sgs = 0;

    int cnt = 0;
    Density = cnt++;
    Xmom = cnt++;
    Ymom = cnt++;
    Zmom = cnt++;
    Eden = cnt++;
    Eint = cnt++;
#ifdef SGS
    Esgs = cnt++;
    use_sgs = 1;
#endif
    Temp = cnt++;

    NumAdv = 0;

    if (NumAdv > 0)
    {
        FirstAdv = cnt;
        cnt += NumAdv;
    }

    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    get_num_spec(&NumSpec);

    if (NumSpec > 0)
    {
        FirstSpec = cnt;
        cnt += NumSpec;
    }

    // Get the number of auxiliary quantities from the network model.
    get_num_aux(&NumAux);

    if (NumAux > 0)
    {
        FirstAux = cnt;
        cnt += NumAux;
    }

    NUM_STATE = cnt;

    // Define NUM_GROW from the f90 module.
    get_method_params(&NUM_GROW);

    const Real run_strt = ParallelDescriptor::second() ; 

#ifndef DIFFUSION
    static Real diffuse_cutoff_density = -1.e200;
#endif

    // we want const_grav in F90, get it here from parmparse, since it
    // it not in the Castro namespace
    ParmParse pp("gravity");
    Real const_grav = 0;
    pp.query("const_grav", const_grav);

    // Pass in the name of the gravity type we're using.
    std::string gravity_type = "none";
    pp.query("gravity_type", gravity_type);    
    int gravity_type_length = gravity_type.length();
    Array<int> gravity_type_name(gravity_type_length);

    for (int i = 0; i < gravity_type_length; i++)
      gravity_type_name[i] = gravity_type[i];    

    int get_g_from_phi = 0;
    pp.query("get_g_from_phi", get_g_from_phi);
    
    set_method_params(dm, Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux, 
		      NumAdv, 
		      gravity_type_name.dataPtr(), &gravity_type_length,
		      get_g_from_phi,
		      use_sgs,
		      diffuse_cutoff_density,
		      const_grav);

#include <castro_call_set_meth.H>

    Real run_stop = ParallelDescriptor::second() - run_strt;
 
    ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
 
    if (ParallelDescriptor::IOProcessor())
        std::cout << "\nTime in set_method_params: " << run_stop << '\n' ;

    int coord_type = Geometry::Coord();

    // Get the center variable from the inputs and pass it directly to Fortran.
    Array<Real> center(BL_SPACEDIM, 0.0);
    ParmParse ppc("castro");
    ppc.queryarr("center",center,0,BL_SPACEDIM);

    set_problem_params(dm,phys_bc.lo(),phys_bc.hi(),
		       Interior,Inflow,Outflow,Symmetry,SlipWall,NoSlipWall,coord_type,
		       Geometry::ProbLo(),Geometry::ProbHi(),center.dataPtr());

    // Read in the parameters for the tagging criteria
    // and store them in the Fortran module.

    int probin_file_length = probin_file.length();
    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
      probin_file_name[i] = probin_file[i];

    get_tagging_params(probin_file_name.dataPtr(),&probin_file_length);

    // Read in the parameters for the sponge
    // and store them in the Fortran module.
    
    get_sponge_params(probin_file_name.dataPtr(),&probin_file_length);    

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
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,NUM_GROW,3,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

    // Source terms. Currently this holds dS/dt for each of the NVAR state variables.

    store_in_checkpoint = true;
    desc_lst.addDescriptor(Source_Type, IndexType::TheCellType(),
			   StateDescriptor::Point,NUM_GROW,NUM_STATE,
			   &cell_cons_interp, state_data_extrap,store_in_checkpoint);
    
#ifdef ROTATION
    store_in_checkpoint = false;
    desc_lst.addDescriptor(PhiRot_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Rotation_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,NUM_GROW,3,
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
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
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
    std::vector<std::string> spec_names;
    for (int i = 0; i < NumSpec; i++) {
          int len = 20;
          Array<int> int_spec_names(len);
          // This call return the actual length of each string in "len" 
          get_spec_names(int_spec_names.dataPtr(),&i,&len);
          char char_spec_names[len+1];
          for (int j = 0; j < len; j++) 
             char_spec_names[j] = int_spec_names[j];
          char_spec_names[len] = '\0';
	  spec_names.push_back(char_spec_names);
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
        name[cnt] = "rho_" + spec_names[i];
    }

    // Get the auxiliary names from the network model.
    std::vector<std::string> aux_names;
    for (int i = 0; i < NumAux; i++) {
          int len = 20;
          Array<int> int_aux_names(len);
          // This call return the actual length of each string in "len"
          get_aux_names(int_aux_names.dataPtr(),&i,&len);
	  char char_aux_names[len+1];
          for (int j = 0; j < len; j++)
             char_aux_names[j] = int_aux_names[j];
          char_aux_names[len] = '\0';
	  aux_names.push_back(char_aux_names);
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
        name[cnt] = "rho_" + aux_names[i];
    }

    desc_lst.setComponent(State_Type,
                          Density,
                          name,
                          bcs,
                          BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL,ca_denfill),
                                    BL_FORT_PROC_CALL(CA_HYPFILL,ca_hypfill)));

#ifdef GRAVITY
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(PhiGrav_Type,0,"phiGrav",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_PHIGRAVFILL,ca_phigravfill)));
    set_x_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Gravity_Type,0,"grav_x",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVXFILL,ca_gravxfill)));
    set_y_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Gravity_Type,1,"grav_y",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVYFILL,ca_gravyfill)));
    set_z_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Gravity_Type,2,"grav_z",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVZFILL,ca_gravzfill)));
#endif

    // For rotation we'll use the same boundary condition routines as for gravity, 
    // since we use the rotation in the same manner as in the gravity.

#ifdef ROTATION
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(PhiRot_Type,0,"phiRot",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_PHIGRAVFILL,ca_phigravfill)));
    set_x_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Rotation_Type,0,"rot_x",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVXFILL,ca_gravxfill)));
    set_y_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Rotation_Type,1,"rot_y",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVYFILL,ca_gravyfill)));
    set_z_vel_bc(bc,phys_bc);
    desc_lst.setComponent(Rotation_Type,2,"rot_z",bc,
			  BndryFunc(BL_FORT_PROC_CALL(CA_GRAVZFILL,ca_gravzfill)));
#endif

    // Source term array will use standard hyperbolic fill.

    Array<std::string> sources_name(NUM_STATE);    
    
    for (int i = 0; i < NUM_STATE; i++)
      sources_name[i] = name[i] + "_source";
    
    desc_lst.setComponent(Source_Type,Density,sources_name,bcs,
			  BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL,ca_denfill),
				    BL_FORT_PROC_CALL(CA_HYPFILL,ca_hypfill)));       
    
#ifdef LEVELSET
    desc_lst.setComponent(LS_State_Type,0,"LSphi",bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_PHIFILL,ca_phifill)));
#endif

#ifdef REACTIONS
    std::string name_react;
    for (int i=0; i<NumSpec; ++i)
    {
       set_scalar_bc(bc,phys_bc);
       name_react = "omegadot_" + spec_names[i];
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
#endif

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure",IndexType::TheCellType(),1,ca_derpres,the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng",IndexType::TheCellType(),1,ca_derkineng,the_same_box);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Xmom,3);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed",IndexType::TheCellType(),1,ca_dersoundspeed,the_same_box);
    derive_lst.addComponent("soundspeed",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber",IndexType::TheCellType(),1,ca_dermachnumber,the_same_box);
    derive_lst.addComponent("MachNumber",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM == 1)
    //
    // Wave speed u+c
    //
    derive_lst.add("uplusc",IndexType::TheCellType(),1,ca_deruplusc,the_same_box);
    derive_lst.addComponent("uplusc",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Wave speed u-c
    //
    derive_lst.add("uminusc",IndexType::TheCellType(),1,ca_deruminusc,the_same_box);
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
    derive_lst.add("entropy",IndexType::TheCellType(),1,ca_derentropy,the_same_box);
    derive_lst.addComponent("entropy",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Vorticity
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,ca_dermagvort,grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom != Density+1)
       BoxLib::Error("We are assuming Xmom = Density + 1 in Castro_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,4);

    //
    // Div(u)
    //
    derive_lst.add("divu",IndexType::TheCellType(),1,ca_derdivu,grow_box_by_one);
    derive_lst.addComponent("divu",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("divu",desc_lst,State_Type,Xmom,3);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E",IndexType::TheCellType(),1,ca_dereint1,the_same_box);
    derive_lst.addComponent("eint_E",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e",IndexType::TheCellType(),1,ca_dereint2,the_same_box);
    derive_lst.addComponent("eint_e",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden",IndexType::TheCellType(),1,ca_derlogden,the_same_box);
    derive_lst.addComponent("logden",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("StateErr",IndexType::TheCellType(),3,ca_derstate,grow_box_by_one);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Temp,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,FirstSpec,1);

    //
    // X from rhoX
    //
    for (int i = 0; i < NumSpec; i++){
      std::string spec_string = "X("+spec_names[i]+")";
      derive_lst.add(spec_string,IndexType::TheCellType(),1,ca_derspec,the_same_box);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,FirstSpec+i,1);
    }
    //
    // Velocities
    //
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Xmom,1);

    derive_lst.add("y_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Ymom,1);

    derive_lst.add("z_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Zmom,1);

#ifdef SGS
    derive_lst.add("K",IndexType::TheCellType(),1,ca_dervel,the_same_box);
    derive_lst.addComponent("K",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("K",desc_lst,State_Type,Esgs,1);

    derive_lst.add("forcex",IndexType::TheCellType(),1,ca_derforcex,the_same_box);
    derive_lst.addComponent("forcex",desc_lst,State_Type,Density,1);

    derive_lst.add("forcey",IndexType::TheCellType(),1,ca_derforcey,the_same_box);
    derive_lst.addComponent("forcey",desc_lst,State_Type,Density,1);

    derive_lst.add("forcez",IndexType::TheCellType(),1,ca_derforcez,the_same_box);
    derive_lst.addComponent("forcez",desc_lst,State_Type,Density,1);
#endif

    derive_lst.add("magvel",IndexType::TheCellType(),1,ca_dermagvel,the_same_box);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Xmom,3);

    derive_lst.add("radvel",IndexType::TheCellType(),1,ca_derradialvel,the_same_box);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Xmom,3);

    derive_lst.add("magmom",IndexType::TheCellType(),1,ca_dermagmom,the_same_box);
    derive_lst.addComponent("magmom",desc_lst,State_Type,Xmom,3);

#ifdef GRAVITY
    derive_lst.add("maggrav",IndexType::TheCellType(),1,ca_dermaggrav,the_same_box);
    derive_lst.addComponent("maggrav",desc_lst,Gravity_Type,0,3);
#endif

#ifdef PARTICLES
    //
    // We want a derived type that corresponds to the number of particles
    // in each cell.  We only intend to use it in plotfiles for debugging
    // purposes.  We'll just use the DERNULL since don't do anything in
    // fortran for now.  We'll actually set the values in writePlotFile().
    //
    derive_lst.add("particle_count",IndexType::TheCellType(),1,ca_dernull,the_same_box);
    derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);

    derive_lst.add("particle_mass_density",IndexType::TheCellType(),1,ca_dernull,grow_box_by_one);
    derive_lst.addComponent("particle_mass_density",desc_lst,State_Type,Density,1);

    derive_lst.add("total_particle_count",IndexType::TheCellType(),1,ca_dernull,the_same_box);
    derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);
#endif

#ifdef RADIATION
    if (Radiation::do_multigroup) {
      derive_lst.add("Ertot", IndexType::TheCellType(),1,ca_derertot,the_same_box);
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
          derive_lst.add(rad_name,IndexType::TheCellType(),1,ca_derneut,the_same_box);
          derive_lst.addComponent(rad_name,desc_lst,Rad_Type,indx,1);
          indx++;
        }
      }
    }

    if (Radiation::nNeutrinoSpecies > 0 &&
        Radiation::nNeutrinoGroups[0] > 0) {
      derive_lst.add("Enue", IndexType::TheCellType(),1,ca_derenue,the_same_box);
      derive_lst.addComponent("Enue",desc_lst,Rad_Type,0,Radiation::nGroups);
      derive_lst.add("Enuae", IndexType::TheCellType(),1,ca_derenuae,the_same_box);
      derive_lst.addComponent("Enuae",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // rho_Yl = rho(Ye + Ynue - Ynuebar)
      //
      derive_lst.add("rho_Yl",IndexType::TheCellType(),1,ca_derrhoyl,the_same_box);
      // Don't actually need density for rho * Yl
      derive_lst.addComponent("rho_Yl",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("rho_Yl",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("rho_Yl",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Yl = (Ye + Ynue - Ynuebar)
      //
      derive_lst.add("Yl",IndexType::TheCellType(),1,ca_deryl,the_same_box);
      derive_lst.addComponent("Yl",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Yl",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Yl",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Ynue
      //
      derive_lst.add("Ynue",IndexType::TheCellType(),1,ca_derynue,the_same_box);
      derive_lst.addComponent("Ynue",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Ynue",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Ynue",desc_lst,Rad_Type,0,Radiation::nGroups);
      //
      // Ynuebar
      //
      derive_lst.add("Ynuae",IndexType::TheCellType(),1,ca_derynuae,the_same_box);
      derive_lst.addComponent("Ynuae",desc_lst,State_Type,Density,1);
      // FirstAux is (rho * Ye)
      derive_lst.addComponent("Ynuae",desc_lst,State_Type,FirstAux,1);
      derive_lst.addComponent("Ynuae",desc_lst,Rad_Type,0,Radiation::nGroups);
    }
#endif

    for (int i = 0; i < NumSpec; i++)  {
      derive_lst.add(spec_names[i],IndexType::TheCellType(),1,ca_derspec,the_same_box);
      derive_lst.addComponent(spec_names[i],desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_names[i],desc_lst,State_Type,FirstSpec+i,1);
    }

    for (int i = 0; i < NumAux; i++)  {
      derive_lst.add(aux_names[i],IndexType::TheCellType(),1,ca_derspec,the_same_box);
      derive_lst.addComponent(aux_names[i],desc_lst,State_Type,Density,1);
      derive_lst.addComponent(aux_names[i],desc_lst,State_Type,FirstAux+i,1);
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

}
