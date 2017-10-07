
#ifdef ROTATION
    //
    // Momemtum
    //
    derive_lst.add("inertial_momentum_x",IndexType::TheCellType(),1,ca_derinertialmomentumx,the_same_box);
    derive_lst.addComponent("inertial_momentum_x",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_x",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_momentum_y",IndexType::TheCellType(),1,ca_derinertialmomentumy,the_same_box);
    derive_lst.addComponent("inertial_momentum_y",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_y",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_momentum_z",IndexType::TheCellType(),1,ca_derinertialmomentumz,the_same_box);
    derive_lst.addComponent("inertial_momentum_z",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_momentum_z",desc_lst,State_Type,Xmom,3);

    //
    // Angular momentum
    //
    derive_lst.add("inertial_angular_momentum_x",IndexType::TheCellType(),1,ca_derinertialangmomx,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_x",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_x",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_angular_momentum_y",IndexType::TheCellType(),1,ca_derinertialangmomy,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_y",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_y",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_angular_momentum_z",IndexType::TheCellType(),1,ca_derinertialangmomz,the_same_box);
    derive_lst.addComponent("inertial_angular_momentum_z",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_angular_momentum_z",desc_lst,State_Type,Xmom,3);

    //
    // Radial momentum
    //
    derive_lst.add("inertial_radial_momentum_x",IndexType::TheCellType(),1,ca_derinertialradmomx,the_same_box);
    derive_lst.addComponent("inertial_radial_momentum_x",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_radial_momentum_x",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_radial_momentum_y",IndexType::TheCellType(),1,ca_derinertialradmomy,the_same_box);
    derive_lst.addComponent("inertial_radial_momentum_y",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_radial_momentum_y",desc_lst,State_Type,Xmom,3);

    derive_lst.add("inertial_radial_momentum_z",IndexType::TheCellType(),1,ca_derinertialradmomz,the_same_box);
    derive_lst.addComponent("inertial_radial_momentum_z",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("inertial_radial_momentum_z",desc_lst,State_Type,Xmom,3);
#endif

#ifdef GRAVITY
#ifdef ROTATION
    //
    // Effective potential (gravitational + rotational)
    //
    derive_lst.add("phiEff",IndexType::TheCellType(),1,ca_derphieff,the_same_box);
    derive_lst.addComponent("phiEff",desc_lst,PhiGrav_Type,0,1);
    derive_lst.addComponent("phiEff",desc_lst,PhiRot_Type,0,1);

    //
    // Effective potential from a point-mass approximation to the potential from the primary
    //
    derive_lst.add("phiEffPM_P",IndexType::TheCellType(),1,ca_derphieffpm_p,the_same_box);
    derive_lst.addComponent("phiEffPM_P",desc_lst,PhiRot_Type,0,1);

    //
    // Effective potential from a point-mass approximation to the potential from the secondary
    //
    derive_lst.add("phiEffPM_S",IndexType::TheCellType(),1,ca_derphieffpm_s,the_same_box);
    derive_lst.addComponent("phiEffPM_S",desc_lst,PhiRot_Type,0,1);
#endif
#endif

#ifdef GRAVITY
    //
    // rho * phi
    //
    derive_lst.add("rho_phiGrav",IndexType::TheCellType(),1,ca_derrhophiGrav,the_same_box);
    derive_lst.addComponent("rho_phiGrav",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("rho_phiGrav",desc_lst,PhiGrav_Type,0,1);
#endif

#ifdef ROTATION
    //
    // rho * phi
    //
    derive_lst.add("rho_phiRot",IndexType::TheCellType(),1,ca_derrhophiRot,the_same_box);
    derive_lst.addComponent("rho_phiRot",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("rho_phiRot",desc_lst,Rotation_Type,0,1);
#endif

#ifdef GRAVITY
#ifdef ROTATION
    //
    // Mask determining which zones are inside the primary
    //
    derive_lst.add("primarymask",IndexType::TheCellType(),1,ca_derprimarymask,the_same_box);
    derive_lst.addComponent("primarymask",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("primarymask",desc_lst,PhiRot_Type,0,1);


    //
    // Mask determining which zones are inside the secondary
    //
    derive_lst.add("secondarymask",IndexType::TheCellType(),1,ca_dersecondarymask,the_same_box);
    derive_lst.addComponent("secondarymask",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("secondarymask",desc_lst,PhiRot_Type,0,1);
#endif
#endif

#ifdef REACTIONS
    //
    // Critical radius for spontaneous detonation
    //
    derive_lst.add("ignition_radius",IndexType::TheCellType(),1,ca_derignitionradius,grow_box_by_one);
    derive_lst.addComponent("ignition_radius",desc_lst,State_Type,Density,NUM_STATE);
#endif
