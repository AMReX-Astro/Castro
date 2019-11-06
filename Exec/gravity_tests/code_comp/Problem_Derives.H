// problem-specific derived variables would be put here
derive_lst.add("external_heating",IndexType::TheCellType(),1,ca_derextheating,the_same_box);
derive_lst.addComponent("external_heating",desc_lst,State_Type,Density,NUM_STATE);

derive_lst.add("adiabatic_invariant",IndexType::TheCellType(),1,ca_deradinvariant,the_same_box);
derive_lst.addComponent("adiabatic_invariant",desc_lst,State_Type,Density,NUM_STATE);

derive_lst.add("enthalpy_flux_y",IndexType::TheCellType(),1,ca_derenthalpyfluxy,the_same_box);
derive_lst.addComponent("enthalpy_flux_y",desc_lst,State_Type,Density,NUM_STATE);

derive_lst.add("kineng_flux_y",IndexType::TheCellType(),1,ca_derkinengfluxy,the_same_box);
// Here we exploit the fact that Xmom = Density + 1
//   in order to use the correct interpolation.
if (Xmom != Density+1)
    amrex::Error("We are assuming Xmom = Density + 1 in Castro_setup.cpp");
derive_lst.addComponent("kineng_flux_y",desc_lst,State_Type,Density,4);
