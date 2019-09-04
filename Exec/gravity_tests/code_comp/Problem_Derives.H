// problem-specific derived variables would be put here
derive_lst.add("external_heating",IndexType::TheCellType(),1,ca_derextheating,the_same_box);
derive_lst.addComponent("external_heating",desc_lst,State_Type,Density,NUM_STATE);
