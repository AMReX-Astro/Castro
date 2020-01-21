
    //
    // Analytic solution and error
    //
    derive_lst.add("Texact",IndexType::TheCellType(),1,ca_dertexact,the_same_box);
    derive_lst.addComponent("Texact",desc_lst,State_Type,Temp,1);

    derive_lst.add("Terror",IndexType::TheCellType(),1,ca_derterror,the_same_box);
    derive_lst.addComponent("Terror",desc_lst,State_Type,Temp,1);
