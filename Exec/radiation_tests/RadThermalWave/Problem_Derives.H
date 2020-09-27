
    //
    // Analytic solution and error
    //
    derive_lst.add("Texact",IndexType::TheCellType(),1,dertexact,the_same_box);
    derive_lst.addComponent("Texact",desc_lst,State_Type,UTEMP,1);

    derive_lst.add("Terror",IndexType::TheCellType(),1,derterror,the_same_box);
    derive_lst.addComponent("Terror",desc_lst,State_Type,UTEMP,1);
