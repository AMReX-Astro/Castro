    //
    // X(ash) from rhoX - sum of all mass fractions for the elements beyond oxygen,
    // excluding the species that make up the star
    //
    derive_lst.add("X(ash)",IndexType::TheCellType(),1,ca_derxash,the_same_box);
    derive_lst.addComponent("X(ash)",desc_lst,State_Type,URHO,1);
    derive_lst.addComponent("X(ash)",desc_lst,State_Type,UFS,NumSpec);
