
//
// Analytic solution
//
derive_lst.add("rad_analytic", IndexType::TheCellType(), Radiation::nGroups,
               ca_deranalytic, the_same_box);
derive_lst.addComponent("rad_analytic", desc_lst, Rad_Type, 0,
                        Radiation::nGroups);
