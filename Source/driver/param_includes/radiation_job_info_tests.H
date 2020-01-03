jobInfoFile << (Radiation::do_real_eos == 1 ? "    " : "[*] ") << "radiation.do_real_eos = " << Radiation::do_real_eos << std::endl;
jobInfoFile << (Radiation::const_c_v == -1.0 ? "    " : "[*] ") << "radiation.const_c_v = " << Radiation::const_c_v << std::endl;
jobInfoFile << (Radiation::c_v_exp_m == 0.0 ? "    " : "[*] ") << "radiation.c_v_exp_m = " << Radiation::c_v_exp_m << std::endl;
jobInfoFile << (Radiation::c_v_exp_n == 0.0 ? "    " : "[*] ") << "radiation.c_v_exp_n = " << Radiation::c_v_exp_n << std::endl;
jobInfoFile << (Radiation::prop_temp_floor == 0.0 ? "    " : "[*] ") << "radiation.prop_temp_floor = " << Radiation::prop_temp_floor << std::endl;
