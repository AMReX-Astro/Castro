jobInfoFile << (Diffusion::verbose == 0 ? "    " : "[*] ") << "diffusion.verbose = " << Diffusion::verbose << std::endl;
jobInfoFile << (Diffusion::use_mlmg_solver == 1 ? "    " : "[*] ") << "diffusion.use_mlmg_solver = " << Diffusion::use_mlmg_solver << std::endl;
jobInfoFile << (Diffusion::mlmg_maxorder == 4 ? "    " : "[*] ") << "diffusion.mlmg_maxorder = " << Diffusion::mlmg_maxorder << std::endl;
