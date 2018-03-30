jobInfoFile << (particle_verbose == 0 ? "    " : "[*] ") << "particle_verbose = " << particle_verbose << std::endl;
jobInfoFile << (particle_init_file == "" ? "    " : "[*] ") << "particle_init_file = " << particle_init_file << std::endl;
jobInfoFile << (particle_restart_file == "" ? "    " : "[*] ") << "particle_restart_file = " << particle_restart_file << std::endl;
jobInfoFile << (restart_from_nonparticle_chkfile == 0 ? "    " : "[*] ") << "restart_from_nonparticle_chkfile = " << restart_from_nonparticle_chkfile << std::endl;
jobInfoFile << (particle_output_file == "" ? "    " : "[*] ") << "particle_output_file = " << particle_output_file << std::endl;
jobInfoFile << (timestamp_dir == "" ? "    " : "[*] ") << "timestamp_dir = " << timestamp_dir << std::endl;
jobInfoFile << (timestamp_density == 1 ? "    " : "[*] ") << "timestamp_density = " << timestamp_density << std::endl;
jobInfoFile << (timestamp_temperature == 0 ? "    " : "[*] ") << "timestamp_temperature = " << timestamp_temperature << std::endl;
