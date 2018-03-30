jobInfoFile << (Castro::particle_verbose == 0 ? "    " : "[*] ") << "particles.particle_verbose = " << Castro::particle_verbose << std::endl;
jobInfoFile << (Castro::particle_init_file == "" ? "    " : "[*] ") << "particles.particle_init_file = " << Castro::particle_init_file << std::endl;
jobInfoFile << (Castro::particle_restart_file == "" ? "    " : "[*] ") << "particles.particle_restart_file = " << Castro::particle_restart_file << std::endl;
jobInfoFile << (Castro::restart_from_nonparticle_chkfile == 0 ? "    " : "[*] ") << "particles.restart_from_nonparticle_chkfile = " << Castro::restart_from_nonparticle_chkfile << std::endl;
jobInfoFile << (Castro::particle_output_file == "" ? "    " : "[*] ") << "particles.particle_output_file = " << Castro::particle_output_file << std::endl;
jobInfoFile << (Castro::timestamp_dir == "" ? "    " : "[*] ") << "particles.timestamp_dir = " << Castro::timestamp_dir << std::endl;
jobInfoFile << (Castro::timestamp_density == 1 ? "    " : "[*] ") << "particles.timestamp_density = " << Castro::timestamp_density << std::endl;
jobInfoFile << (Castro::timestamp_temperature == 0 ? "    " : "[*] ") << "particles.timestamp_temperature = " << Castro::timestamp_temperature << std::endl;
