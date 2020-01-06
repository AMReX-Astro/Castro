jobInfoFile << (RadSolve::level_solver_flag == 1 ? "    " : "[*] ") << "radsolve.level_solver_flag = " << RadSolve::level_solver_flag << std::endl;
jobInfoFile << (RadSolve::use_hypre_nonsymmetric_terms == 0 ? "    " : "[*] ") << "radsolve.use_hypre_nonsymmetric_terms = " << RadSolve::use_hypre_nonsymmetric_terms << std::endl;
jobInfoFile << (RadSolve::reltol == 1.e-10 ? "    " : "[*] ") << "radsolve.reltol = " << RadSolve::reltol << std::endl;
jobInfoFile << (RadSolve::abstol == 1.e-10 ? "    " : "[*] ") << "radsolve.abstol = " << RadSolve::abstol << std::endl;
jobInfoFile << (RadSolve::maxiter == 40 ? "    " : "[*] ") << "radsolve.maxiter = " << RadSolve::maxiter << std::endl;
jobInfoFile << (RadSolve::alpha == 1.0 ? "    " : "[*] ") << "radsolve.alpha = " << RadSolve::alpha << std::endl;
jobInfoFile << (RadSolve::beta == 1.0 ? "    " : "[*] ") << "radsolve.beta = " << RadSolve::beta << std::endl;
jobInfoFile << (RadSolve::verbose == 0 ? "    " : "[*] ") << "radsolve.verbose = " << RadSolve::verbose << std::endl;
