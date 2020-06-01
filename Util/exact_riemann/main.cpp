#include <iostream>
#include <cstring>
#include <vector>

#include "riemann_F.H"

int main(int argc, char *argv[]) {

  std::cout << "starting the exact Riemann solver..." << std::endl;
  std::cout << argv[1] << std::endl;
  std::cout << strlen(argv[1]) << std::endl;

  const int probin_file_length = strlen(argv[1]);
  std::vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = argv[1][i];

  std::cout << "here" << std::endl;

  ca_extern_init(probin_file_name.data(), &probin_file_length);

  riemann_exact();

}
