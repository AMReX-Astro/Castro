name: Sedov

on: [pull_request]
jobs:
  Sedov-3d:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
        with:
          fetch-depth: 0

      - name: Get AMReX
        run: |
          mkdir external
          cd external
          git clone https://github.com/AMReX-Codes/amrex.git
          cd amrex
          git checkout development
          echo 'AMREX_HOME=$(GITHUB_WORKSPACE)/external/amrex' >> $GITHUB_ENV
          echo $AMREX_HOME
          if [[ -n "${AMREX_HOME}" ]]; then exit 1; fi
          cd ../..

      - name: Get Microphysics
        run: |
          cd external
          git clone https://github.com/AMReX-Codes/Microphysics.git
          cd Microphysics
          git checkout development
          echo 'MICROPHYSICS_HOME=$(GITHUB_WORKSPACE)/external/Microphysics' >> $GITHUB_ENV
          echo $MICROPHYSICS_HOME
          if [[ -n "${MICROPHYSICS_HOME}" ]]; then exit 1; fi
          cd ../..

      - name: Install dependencies
        run: |
          sudo apt-get update -y -qq
          sudo apt-get -qq -y install curl cmake jq clang g++>=9.3.0 gfortran>=9.3.0

      - name: Compile Sedov
        run: |
          cd Exec/hydro_tests/Sedov
          make DEBUG=TRUE USE_MPI=FALSE -j 4

      - name: Run Sedov
        run: |
          cd Exec/hydro_tests/Sedov
          ./Castro3d.gnu.DEBUG.ex inputs.3d.sph max_step=20 amr.max_level=2 amr.plot_files_output=0 amr.checkpoint_files_output=0

      - name: Compare to stored output
        run: |
          cd Exec/hydro_tests/Sedov
          diff grid_diag.out ci-benchmarks/grid_diag.out


