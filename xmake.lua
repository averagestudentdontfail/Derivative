set_project("MFEM_Example8")
set_version("1.0.0")
set_languages("c++17")

-- Spack
local mfem_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/mfem-4.7.0-rsop3l7ix26ze4hyqvyaitixboqdko7g"
local mpi_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/openmpi-5.0.6-b2hzcc2nbkd4z66hutm3ood6blytv55b"
local hypre_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/hypre-2.32.0-fwr5lltmbnhnkwdlvnxtvchj7anp52qn"
local sundials_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/sundials-7.2.1-vxaqjkeh4irnfm62a2lwsaf7mbkxbdek"
local petsc_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/petsc-3.22.4-hqk7qbo3du7qp53l2xps2ibmdmweyx7q"
local conduit_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/conduit-0.9.3-d3fvlbb4psq4q35r2wumgphevcotejm7"
local eigen_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/eigen-3.4.0-enqwlmq352incu7he6lkthjyvgs2st3s"
local metis_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/metis-5.1.0-hhs4j24hpqjmbtnwdglcyeca6vgmmr2b"
local superlu_dist_dir = "/home/send2/spack/opt/spack/linux-ubuntu24.04-x86_64_v4/gcc-13.3.0/superlu-dist-9.1.0-mz6dtycjksxyh6clwbhnxbe27b67gqxi"

-- Target
target("ex8")
    set_kind("binary")
    add_files("src/ex8.cpp")

    -- Directories
    add_includedirs("/usr/include")
    add_includedirs("/usr/include/suitesparse")
    add_includedirs(mfem_dir .. "/include")
    add_includedirs(mpi_dir .. "/include")
    add_includedirs(hypre_dir .. "/include")
    add_includedirs(sundials_dir .. "/include")
    add_includedirs(petsc_dir .. "/include")
    add_includedirs(conduit_dir .. "/include/conduit")
    add_includedirs(eigen_dir .. "/include/eigen3")
    add_includedirs("src")

    -- Directories
    add_linkdirs("/usr/lib/x86_64-linux-gnu")
    add_linkdirs(mfem_dir .. "/lib")
    add_linkdirs(mpi_dir .. "/lib")
    add_linkdirs(hypre_dir .. "/lib")
    add_linkdirs(sundials_dir .. "/lib")
    add_linkdirs(petsc_dir .. "/lib")
    add_linkdirs(conduit_dir .. "/lib")
    add_linkdirs(metis_dir .. "/lib")
    add_linkdirs(superlu_dist_dir .. "/lib")

    -- Libraries
    add_links(
        "umfpack", "cholmod", "amd", "colamd", "suitesparseconfig"
    )
    add_links(
        "HYPRE", "metis", "petsc", "superlu_dist", 
        "sundials_cvode", "sundials_nvecserial", "sundials_nvecparallel", 
        "conduit", "conduit_relay", "conduit_blueprint", 
        "lapack", "blas", "z", "mpi" -
    )
    add_links("mfem") 

    -- Compiler Flags
    add_cxflags("-Wall", "-fopenmp")
    if is_mode("debug") then
        add_cxflags("-g", "-O0")
        add_defines("DEBUG")
    else
        add_cxflags("-O3")
        add_defines("NDEBUG")
    end

    -- Linker Flags
    add_ldflags("-fopenmp")
    local run_lib_paths = {
        mfem_dir .. "/lib",
        mpi_dir .. "/lib",
        hypre_dir .. "/lib",
        sundials_dir .. "/lib",
        petsc_dir .. "/lib",
        conduit_dir .. "/lib",
        metis_dir .. "/lib",
        superlu_dist_dir .. "/lib"
    }

    add_runenvs("LD_LIBRARY_PATH", table.concat(run_lib_paths, ":") .. ":${LD_LIBRARY_PATH}")