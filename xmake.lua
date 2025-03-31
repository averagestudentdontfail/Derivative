set_project("Derivative")
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
 
-- Directories
add_includedirs("/usr/include")
add_includedirs("/usr/include/suitesparse")

add_linkdirs("/usr/lib/x86_64-linux-gnu")

add_includedirs(mfem_dir .. "/include")
add_includedirs(mpi_dir .. "/include")
add_includedirs(hypre_dir .. "/include")
add_includedirs(sundials_dir .. "/include")
add_includedirs(petsc_dir .. "/include")
add_includedirs(conduit_dir .. "/include/conduit")
add_includedirs(eigen_dir .. "/include/eigen3")


add_linkdirs(mfem_dir .. "/lib")
add_linkdirs(mpi_dir .. "/lib")
add_linkdirs(hypre_dir .. "/lib")
add_linkdirs(sundials_dir .. "/lib")
add_linkdirs(petsc_dir .. "/lib")
add_linkdirs(conduit_dir .. "/lib")
add_linkdirs(conduit_dir .. "/share")

-- Dependencies
add_links("mfem")
add_links("hypre")
add_links("metis")
add_links("petsc")
add_links("eigen")
add_links("eigen_dense")
add_links("superlu_dist")
add_links("sundials_cvode")
add_links("sundials_nvecserial")
add_links("sundials_nvecparallel")
add_links("umfpack")
add_links("cholmod")
add_links("conduit")
add_links("conduit_relay")
add_links("conduit_blueprint")
add_links("amd")
add_links("colamd")
add_links("suitesparseconfig")
add_links("lapack")
add_links("blas")
add_links("z")
add_links("mpi")

-- Target
target("Derivative")
    set_kind("binary")
    
    -- Source
    add_files("src/*.cpp")
    
    -- Header - This line was missing
    add_headerfiles("src/*.hpp")
    
    -- Directories
    add_includedirs("src")
    
    -- Compiler Flags
    add_cxflags("-Wall", "-fopenmp")
    if is_mode("debug") then
        add_cxflags("-g", "-O0")
    else
        add_cxflags("-O3")
    end
end