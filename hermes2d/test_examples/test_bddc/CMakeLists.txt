project(test_bddc)

add_executable(${PROJECT_NAME} main.cpp bddcml_wrapper.cpp definitions.cpp)

add_definitions(-DAdd_)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS} PROPERTY LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})

FIND_PACKAGE(MPI REQUIRED)

set(BDDCML_ROOT "/home/lukas/hpfem/bddcml")

find_package(BDDCML REQUIRED)

include_directories(${BDDCML_INCLUDES})

add_custom_command(TARGET ${PROJECT_NAME}
      POST_BUILD
			COMMAND "/usr/bin/c++   -Wuninitialized -Wvla -Wsign-compare -fopenmp -w     CMakeFiles/test_bddc.dir/main.cpp.o CMakeFiles/test_bddc.dir/bddcml_wrapper.cpp.o CMakeFiles/test_bddc.dir/definitions.cpp.o  -o test_bddc -rdynamic /usr/lib/openmpi/lib/libmpi_cxx.so /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libopen-rte.so /usr/lib/openmpi/lib/libopen-pal.so -ldl -lnsl -lutil -lm -ldl /usr/lib/openmpi/lib/libmpi_f90.so /usr/lib/openmpi/lib/libmpi_f77.so /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libopen-rte.so /usr/lib/openmpi/lib/libopen-pal.so -ldl -lnsl -lutil -lm -ldl /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libopen-rte.so /usr/lib/openmpi/lib/libopen-pal.so -ldl -lnsl -lutil -lm -ldl ../../libhermes2d.so /home/lukas/hpfem/bddcml/lib/libbddcml.a -ldmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi -lblacs-openmpi -llapack -lblas /usr/lib/openmpi/lib/libmpi_f77.so /home/lukas/hpfem/bddcml/lib/libbddcml.a -ldmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi -lblacs-openmpi /home/lukas/blopex-read-only/blopex_serial_double/multivector/multi_vector.o /home/lukas/blopex-read-only/blopex_serial_double/pcg_multi/pcg_multi.o /home/lukas/blopex-read-only/blopex_serial_double/matmultivec/matmultivec.o /home/lukas/blopex-read-only/blopex_abstract/lib/libBLOPEX.a ../../../hermes_common/libhermes_common.so -lumfpack -lamd -llapack -lblas -lpthread -lrt -lglut -lGLEW -lxerces-c -Wl,-rpath,/usr/lib/openmpi/lib:/home/lukas/hpfem/hermes/l-korous/hermes/hermes2d:/home/lukas/hpfem/hermes/l-korous/hermes/hermes_common")
