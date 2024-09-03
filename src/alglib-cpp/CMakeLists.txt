add_library(alglib STATIC 
    src/optimization.cpp 
    src/alglibinternal.cpp 
    src/ap.cpp 
    src/dataanalysis.cpp 
    src/alglibmisc.cpp 
    src/diffequations.cpp 
    src/fasttransforms.cpp 
    src/integration.cpp 
    src/interpolation.cpp 
    src/kernels_avx2.cpp 
    src/kernels_fma.cpp 
    src/kernels_sse2.cpp 
    src/linalg.cpp 
    src/solvers.cpp 
    src/specialfunctions.cpp 
    src/statistics.cpp
)
target_include_directories(alglib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(alglib PUBLIC cxx_std_20)