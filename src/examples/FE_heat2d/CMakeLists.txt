


add_executable("FE_heat2d_sdc" FE_heat2d_sdc.cpp)
target_link_dune_default_libraries("FE_heat2d_sdc")

add_executable("FE_heat2d_mlsdc" FE_heat2d_mlsdc.cpp)
target_link_dune_default_libraries("FE_heat2d_mlsdc")

add_executable("FE_heat2d_pfasst" FE_heat2d_pfasst.cpp)
target_link_dune_default_libraries("FE_heat2d_pfasst")
