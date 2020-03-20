add_executable(genIon genIon.cc ${LIBRARIES}/cross_sections/NNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ionGenerator.cc)
target_link_libraries(genIon ${ROOT_LIBRARIES})
