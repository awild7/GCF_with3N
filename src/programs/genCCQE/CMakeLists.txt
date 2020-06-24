add_executable(genCCQE genCCQE.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc CCQEGenerator.cc)
target_link_libraries(genCCQE ${ROOT_LIBRARIES})

