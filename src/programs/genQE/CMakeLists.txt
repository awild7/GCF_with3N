add_executable(genQE genQE.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc QEGenerator.cc)
target_link_libraries(genQE ${ROOT_LIBRARIES})

add_executable(genQE_EIC genQE_EIC.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc QEGenerator.cc)
target_link_libraries(genQE_EIC ${ROOT_LIBRARIES})
