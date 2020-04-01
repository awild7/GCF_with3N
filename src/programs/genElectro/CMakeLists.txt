add_executable(genElectro genElectro.cc ${LIBRARIES}/cross_sections/photoCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc electroGenerator.cc)
target_link_libraries(genElectro ${ROOT_LIBRARIES})
