add_executable(genQE_3N genQE_3N.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
target_link_libraries(genQE_3N ${ROOT_LIBRARIES})

#add_executable(read_in read_in.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
#target_link_libraries(read_in ${ROOT_LIBRARIES})
#
#add_executable(read_in_withclass read_in_withclass.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
#target_link_libraries(read_in_withclass ${ROOT_LIBRARIES})
#
#add_executable(k1k2theta12 k1k2theta12.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
#target_link_libraries(k1k2theta12 ${ROOT_LIBRARIES})
#
#add_executable(Sf1f2 Sf1f2.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
#target_link_libraries(Sf1f2 ${ROOT_LIBRARIES})
#
#add_executable(fS fS.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
#target_link_libraries(fS ${ROOT_LIBRARIES})
#
add_executable(Ternary Ternary.cc ${LIBRARIES}/cross_sections/eNCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc QEGenerator_3N.cc)
target_link_libraries(Ternary ${ROOT_LIBRARIES})


