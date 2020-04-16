add_executable(genDIS genDIS.cc ${LIBRARIES}/cross_sections/DISCrossSection.cc  ${LIBRARIES}/cross_sections/CTEQ/ct11pdf.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc DISGenerator.cc)
target_link_libraries(genDIS ${ROOT_LIBRARIES})
