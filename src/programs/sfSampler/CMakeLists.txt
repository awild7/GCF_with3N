add_executable(sfSampler sfSampler.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc)
target_link_libraries(sfSampler ${ROOT_LIBRARIES})

add_executable(sfSamplerLC sfSamplerLC.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc ${LIBRARIES}/nucleus/gcfSRC.cc)
target_link_libraries(sfSamplerLC ${ROOT_LIBRARIES})
