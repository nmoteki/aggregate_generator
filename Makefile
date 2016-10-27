TARGET = aggregate_gen
CXX = g++-6
CXXFLAGS = -O3
CPPFLAGS = -DEIGEN_NO_DEBUG -DNDEBUG
INCLUDES = -I/Users/moteki/eigen_3_2_10
$(TARGET) : aggregate_gen_main.o aggregate_gen_CCA.o aggregate_gen_SA.o aggregate_gen.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -o $@ aggregate_gen_main.o aggregate_gen_CCA.o aggregate_gen_SA.o

aggregate_gen_main.o : aggregate_gen_main.cpp aggregate_gen.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c aggregate_gen_main.cpp

aggregate_gen_CCA.o : aggregate_gen_CCA.cpp aggregate_gen.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c aggregate_gen_CCA.cpp

aggregate_gen_SA.o : aggregate_gen_SA.cpp aggregate_gen.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c aggregate_gen_SA.cpp

clean:
	rm -f *.o $(TARGET)
