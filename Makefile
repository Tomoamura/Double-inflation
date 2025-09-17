MODEL = STOLAS
NOISE = noisemap
FUNC = functions
AVE = curvature_av

CXX := g++

CXXFLAGS := -std=c++14 -O2 -Xpreprocessor -fopenmp
CPPFLAGS := -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/boost/include

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
LDLIBS := -lfftw3_threads -lfftw3 -lm -lomp

all: $(NOISE) $(MODEL)# $(FUNC) $(AVE)

$(NOISE): $(NOISE).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(FUNC): $(FUNC).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(MODEL): $(MODEL).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(AVE): $(AVE).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


$(NOISE).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(FUNC).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(MODEL).o: src/vec_op.hpp STOLAS.hpp parameters.hpp model.hpp
$(AVE).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp


clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(NOISE)
	$(RM) $(FUNC)
	$(RM) $(MODEL)
	$(RM) $(AVE)
