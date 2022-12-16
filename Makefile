SRCDIR = src
OBJDIR = obj
BINDIR = bin
INCLUDEDIR = ./include

include local.mk

main.out : $(SRCDIR)/main.cpp $(OBJDIR)/NVT_MC_core.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)$@ $?

test_ex.out : $(SRCDIR)/test_ex.cpp $(OBJDIR)/test_mod.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $?

run_test : test_ex.out
	$(BINDIR)/test_ex.out

NVT_MC_core.o : $(OBJDIR)/my_rand.o
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/$* $? -o $(OBJDIR)/$@

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

clean :
	rm $(OBJDIR)/*.o $(BINDIR)/*.out