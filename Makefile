# Lista głównych programów wykonywalnych
PROGS = main

# Ścieżki do katalogów
CAPDBINDIR =~/CAPD/build/bin/
SRCDIR = src
INCLUDEDIR = include
OBJDIR = .obj

# Pliki źródłowe i nagłówkowe
SOURCES = $(filter-out $(SRCDIR)/main.cpp, $(wildcard $(SRCDIR)/*.cpp))
HEADERS = $(wildcard $(INCLUDEDIR)/*.h)

# Pliki obiektowe (bez main.cpp)
OBJ_FILES = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Kompilator i flagi
CAPDCXX := $(shell $(CAPDBINDIR)capd-config --variable=capd_cxx)
CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CXXFLAGS += $(CAPDFLAGS) -I$(INCLUDEDIR)

# Główna reguła
.PHONY: all
all: $(PROGS)

# Linkowanie plików wykonywalnych
$(PROGS): % : $(OBJDIR)/%.o $(OBJ_FILES)
	$(CAPDCXX) -o $@ $< $(OBJ_FILES) $(CAPDLIBS)

# Kompilacja plików źródłowych na obiekty
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CAPDCXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $(@:.o=.d) -c -o $@ $<

# Włączanie plików zależności
-include $(OBJ_FILES:.o=.d)

# Czyszczenie
.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(PROGS)